//! xQuest job planning from prefilter outputs.

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use crate::cli::settings::Settings;
use crate::crosslinker::CrosslinkerProfile;
use crate::glyco::GlycanLibrary;
use crate::prefilter::{FilteredSpectrum, GlycanPruningRow, IsotopePair, PrefilterResult};

const WATER_LOSS_DA: f64 = 18.010565;

#[derive(Debug, Clone, PartialEq)]
pub struct GlycanVariant {
    pub glycan_name: String,
    pub composition: String,
    pub mass: f64,
    pub loss_label: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PlannedJob {
    pub job_id: String,
    pub variant: GlycanVariant,
    pub spectrum_keys: Vec<SpectrumKey>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct SpectrumKey {
    pub source_file: PathBuf,
    pub scan_number: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct JobPlan {
    pub jobs: Vec<PlannedJob>,
    pub total_comparisons: u64,
}

impl JobPlan {
    pub fn build(
        prefilter: &PrefilterResult,
        library: &GlycanLibrary,
        settings: &Settings,
        _crosslinker: &CrosslinkerProfile,
    ) -> Result<Self, String> {
        let mut variant_spectra: HashMap<String, HashSet<SpectrumKey>> = HashMap::new();
        let mut variant_meta: HashMap<String, GlycanVariant> = HashMap::new();

        for row in &prefilter.pruning {
            for variant in glycan_variants_for_row(row, library)? {
                let key = variant_key(&variant);
                variant_spectra
                    .entry(key.clone())
                    .or_default()
                    .insert(SpectrumKey {
                        source_file: row.source_file.clone(),
                        scan_number: row.scan_number,
                    });
                variant_meta.entry(key).or_insert(variant);
            }
        }

        let mut jobs = Vec::new();
        let mut total_comparisons = 0u64;

        for (key, spectra) in variant_spectra {
            let variant = variant_meta
                .remove(&key)
                .expect("variant metadata present");
            let spectrum_keys: Vec<_> = spectra.into_iter().collect();
            let comparisons = estimate_comparisons(&spectrum_keys, prefilter);
            total_comparisons += comparisons;

            jobs.push(PlannedJob {
                job_id: sanitize_job_id(&variant),
                variant,
                spectrum_keys,
            });
        }

        jobs.sort_by(|a, b| a.job_id.cmp(&b.job_id));

        if settings.max_jobs > 0 && jobs.len() > settings.max_jobs as usize {
            return Err(format!(
                "planned xQuest jobs ({}) exceed max_jobs ({})",
                jobs.len(),
                settings.max_jobs
            ));
        }
        if settings.max_total_job_spectrum_comparisons > 0
            && total_comparisons > settings.max_total_job_spectrum_comparisons
        {
            return Err(format!(
                "estimated spectrum comparisons ({total_comparisons}) exceed max_total_job_spectrum_comparisons ({})",
                settings.max_total_job_spectrum_comparisons
            ));
        }

        Ok(Self {
            jobs,
            total_comparisons,
        })
    }
}

fn glycan_variants_for_row(
    row: &GlycanPruningRow,
    library: &GlycanLibrary,
) -> Result<Vec<GlycanVariant>, String> {
    let entry = library
        .entries
        .iter()
        .find(|e| e.name == row.glycan_name || e.composition == row.composition)
        .ok_or_else(|| format!("pruned glycan not in library: {}", row.glycan_name))?;

    Ok(vec![
        GlycanVariant {
            glycan_name: entry.name.clone(),
            composition: entry.composition.clone(),
            mass: entry.monoisotopic_mass,
            loss_label: String::new(),
        },
        GlycanVariant {
            glycan_name: entry.name.clone(),
            composition: entry.composition.clone(),
            mass: entry.monoisotopic_mass - WATER_LOSS_DA,
            loss_label: "-H2O".into(),
        },
    ])
}

fn variant_key(variant: &GlycanVariant) -> String {
    format!("{}:{}", variant.glycan_name, variant.loss_label)
}

fn sanitize_job_id(variant: &GlycanVariant) -> String {
    let base = variant
        .composition
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect::<String>();
    if variant.loss_label.is_empty() {
        base
    } else {
        format!("{base}_{}", variant.loss_label.trim_start_matches('-'))
    }
}

fn estimate_comparisons(keys: &[SpectrumKey], prefilter: &PrefilterResult) -> u64 {
    let mut count = 0u64;
    for key in keys {
        if let Some(spec) = prefilter
            .filtered
            .iter()
            .find(|f| f.source_file == key.source_file && f.scan_number == key.scan_number)
        {
            let charge = spec.precursor_charge.unwrap_or(2) as u64;
            count += charge.max(1);
        }
    }
    count.max(1)
}

pub fn filtered_for_key<'a>(
    prefilter: &'a PrefilterResult,
    key: &SpectrumKey,
) -> Option<&'a FilteredSpectrum> {
    prefilter.filtered.iter().find(|row| {
        row.source_file == key.source_file && row.scan_number == key.scan_number
    })
}

pub fn isotope_pair_for_scan<'a>(
    pairs: &'a [IsotopePair],
    file: &PathBuf,
    scan: u32,
) -> Option<&'a IsotopePair> {
    pairs.iter().find(|pair| {
        (pair.light_file == *file && pair.light_scan == scan)
            || (pair.heavy_file == *file && pair.heavy_scan == scan)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crosslinker::CrosslinkerProfile;
    use crate::glyco::load_glycan_database;
    use crate::prefilter::run_prefilter;
    use std::path::PathBuf;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/mzxml")
            .join(name)
    }

    #[test]
    fn builds_jobs_for_dss_fixture() {
        let library = load_glycan_database("nglyc309").unwrap();
        let settings = Settings::defaults();
        let crosslinker = CrosslinkerProfile::resolve(&settings, Some("dss")).unwrap();
        let files = vec![fixture("dss_pair.mzXML")];
        let prefilter = run_prefilter(&files, &library, &settings, &crosslinker).unwrap();
        let plan = JobPlan::build(&prefilter, &library, &settings, &crosslinker).unwrap();
        assert!(!plan.jobs.is_empty());
        assert!(plan.total_comparisons > 0);
    }
}
