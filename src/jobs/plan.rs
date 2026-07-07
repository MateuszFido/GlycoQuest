//! xQuest job planning from prefilter outputs.

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use crate::cli::settings::Settings;
use crate::crosslinker::CrosslinkerProfile;
use crate::glyco::GlycanLibrary;
use crate::prefilter::{FilteredSpectrum, GlycanPruningRow, PrefilterResult};

const WATER_LOSS_DA: f64 = 18.010565;
const OXIDATION_DA: f64 = 15.994915;

/// xQuest exposes exactly four variable-modification pseudo-residues.
const PSEUDO_RESIDUES: [char; 4] = ['X', 'U', 'B', 'J'];

#[derive(Debug, Clone, PartialEq)]
pub struct GlycanVariant {
    pub glycan_name: String,
    pub composition: String,
    pub mass: f64,
    pub loss_label: String,
    /// Source residues a glycan of this variant may attach to (e.g. `N`, or `S`/`T`).
    pub residue_targets: Vec<char>,
}

/// What a variable-modification pseudo-residue represents.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VarModKind {
    Glycan,
    Oxidation,
}

/// One xQuest variable modification, mapped to a pseudo-residue (`X`/`U`/`B`/`J`).
#[derive(Debug, Clone, PartialEq)]
pub struct VarModEntry {
    pub pseudo: char,
    pub source_residue: char,
    pub mass: f64,
    pub kind: VarModKind,
}

/// Ordered variable-modification layout for a single xQuest job.
///
/// `entries[i]` maps to xQuest pseudo-residue `PSEUDO_RESIDUES[i]`.
#[derive(Debug, Clone, PartialEq)]
pub struct VarModPlan {
    pub entries: Vec<VarModEntry>,
}

impl VarModPlan {
    /// The comma-separated `variable_mod` value written to `xquest.def`.
    pub fn variable_mod_value(&self) -> String {
        self.entries
            .iter()
            .map(|entry| format!("{},{:.6}", entry.source_residue, entry.mass))
            .collect::<Vec<_>>()
            .join(",")
    }

    /// Number of variable modifications (`nvariable_mod`).
    pub fn nvariable_mod(&self) -> usize {
        self.entries.len()
    }

    /// The set of pseudo-residues that represent a glycan attachment.
    pub fn glycan_pseudos(&self) -> Vec<char> {
        self.entries
            .iter()
            .filter(|entry| entry.kind == VarModKind::Glycan)
            .map(|entry| entry.pseudo)
            .collect()
    }

    /// Map a pseudo-residue back to the entry it represents.
    pub fn entry_for_pseudo(&self, pseudo: char) -> Option<&VarModEntry> {
        self.entries.iter().find(|entry| entry.pseudo == pseudo)
    }
}

/// Build the variable-modification layout for a glycan variant.
///
/// One glycan entry per configured attachment residue, followed by an optional
/// oxidation entry. Fails if more than four entries are required, because xQuest
/// only exposes four pseudo-residues.
pub fn build_varmod_plan(
    variant: &GlycanVariant,
    settings: &Settings,
) -> Result<VarModPlan, String> {
    let mut entries = Vec::new();

    for &residue in &variant.residue_targets {
        entries.push(VarModEntry {
            pseudo: '?',
            source_residue: residue,
            mass: variant.mass,
            kind: VarModKind::Glycan,
        });
    }

    if entries.is_empty() {
        return Err(format!(
            "glycan {} has no attachment residue targets",
            variant.glycan_name
        ));
    }

    if settings.variable_oxidation {
        entries.push(VarModEntry {
            pseudo: '?',
            source_residue: 'M',
            mass: OXIDATION_DA,
            kind: VarModKind::Oxidation,
        });
    }

    if entries.len() > PSEUDO_RESIDUES.len() {
        return Err(format!(
            "glycan {} requires {} variable modifications but xQuest supports at most {}",
            variant.glycan_name,
            entries.len(),
            PSEUDO_RESIDUES.len()
        ));
    }

    for (index, entry) in entries.iter_mut().enumerate() {
        entry.pseudo = PSEUDO_RESIDUES[index];
    }

    Ok(VarModPlan { entries })
}

#[derive(Debug, Clone, PartialEq)]
pub struct PlannedJob {
    pub job_id: String,
    pub variant: GlycanVariant,
    pub spectrum_keys: Vec<SpectrumKey>,
}

/// Records what each generated xQuest job searched for, so hits can be annotated
/// with the glycan and attachment residue after xQuest runs.
#[derive(Debug, Clone, PartialEq)]
pub struct JobManifestEntry {
    pub job_id: String,
    pub variant: GlycanVariant,
    pub varmod_plan: VarModPlan,
    pub source_file: PathBuf,
    pub spectrum_keys: Vec<SpectrumKey>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct JobManifest {
    pub entries: Vec<JobManifestEntry>,
}

impl JobManifest {
    pub fn by_job_id(&self, job_id: &str) -> Option<&JobManifestEntry> {
        self.entries.iter().find(|entry| entry.job_id == job_id)
    }
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

    let residue_targets = residue_targets_to_chars(&entry.residue_targets)?;

    let make = |mass: f64, loss_label: &str| GlycanVariant {
        glycan_name: entry.name.clone(),
        composition: entry.composition.clone(),
        mass,
        loss_label: loss_label.to_string(),
        residue_targets: residue_targets.clone(),
    };

    Ok(vec![
        make(entry.monoisotopic_mass, ""),
        make(entry.monoisotopic_mass - WATER_LOSS_DA, "-H2O"),
        make(entry.monoisotopic_mass - 2.0 * WATER_LOSS_DA, "-2H2O"),
    ])
}

fn residue_targets_to_chars(targets: &[String]) -> Result<Vec<char>, String> {
    let mut chars = Vec::with_capacity(targets.len());
    for target in targets {
        let trimmed = target.trim();
        let mut iter = trimmed.chars();
        match (iter.next(), iter.next()) {
            (Some(ch), None) => chars.push(ch),
            _ => {
                return Err(format!(
                    "glycan residue target must be a single residue letter, got '{trimmed}'"
                ));
            }
        }
    }
    Ok(chars)
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
