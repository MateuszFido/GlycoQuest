// Copyright (c) ETH Zurich, Mateusz Fido

//! xQuest job planning from prefilter outputs.

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

use rayon::prelude::*;

use crate::cli::settings::Settings;
use crate::crosslinker::CrosslinkerProfile;
use crate::glyco::{GlycanEntry, GlycanLibrary};
use crate::prefilter::{FilteredSpectrum, GlycanPruningRow, PrefilterResult};
use crate::progress::PhaseProgress;

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
    /// Shared by the intact, -H2O, and -2H2O variants of one glycan.
    pub spectrum_keys: Arc<[SpectrumKey]>,
    /// Number of spectra assigned to this xQuest job.
    pub estimated_comparisons: u64,
}

/// Records what each generated xQuest job searched for, so hits can be annotated
/// with the glycan and attachment residue after xQuest runs.
#[derive(Debug, Clone, PartialEq)]
pub struct JobManifestEntry {
    pub job_id: String,
    pub variant: GlycanVariant,
    pub varmod_plan: VarModPlan,
    pub source_file: PathBuf,
    pub spectrum_keys: Arc<[SpectrumKey]>,
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
    #[cfg(test)]
    pub fn build(
        prefilter: &PrefilterResult,
        library: &GlycanLibrary,
        settings: &Settings,
        crosslinker: &CrosslinkerProfile,
    ) -> Result<Self, String> {
        Self::build_with_progress(prefilter, library, settings, crosslinker, None)
    }

    pub(crate) fn build_with_progress(
        prefilter: &PrefilterResult,
        library: &GlycanLibrary,
        settings: &Settings,
        _crosslinker: &CrosslinkerProfile,
        progress: Option<&PhaseProgress>,
    ) -> Result<Self, String> {
        if let Some(progress) = progress {
            progress.set_message("indexing glycans and retained spectra");
        }
        let indexes = PlanningIndexes::new(prefilter, library);
        let bucket_count = library.entries.len();
        let progress = progress.cloned();

        if let Some(progress) = &progress {
            progress.set_message(format!(
                "grouping {} spectrum-glycan candidates",
                prefilter.pruning.len()
            ));
        }

        // Each Rayon task groups into integer spectrum indexes. This avoids cloning
        // PathBufs for every pruning row and lets tasks merge without shared locks.
        let mut spectra_by_glycan = prefilter
            .pruning
            .par_chunks(4_096)
            .map(|rows| -> Result<Vec<Vec<usize>>, String> {
                let mut local = vec![Vec::new(); bucket_count];
                for row in rows {
                    let (glycan_index, spectrum_index) = indexes.resolve(row)?;
                    local[glycan_index].push(spectrum_index);
                }
                if let Some(progress) = &progress {
                    progress.inc(rows.len() as u64);
                }
                Ok(local)
            })
            .try_reduce(
                || vec![Vec::new(); bucket_count],
                |mut left, mut right| {
                    for (left_bucket, right_bucket) in left.iter_mut().zip(&mut right) {
                        left_bucket.append(right_bucket);
                    }
                    Ok(left)
                },
            )?;

        if let Some(progress) = &progress {
            progress.set_message("materializing glycan job variants");
        }

        let job_groups: Result<Vec<_>, String> = spectra_by_glycan
            .par_iter_mut()
            .enumerate()
            .map(|(glycan_index, spectrum_indexes)| {
                spectrum_indexes.sort_unstable();
                spectrum_indexes.dedup();
                build_jobs_for_glycan(
                    &library.entries[glycan_index],
                    spectrum_indexes,
                    &prefilter.filtered,
                )
            })
            .collect();

        let mut jobs = Vec::new();
        let mut total_comparisons = 0u64;
        for (mut glycan_jobs, comparisons) in job_groups? {
            total_comparisons += comparisons;
            jobs.append(&mut glycan_jobs);
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
                "spectrum-job assignments ({total_comparisons}) exceed max_total_job_spectrum_comparisons ({})",
                settings.max_total_job_spectrum_comparisons
            ));
        }

        Ok(Self {
            jobs,
            total_comparisons,
        })
    }
}

struct PlanningIndexes {
    glycan_by_name: HashMap<String, usize>,
    glycan_by_composition: HashMap<String, usize>,
    file_ids: HashMap<PathBuf, usize>,
    spectrum_by_scan: HashMap<(usize, u32), usize>,
}

impl PlanningIndexes {
    fn new(prefilter: &PrefilterResult, library: &GlycanLibrary) -> Self {
        let mut glycan_by_name = HashMap::with_capacity(library.entries.len());
        let mut glycan_by_composition = HashMap::with_capacity(library.entries.len());
        for (index, entry) in library.entries.iter().enumerate() {
            glycan_by_name.entry(entry.name.clone()).or_insert(index);
            glycan_by_composition
                .entry(entry.composition.clone())
                .or_insert(index);
        }

        let mut file_ids = HashMap::new();
        let mut spectrum_by_scan = HashMap::with_capacity(prefilter.filtered.len());
        for (spectrum_index, spectrum) in prefilter.filtered.iter().enumerate() {
            let next_file_id = file_ids.len();
            let file_id = *file_ids
                .entry(spectrum.source_file.clone())
                .or_insert(next_file_id);
            spectrum_by_scan
                .entry((file_id, spectrum.scan_number))
                .or_insert(spectrum_index);
        }

        Self {
            glycan_by_name,
            glycan_by_composition,
            file_ids,
            spectrum_by_scan,
        }
    }

    fn resolve(&self, row: &GlycanPruningRow) -> Result<(usize, usize), String> {
        // The old linear search returned the earliest entry matching either field.
        // Preserve that behavior when a name and composition resolve differently.
        let by_name = self.glycan_by_name.get(&row.glycan_name).copied();
        let by_composition = self.glycan_by_composition.get(&row.composition).copied();
        let glycan_index = match (by_name, by_composition) {
            (Some(name), Some(composition)) => name.min(composition),
            (Some(index), None) | (None, Some(index)) => index,
            (None, None) => {
                return Err(format!("pruned glycan not in library: {}", row.glycan_name));
            }
        };

        let file_id = self.file_ids.get(&row.source_file).ok_or_else(|| {
            format!(
                "pruned spectrum file not in filtered spectra: {}",
                row.source_file.display()
            )
        })?;
        let spectrum_index = self
            .spectrum_by_scan
            .get(&(*file_id, row.scan_number))
            .copied()
            .ok_or_else(|| {
                format!(
                    "pruned spectrum missing from filtered spectra: {} scan {}",
                    row.source_file.display(),
                    row.scan_number
                )
            })?;
        Ok((glycan_index, spectrum_index))
    }
}

fn glycan_variants(entry: &GlycanEntry) -> Result<[GlycanVariant; 3], String> {
    let residue_targets = residue_targets_to_chars(&entry.residue_targets)?;

    let make = |mass: f64, loss_label: &str| GlycanVariant {
        glycan_name: entry.name.clone(),
        composition: entry.composition.clone(),
        mass,
        loss_label: loss_label.to_string(),
        residue_targets: residue_targets.clone(),
    };

    Ok([
        make(entry.monoisotopic_mass, ""),
        make(entry.monoisotopic_mass - WATER_LOSS_DA, "-H2O"),
        make(entry.monoisotopic_mass - 2.0 * WATER_LOSS_DA, "-2H2O"),
    ])
}

fn build_jobs_for_glycan(
    entry: &GlycanEntry,
    spectrum_indexes: &[usize],
    filtered: &[FilteredSpectrum],
) -> Result<(Vec<PlannedJob>, u64), String> {
    if spectrum_indexes.is_empty() {
        return Ok((Vec::new(), 0));
    }

    let spectrum_keys: Arc<[SpectrumKey]> = spectrum_indexes
        .iter()
        .map(|&index| SpectrumKey {
            source_file: filtered[index].source_file.clone(),
            scan_number: filtered[index].scan_number,
        })
        .collect::<Vec<_>>()
        .into();
    // xQuest's native progress is one unit per searched spectrum. Counting
    // precursor charges here inflated the progress total even though xQuest
    // does not emit a separate progress event for each charge hypothesis.
    let comparisons = spectrum_indexes.len() as u64;
    let jobs = glycan_variants(entry)?
        .into_iter()
        .map(|variant| PlannedJob {
            job_id: sanitize_job_id(&variant),
            variant,
            spectrum_keys: Arc::clone(&spectrum_keys),
            estimated_comparisons: comparisons,
        })
        .collect();
    Ok((jobs, comparisons * 3))
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crosslinker::CrosslinkerProfile;
    use crate::glyco::{GlycanEntry, load_glycan_database};
    use crate::prefilter::{PrefilterStats, run_prefilter};
    use crate::{ProgressMode, progress::ProgressReporter};
    use std::path::PathBuf;
    use std::time::{Duration, Instant};

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
        assert_eq!(
            plan.total_comparisons,
            plan.jobs
                .iter()
                .map(|job| job.estimated_comparisons)
                .sum::<u64>()
        );
    }

    #[test]
    fn large_candidate_matrix_plans_without_quadratic_spectrum_lookups() {
        const SPECTRA: usize = 1_200;
        const GLYCANS: usize = 80;
        let source = PathBuf::from("synthetic.mzXML");
        let filtered: Vec<_> = (1..=SPECTRA)
            .map(|scan| FilteredSpectrum {
                source_file: source.clone(),
                scan_number: scan as u32,
                retention_time_min: scan as f64 / 60.0,
                precursor_mz: 800.0,
                precursor_charge: Some(2),
                matched_families: vec!["HexNAc".into()],
                matched_ions: vec![],
            })
            .collect();
        let library = GlycanLibrary {
            database_id: "synthetic".into(),
            entries: (0..GLYCANS)
                .map(|index| GlycanEntry {
                    name: format!("glycan-{index}"),
                    composition: format!("HexNAc(2)Hex({})", index + 1),
                    monoisotopic_mass: 1_000.0 + index as f64,
                    diagnostic_ions: vec![],
                    residue_targets: vec!["N".into()],
                })
                .collect(),
        };
        let pruning = filtered
            .iter()
            .flat_map(|spectrum| {
                library.entries.iter().map(|glycan| GlycanPruningRow {
                    source_file: spectrum.source_file.clone(),
                    scan_number: spectrum.scan_number,
                    glycan_name: glycan.name.clone(),
                    composition: glycan.composition.clone(),
                    matched_families: vec!["HexNAc".into()],
                })
            })
            .collect();
        let prefilter = PrefilterResult {
            filtered,
            isotope_pairs: vec![],
            rejected: vec![],
            pruning,
            stats: PrefilterStats {
                scans_total: SPECTRA,
                diagnostic_positive: SPECTRA,
                filtered_scans: SPECTRA,
                ..PrefilterStats::default()
            },
        };
        let settings = Settings::defaults();
        let crosslinker = CrosslinkerProfile::resolve(&settings, Some("dmtmm")).unwrap();

        let progress = ProgressReporter::new(ProgressMode::Never).determinate(
            2,
            4,
            "Preparing xQuest jobs",
            prefilter.pruning.len() as u64,
        );
        let started = Instant::now();
        let plan = JobPlan::build_with_progress(
            &prefilter,
            &library,
            &settings,
            &crosslinker,
            Some(&progress),
        )
        .unwrap();
        let elapsed = started.elapsed();

        assert_eq!(plan.jobs.len(), GLYCANS * 3);
        assert!(
            plan.jobs
                .iter()
                .all(|job| job.spectrum_keys.len() == SPECTRA)
        );
        assert_eq!(plan.total_comparisons, (GLYCANS * 3 * SPECTRA) as u64);
        assert!(
            elapsed < Duration::from_secs(3),
            "planning {SPECTRA} spectra × {GLYCANS} glycans took {elapsed:?}"
        );
        assert_eq!(
            progress.test_snapshot(),
            (
                (SPECTRA * GLYCANS) as u64,
                "materializing glycan job variants".into()
            )
        );
    }
}
