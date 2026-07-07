//! GlycoQuest post-filters and glycan annotation for xQuest hits.

use std::path::PathBuf;

use crate::cli::settings::Settings;
use crate::crosslinker::CrosslinkerProfile;
use crate::jobs::{JobManifest, VarModPlan};
use crate::prefilter::{FilteredSpectrum, PrefilterResult};
use crate::results::extract::XQuestHit;

/// A hit annotated with its glycan, originating spectrum, and post-filter outcome.
#[derive(Debug, Clone, PartialEq)]
pub struct AnnotatedHit {
    pub hit: XQuestHit,
    pub job_id: String,
    pub source_file: Option<PathBuf>,
    pub scan: Option<u32>,
    pub glycan_name: Option<String>,
    pub glycan_composition: Option<String>,
    pub glycan_mass: Option<f64>,
    pub loss_label: Option<String>,
    pub glyco_residue: Option<char>,
    pub glyco_peptide: Option<u8>,
    pub n_glycan_pseudo: usize,
    pub matched_families: Vec<String>,
    pub matched_ion_count: usize,
    pub sequon_present: Option<bool>,
    pub charge_plausible: bool,
    pub hard_status: HardStatus,
    pub soft_score: f64,
    pub postfilter_status: PostfilterStatus,
}

/// Result of the hard (pass/fail) post-filter requirements.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HardStatus {
    Pass,
    FailNoXlink,
    FailGlycanCount,
    FailNoDiagnostic,
    FailPrecursorError,
    FailScore,
}

impl HardStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Pass => "pass",
            Self::FailNoXlink => "fail_no_xlink",
            Self::FailGlycanCount => "fail_glycan_count",
            Self::FailNoDiagnostic => "fail_no_diagnostic",
            Self::FailPrecursorError => "fail_precursor_error",
            Self::FailScore => "fail_score",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PostfilterStatus {
    Pass,
    Fail,
}

impl PostfilterStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Pass => "pass",
            Self::Fail => "fail",
        }
    }
}

/// Annotate and post-filter xQuest hits.
///
/// `hits` are `(job_id, hit)` pairs; the manifest supplies the glycan searched by
/// each job, and the prefilter supplies diagnostic-ion evidence per scan. When the
/// manifest is absent (e.g. `--resume`), annotation falls back to the job id and
/// prefilter-dependent checks are skipped.
pub fn apply_postfilters(
    hits: Vec<(String, XQuestHit)>,
    settings: &Settings,
    _crosslinker: &CrosslinkerProfile,
    prefilter: &PrefilterResult,
    manifest: Option<&JobManifest>,
) -> Vec<AnnotatedHit> {
    hits.into_iter()
        .map(|(job_id, hit)| annotate_hit(job_id, hit, settings, prefilter, manifest))
        .collect()
}

fn annotate_hit(
    job_id: String,
    hit: XQuestHit,
    settings: &Settings,
    prefilter: &PrefilterResult,
    manifest: Option<&JobManifest>,
) -> AnnotatedHit {
    let entry = manifest.and_then(|m| m.by_job_id(&job_id));

    let (glycan_name, glycan_composition, glycan_mass, loss_label) = match entry {
        Some(entry) => (
            Some(entry.variant.glycan_name.clone()),
            Some(entry.variant.composition.clone()),
            Some(entry.variant.mass),
            Some(loss_label_or_none(&entry.variant.loss_label)),
        ),
        None => {
            let (name, loss) = glycan_from_job_id(&job_id);
            (name.clone(), name, None, loss)
        }
    };

    let plan = entry.map(|e| &e.varmod_plan);
    let n_glycan_pseudo = plan
        .map(|plan| count_glycan_pseudos(&hit.seq1, &hit.seq2, plan))
        .unwrap_or_else(|| fallback_glycan_pseudo_count(&hit.seq1, &hit.seq2));

    let site = plan.and_then(|plan| glyco_site(&hit.seq1, &hit.seq2, plan));
    let glyco_peptide = site.map(|(peptide, _, _)| peptide);
    let glyco_residue = site.map(|(_, _, residue)| residue);

    let scan = parse_scan(&hit.spectrum_id);
    let source_file = entry.map(|e| e.source_file.clone());
    let spectrum = lookup_spectrum(prefilter, source_file.as_deref(), scan);

    let matched_families = spectrum
        .map(|spec| spec.matched_families.clone())
        .unwrap_or_default();
    let matched_ion_count = spectrum.map(|spec| spec.matched_ions.len()).unwrap_or(0);

    let sequon_present = site.and_then(|(peptide, pos, residue)| {
        if residue == 'N' {
            let seq = if peptide == 1 { &hit.seq1 } else { &hit.seq2 };
            Some(has_sequon(seq, pos, plan))
        } else {
            None
        }
    });

    let charge_plausible = charge_plausible(hit.charge);

    let hard_status = hard_status(
        &hit,
        settings,
        n_glycan_pseudo,
        spectrum.is_some(),
        manifest.is_some(),
    );

    let soft_score = soft_score(
        &hit,
        sequon_present,
        charge_plausible,
        matched_ion_count,
    );

    let postfilter_status = if hard_status == HardStatus::Pass {
        PostfilterStatus::Pass
    } else {
        PostfilterStatus::Fail
    };

    AnnotatedHit {
        hit,
        job_id,
        source_file,
        scan,
        glycan_name,
        glycan_composition,
        glycan_mass,
        loss_label,
        glyco_residue,
        glyco_peptide,
        n_glycan_pseudo,
        matched_families,
        matched_ion_count,
        sequon_present,
        charge_plausible,
        hard_status,
        soft_score,
        postfilter_status,
    }
}

fn hard_status(
    hit: &XQuestHit,
    settings: &Settings,
    n_glycan_pseudo: usize,
    diagnostic_positive: bool,
    have_manifest: bool,
) -> HardStatus {
    if hit.xlink_position.trim().is_empty() && hit.topology.trim().is_empty() {
        return HardStatus::FailNoXlink;
    }
    // V1 class: peptide-glycopeptide crosslink requires exactly one glycan.
    if n_glycan_pseudo != 1 {
        return HardStatus::FailGlycanCount;
    }
    // Diagnostic-ion evidence in the originating spectrum is a hard requirement,
    // but only checkable when prefilter state is available (a normal run).
    if have_manifest && !diagnostic_positive {
        return HardStatus::FailNoDiagnostic;
    }
    if hit.precursor_error_ppm.abs() > settings.max_precursor_error_ppm {
        return HardStatus::FailPrecursorError;
    }
    if hit.score < settings.min_score {
        return HardStatus::FailScore;
    }
    HardStatus::Pass
}

fn soft_score(
    hit: &XQuestHit,
    sequon_present: Option<bool>,
    charge_plausible: bool,
    matched_ion_count: usize,
) -> f64 {
    let mut score = hit.score;
    if sequon_present == Some(true) {
        score += 1.0;
    }
    if charge_plausible {
        score += 0.5;
    }
    score += (matched_ion_count as f64).min(10.0) * 0.1;
    let mass_penalty = (hit.precursor_error_ppm.abs() / 10.0).min(1.0);
    score -= mass_penalty;
    score
}

/// Occurrences of any glycan pseudo-residue across both peptide sequences.
pub fn count_glycan_pseudos(seq1: &str, seq2: &str, plan: &VarModPlan) -> usize {
    let glycan_pseudos = plan.glycan_pseudos();
    [seq1, seq2]
        .iter()
        .flat_map(|seq| seq.chars())
        .filter(|ch| glycan_pseudos.contains(ch))
        .count()
}

/// The peptide (1 or 2), 0-based position, and source residue of the single glycan site.
pub fn glyco_site(seq1: &str, seq2: &str, plan: &VarModPlan) -> Option<(u8, usize, char)> {
    let glycan_pseudos = plan.glycan_pseudos();
    for (peptide, seq) in [(1u8, seq1), (2u8, seq2)] {
        for (pos, ch) in seq.chars().enumerate() {
            if glycan_pseudos.contains(&ch) {
                let residue = plan
                    .entry_for_pseudo(ch)
                    .map(|entry| entry.source_residue)
                    .unwrap_or(ch);
                return Some((peptide, pos, residue));
            }
        }
    }
    None
}

/// Is there an N-glycosylation sequon (N-X-S/T, X != P) at `pos` in the peptide?
///
/// `pos` refers to the glycosylated N. Pseudo-residues are resolved to their
/// source residue before evaluating the sequon.
pub fn has_sequon(seq: &str, pos: usize, plan: Option<&VarModPlan>) -> bool {
    let resolved: Vec<char> = seq
        .chars()
        .map(|ch| {
            plan.and_then(|plan| plan.entry_for_pseudo(ch))
                .map(|entry| entry.source_residue)
                .unwrap_or(ch)
        })
        .collect();

    let x = match resolved.get(pos + 1) {
        Some(&c) => c,
        None => return false,
    };
    if x == 'P' {
        return false;
    }
    matches!(resolved.get(pos + 2), Some('S') | Some('T'))
}

/// Precursor charges plausible for (glyco)peptide-peptide crosslinks.
pub fn charge_plausible(charge: u8) -> bool {
    (2..=7).contains(&charge)
}

/// Recover the light scan number from an xQuest spectrum id.
///
/// compare_peaks3.pl names each spectrum `basename(col4)_basename(col5)`, and
/// GlycoQuest writes those columns as `{scan}.{stem}` (light first), so the id
/// looks like `2173.sample.c_1896.sample.c`. Taking the leading run of digits
/// yields the light scan regardless of what the source filename contains.
fn parse_scan(spectrum_id: &str) -> Option<u32> {
    let digits: String = spectrum_id
        .trim()
        .chars()
        .take_while(|c| c.is_ascii_digit())
        .collect();
    digits.parse().ok()
}

fn lookup_spectrum<'a>(
    prefilter: &'a PrefilterResult,
    source_file: Option<&std::path::Path>,
    scan: Option<u32>,
) -> Option<&'a FilteredSpectrum> {
    let scan = scan?;
    prefilter.filtered.iter().find(|spec| {
        spec.scan_number == scan
            && source_file
                .map(|file| spec.source_file == file)
                .unwrap_or(true)
    })
}

fn loss_label_or_none(label: &str) -> String {
    if label.trim().is_empty() {
        "none".to_string()
    } else {
        label.to_string()
    }
}

/// Reconstruct a glycan label from a job id, e.g. `HexNAc_1_Hex_5__H2O`.
fn glycan_from_job_id(job_id: &str) -> (Option<String>, Option<String>) {
    if job_id.is_empty() {
        return (None, None);
    }
    for suffix in ["_2H2O", "_H2O"] {
        if let Some(stem) = job_id.strip_suffix(suffix) {
            let loss = format!("-{}", suffix.trim_start_matches('_'));
            return (Some(stem.trim_end_matches('_').to_string()), Some(loss));
        }
    }
    (Some(job_id.trim_end_matches('_').to_string()), Some("none".to_string()))
}

/// Without a manifest we cannot know the pseudo-residue set, so use xQuest's
/// reserved variable-mod letters as a best-effort glycan indicator.
fn fallback_glycan_pseudo_count(seq1: &str, seq2: &str) -> usize {
    const PSEUDO: [char; 4] = ['X', 'U', 'B', 'J'];
    [seq1, seq2]
        .iter()
        .flat_map(|seq| seq.chars())
        .filter(|ch| PSEUDO.contains(ch))
        .count()
}

pub fn write_annotated_csv(path: &std::path::Path, rows: &[AnnotatedHit]) -> Result<(), String> {
    let mut lines = vec![
        [
            "source_file",
            "scan",
            "glycan_name",
            "glycan_composition",
            "glycan_mass",
            "loss_label",
            "glyco_residue",
            "glyco_peptide",
            "n_glycan_pseudo",
            "sequon_present",
            "charge",
            "charge_plausible",
            "matched_families",
            "matched_ion_count",
            "seq1",
            "seq2",
            "prot1",
            "prot2",
            "topology",
            "precursor_mz",
            "mr",
            "precursor_error_ppm",
            "xlink_position",
            "score",
            "hard_status",
            "soft_score",
            "postfilter_status",
        ]
        .join("\t"),
    ];

    for row in rows {
        let hit = &row.hit;
        lines.push(
            [
                opt_path(&row.source_file),
                opt_u32(row.scan),
                row.glycan_name.clone().unwrap_or_default(),
                row.glycan_composition.clone().unwrap_or_default(),
                opt_f64(row.glycan_mass),
                row.loss_label.clone().unwrap_or_default(),
                opt_char(row.glyco_residue),
                opt_u8(row.glyco_peptide),
                row.n_glycan_pseudo.to_string(),
                opt_bool(row.sequon_present),
                hit.charge.to_string(),
                row.charge_plausible.to_string(),
                row.matched_families.join(";"),
                row.matched_ion_count.to_string(),
                hit.seq1.clone(),
                hit.seq2.clone(),
                hit.prot1.clone(),
                hit.prot2.clone(),
                hit.topology.clone(),
                format!("{}", hit.precursor_mz),
                format!("{}", hit.mr),
                format!("{}", hit.precursor_error_ppm),
                hit.xlink_position.clone(),
                format!("{}", hit.score),
                row.hard_status.as_str().to_string(),
                format!("{:.3}", row.soft_score),
                row.postfilter_status.as_str().to_string(),
            ]
            .join("\t"),
        );
    }
    std::fs::write(path, lines.join("\n") + "\n").map_err(|err| err.to_string())
}

fn opt_path(value: &Option<PathBuf>) -> String {
    value
        .as_ref()
        .map(|p| p.display().to_string())
        .unwrap_or_default()
}

fn opt_u32(value: Option<u32>) -> String {
    value.map(|v| v.to_string()).unwrap_or_default()
}

fn opt_u8(value: Option<u8>) -> String {
    value.map(|v| v.to_string()).unwrap_or_default()
}

fn opt_char(value: Option<char>) -> String {
    value.map(|v| v.to_string()).unwrap_or_default()
}

fn opt_f64(value: Option<f64>) -> String {
    value.map(|v| format!("{v:.6}")).unwrap_or_default()
}

fn opt_bool(value: Option<bool>) -> String {
    value.map(|v| v.to_string()).unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crosslinker::CrosslinkerProfile;
    use crate::jobs::{build_varmod_plan, GlycanVariant, JobManifestEntry, SpectrumKey};
    use crate::prefilter::{FilteredSpectrum, MatchedIon};

    fn n_glycan_plan() -> VarModPlan {
        let variant = GlycanVariant {
            glycan_name: "HexNAc(1)".into(),
            composition: "HexNAc(1)".into(),
            mass: 203.079373,
            loss_label: String::new(),
            residue_targets: vec!['N'],
        };
        build_varmod_plan(&variant, &Settings::defaults()).unwrap()
    }

    fn manifest_with(job_id: &str, source_file: PathBuf) -> JobManifest {
        let variant = GlycanVariant {
            glycan_name: "HexNAc(1)".into(),
            composition: "HexNAc(1)".into(),
            mass: 203.079373,
            loss_label: String::new(),
            residue_targets: vec!['N'],
        };
        let varmod_plan = build_varmod_plan(&variant, &Settings::defaults()).unwrap();
        JobManifest {
            entries: vec![JobManifestEntry {
                job_id: job_id.to_string(),
                variant,
                varmod_plan,
                source_file: source_file.clone(),
                spectrum_keys: vec![SpectrumKey {
                    source_file,
                    scan_number: 42,
                }],
            }],
        }
    }

    fn prefilter_with_scan(source_file: PathBuf, scan: u32) -> PrefilterResult {
        PrefilterResult {
            filtered: vec![FilteredSpectrum {
                source_file,
                scan_number: scan,
                retention_time_min: 20.0,
                precursor_mz: 800.0,
                precursor_charge: Some(4),
                matched_families: vec!["HexNAc".into()],
                matched_ions: vec![MatchedIon {
                    family: "HexNAc".into(),
                    expected_mz: 204.0866,
                    observed_mz: 204.0867,
                    loss_label: String::new(),
                }],
            }],
            isotope_pairs: vec![],
            rejected: vec![],
            pruning: vec![],
            stats: Default::default(),
        }
    }

    #[test]
    fn multiple_glycan_pseudo_residues_fail_glycan_count() {
        let settings = Settings::defaults();
        let crosslinker = CrosslinkerProfile::resolve(&settings, Some("dss")).unwrap();
        let source = PathBuf::from("run.mzXML");
        let manifest = manifest_with("HexNAc_1_", source.clone());
        let prefilter = prefilter_with_scan(source, 42);
        let hits = vec![(
            "HexNAc_1_".to_string(),
            XQuestHit {
                seq1: "AXCXDE".into(),
                seq2: "PEPXIDE".into(),
                xlink_position: "2-1".into(),
                charge: 4,
                spectrum_id: "42".into(),
                ..Default::default()
            },
        )];
        let annotated = apply_postfilters(hits, &settings, &crosslinker, &prefilter, Some(&manifest));
        assert_eq!(annotated[0].hard_status, HardStatus::FailGlycanCount);
        assert_eq!(annotated[0].postfilter_status, PostfilterStatus::Fail);
    }

    #[test]
    fn single_glycan_hit_annotated_and_passes() {
        let settings = Settings::defaults();
        let crosslinker = CrosslinkerProfile::resolve(&settings, Some("dss")).unwrap();
        let source = PathBuf::from("run.mzXML");
        let manifest = manifest_with("HexNAc_1_", source.clone());
        let prefilter = prefilter_with_scan(source.clone(), 42);
        let hits = vec![(
            "HexNAc_1_".to_string(),
            XQuestHit {
                seq1: "AXCDE".into(),
                seq2: "PEPKIDE".into(),
                xlink_position: "3-1".into(),
                charge: 4,
                spectrum_id: "42".into(),
                score: 5.0,
                ..Default::default()
            },
        )];
        let annotated = apply_postfilters(hits, &settings, &crosslinker, &prefilter, Some(&manifest));
        let row = &annotated[0];
        assert_eq!(row.hard_status, HardStatus::Pass);
        assert_eq!(row.glycan_composition.as_deref(), Some("HexNAc(1)"));
        assert_eq!(row.glyco_residue, Some('N'));
        assert_eq!(row.glyco_peptide, Some(1));
        assert_eq!(row.scan, Some(42));
        assert_eq!(row.n_glycan_pseudo, 1);
    }

    #[test]
    fn non_glycosylated_pair_fails_glycan_count() {
        let settings = Settings::defaults();
        let crosslinker = CrosslinkerProfile::resolve(&settings, Some("dss")).unwrap();
        let source = PathBuf::from("run.mzXML");
        let manifest = manifest_with("HexNAc_1_", source.clone());
        let prefilter = prefilter_with_scan(source, 42);
        let hits = vec![(
            "HexNAc_1_".to_string(),
            XQuestHit {
                seq1: "DDSEKGK".into(),
                seq2: "KDAGGR".into(),
                xlink_position: "5-1".into(),
                charge: 3,
                spectrum_id: "42".into(),
                ..Default::default()
            },
        )];
        let annotated = apply_postfilters(hits, &settings, &crosslinker, &prefilter, Some(&manifest));
        assert_eq!(annotated[0].hard_status, HardStatus::FailGlycanCount);
    }

    #[test]
    fn sequon_detection_resolves_pseudo_residues() {
        let plan = n_glycan_plan();
        // X is glycosylated N; sequon N-K-T at positions 0,1,2.
        assert!(has_sequon("XKT", 0, Some(&plan)));
        // N-P-T is not a sequon.
        assert!(!has_sequon("XPT", 0, Some(&plan)));
    }
}
