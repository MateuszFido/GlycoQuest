//! GlycoQuest post-filters and glycan annotation for xQuest hits.

use std::path::PathBuf;

use crate::cli::settings::Settings;
use crate::crosslinker::CrosslinkerProfile;
use crate::jobs::{JobManifest, VarModPlan};
use crate::prefilter::{FilteredSpectrum, MatchedIon, PrefilterResult};
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
    pub glyco_sites: Vec<GlycoSite>,
    pub all_sites_plausible: bool,
    pub n_glycan_pseudo: usize,
    pub matched_families: Vec<String>,
    pub matched_ions: Vec<MatchedIon>,
    pub matched_ion_count: usize,
    pub sequon_present: Option<bool>,
    pub charge_plausible: bool,
    pub hard_status: HardStatus,
    pub soft_score: f64,
    pub postfilter_status: PostfilterStatus,
}

/// One decoded glycan attachment site. Peptide positions are 1-based.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GlycoSite {
    pub peptide: u8,
    pub peptide_position: usize,
    pub residue: char,
    pub sequon_present: Option<bool>,
    pub plausible: bool,
}

/// Result of the hard (pass/fail) post-filter requirements.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HardStatus {
    Pass,
    FailNoXlink,
    FailNoGlycan,
    FailGlycanLimit,
    FailNoDiagnostic,
    FailPrecursorError,
    FailScore,
}

impl HardStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Pass => "pass",
            Self::FailNoXlink => "fail_no_xlink",
            Self::FailNoGlycan => "fail_no_glycan",
            Self::FailGlycanLimit => "fail_glycan_limit",
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
/// each job, and the prefilter supplies diagnostic-ion evidence per scan.
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

    let glyco_sites = plan
        .map(|plan| glyco_sites(&hit.seq1, &hit.seq2, plan))
        .unwrap_or_default();
    let glyco_peptide = glyco_sites.first().map(|site| site.peptide);
    let glyco_residue = glyco_sites.first().map(|site| site.residue);
    let all_sites_plausible = glyco_sites.len() == n_glycan_pseudo
        && !glyco_sites.is_empty()
        && glyco_sites.iter().all(|site| site.plausible);

    let scan = parse_scan(&hit.spectrum_id);
    let source_file = entry.map(|e| e.source_file.clone());
    let spectrum = lookup_spectrum(prefilter, source_file.as_deref(), scan);

    let matched_families = spectrum
        .map(|spec| spec.matched_families.clone())
        .unwrap_or_default();
    let matched_ions = spectrum
        .map(|spec| spec.matched_ions.clone())
        .unwrap_or_default();
    let matched_ion_count = matched_ions.len();

    let n_sites: Vec<&GlycoSite> = glyco_sites
        .iter()
        .filter(|site| site.residue == 'N')
        .collect();
    let sequon_present =
        (!n_sites.is_empty()).then(|| n_sites.iter().all(|site| site.sequon_present == Some(true)));
    let plausible_glycan_count = glyco_sites.iter().filter(|site| site.plausible).count();

    let charge_plausible = charge_plausible(hit.charge);

    let hard_status = hard_status(
        &hit,
        settings,
        n_glycan_pseudo,
        &glyco_sites,
        spectrum.is_some(),
        manifest.is_some(),
    );

    let soft_score = soft_score(
        &hit,
        plausible_glycan_count,
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
        glyco_sites,
        all_sites_plausible,
        n_glycan_pseudo,
        matched_families,
        matched_ions,
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
    glyco_sites: &[GlycoSite],
    diagnostic_positive: bool,
    have_manifest: bool,
) -> HardStatus {
    if hit.xlink_position.trim().is_empty() && hit.topology.trim().is_empty() {
        return HardStatus::FailNoXlink;
    }
    if n_glycan_pseudo == 0 {
        return HardStatus::FailNoGlycan;
    }
    let max = settings.max_glycans_per_peptide as usize;
    let peptide_1 = glyco_sites.iter().filter(|site| site.peptide == 1).count();
    let peptide_2 = glyco_sites.iter().filter(|site| site.peptide == 2).count();
    let decoded_all_sites = glyco_sites.len() == n_glycan_pseudo;
    if (decoded_all_sites && (peptide_1 > max || peptide_2 > max))
        || (!decoded_all_sites && n_glycan_pseudo > 2 * max)
    {
        return HardStatus::FailGlycanLimit;
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
    plausible_glycan_count: usize,
    charge_plausible: bool,
    matched_ion_count: usize,
) -> f64 {
    let mut score = hit.score;
    score += plausible_glycan_count as f64;
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

/// Decode every glycan pseudo-residue across both peptide sequences.
pub fn glyco_sites(seq1: &str, seq2: &str, plan: &VarModPlan) -> Vec<GlycoSite> {
    let glycan_pseudos = plan.glycan_pseudos();
    let mut sites = Vec::new();
    for (peptide, seq) in [(1u8, seq1), (2u8, seq2)] {
        for (pos, ch) in seq.chars().enumerate() {
            if glycan_pseudos.contains(&ch) {
                let residue = plan
                    .entry_for_pseudo(ch)
                    .map(|entry| entry.source_residue)
                    .unwrap_or(ch);
                let sequon_present = (residue == 'N').then(|| has_sequon(seq, pos, Some(plan)));
                let plausible = match residue {
                    'N' => sequon_present == Some(true),
                    'S' | 'T' => true,
                    _ => false,
                };
                sites.push(GlycoSite {
                    peptide,
                    peptide_position: pos + 1,
                    residue,
                    sequon_present,
                    plausible,
                });
            }
        }
    }
    sites
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
    (
        Some(job_id.trim_end_matches('_').to_string()),
        Some("none".to_string()),
    )
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

/// Read annotated rows from `glycoquest_xquest.csv`.
#[cfg(test)]
pub fn read_annotated_csv(path: &std::path::Path) -> Result<Vec<AnnotatedHit>, String> {
    let content = std::fs::read_to_string(path).map_err(|err| err.to_string())?;
    let mut lines = content.lines();
    let header = lines
        .next()
        .ok_or_else(|| format!("empty annotated CSV: {}", path.display()))?;
    let cols: Vec<&str> = header.split('\t').collect();
    let index = |name: &str| -> Result<usize, String> {
        cols.iter()
            .position(|col| *col == name)
            .ok_or_else(|| format!("missing column {name} in {}", path.display()))
    };
    let i_source = index("source_file")?;
    let i_scan = index("scan")?;
    let i_glycan_name = index("glycan_name")?;
    let i_glycan_composition = index("glycan_composition")?;
    let i_glycan_mass = index("glycan_mass")?;
    let i_loss = index("loss_label")?;
    let i_glyco_residue = index("glyco_residue")?;
    let i_glyco_peptide = index("glyco_peptide")?;
    let i_glyco_sites = cols.iter().position(|col| *col == "glyco_sites");
    let i_all_sites_plausible = cols.iter().position(|col| *col == "all_sites_plausible");
    let i_n_pseudo = index("n_glycan_pseudo")?;
    let i_sequon = index("sequon_present")?;
    let i_charge = index("charge")?;
    let i_charge_ok = index("charge_plausible")?;
    let i_families = index("matched_families")?;
    let i_ion_count = index("matched_ion_count")?;
    let i_matched_ions = cols.iter().position(|col| *col == "matched_ions");
    let i_link_type = cols.iter().position(|col| *col == "link_type");
    let i_seq1 = index("seq1")?;
    let i_seq2 = index("seq2")?;
    let i_prot1 = index("prot1")?;
    let i_prot2 = index("prot2")?;
    let i_topology = index("topology")?;
    let i_precursor = index("precursor_mz")?;
    let i_mr = index("mr")?;
    let i_ppm = index("precursor_error_ppm")?;
    let i_xlink = index("xlink_position")?;
    let i_xlinker_mass = cols.iter().position(|col| *col == "xlinker_mass");
    let i_xquest_version = cols.iter().position(|col| *col == "xquest_version");
    let i_xlinkions = cols.iter().position(|col| *col == "xlinkions_matched");
    let i_backboneions = cols.iter().position(|col| *col == "backboneions_matched");
    let i_xquest_matched_ions = cols.iter().position(|col| *col == "xquest_matched_ions");
    let i_score = index("score")?;
    let i_hard = index("hard_status")?;
    let i_soft = index("soft_score")?;
    let i_status = index("postfilter_status")?;

    let mut rows = Vec::new();
    for (line_no, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < cols.len() {
            return Err(format!(
                "row {} in {} has {} columns, expected {}",
                line_no + 2,
                path.display(),
                fields.len(),
                cols.len()
            ));
        }
        let composition = fields[i_glycan_composition].to_string();
        let loss_label = fields[i_loss].to_string();
        let job_id = job_id_from_composition_and_loss(&composition, &loss_label);
        let scan = parse_optional_u32(fields[i_scan]);
        let link_type = i_link_type
            .and_then(|idx| fields.get(idx))
            .and_then(|value| non_empty(value))
            .unwrap_or_else(|| infer_link_type(fields[i_seq2], fields[i_prot2]));
        let glyco_sites = i_glyco_sites
            .and_then(|idx| fields.get(idx))
            .map(|value| parse_glyco_sites(value))
            .unwrap_or_default();
        let all_sites_plausible = i_all_sites_plausible
            .and_then(|idx| fields.get(idx))
            .and_then(|value| parse_optional_bool(value))
            .unwrap_or_else(|| {
                !glyco_sites.is_empty() && glyco_sites.iter().all(|site| site.plausible)
            });
        rows.push(AnnotatedHit {
            hit: XQuestHit {
                spectrum_id: scan.map(|s| s.to_string()).unwrap_or_default(),
                search_hit_rank: 0,
                link_type,
                score: parse_f64(fields[i_score]),
                seq1: fields[i_seq1].to_string(),
                seq2: fields[i_seq2].to_string(),
                prot1: fields[i_prot1].to_string(),
                prot2: fields[i_prot2].to_string(),
                topology: fields[i_topology].to_string(),
                charge: parse_u8(fields[i_charge]),
                precursor_mz: parse_f64(fields[i_precursor]),
                mr: parse_f64(fields[i_mr]),
                precursor_error_ppm: parse_f64(fields[i_ppm]),
                xlink_position: fields[i_xlink].to_string(),
                xlinker_mass: i_xlinker_mass
                    .and_then(|idx| fields.get(idx))
                    .and_then(|value| parse_optional_f64(value)),
                xquest_version: i_xquest_version
                    .and_then(|idx| fields.get(idx))
                    .and_then(|value| non_empty(value)),
                xlinkions_matched: i_xlinkions
                    .and_then(|idx| fields.get(idx))
                    .and_then(|value| non_empty(value)),
                backboneions_matched: i_backboneions
                    .and_then(|idx| fields.get(idx))
                    .and_then(|value| non_empty(value)),
                matched_ions: i_xquest_matched_ions
                    .and_then(|idx| fields.get(idx))
                    .map(|value| parse_xquest_matched_ions(value))
                    .unwrap_or_default(),
                ..Default::default()
            },
            job_id,
            source_file: non_empty(fields[i_source]).map(PathBuf::from),
            scan,
            glycan_name: non_empty(fields[i_glycan_name]),
            glycan_composition: non_empty(&composition),
            glycan_mass: parse_optional_f64(fields[i_glycan_mass]),
            loss_label: non_empty(&loss_label),
            glyco_residue: fields[i_glyco_residue].chars().next(),
            glyco_peptide: parse_optional_u8(fields[i_glyco_peptide]),
            glyco_sites,
            all_sites_plausible,
            n_glycan_pseudo: parse_usize(fields[i_n_pseudo]),
            matched_families: fields[i_families]
                .split(';')
                .filter(|part| !part.is_empty())
                .map(str::to_string)
                .collect(),
            matched_ions: i_matched_ions
                .and_then(|idx| fields.get(idx))
                .map(|value| parse_matched_ions(value))
                .unwrap_or_default(),
            matched_ion_count: parse_usize(fields[i_ion_count]),
            sequon_present: parse_optional_bool(fields[i_sequon]),
            charge_plausible: fields[i_charge_ok].eq_ignore_ascii_case("true"),
            hard_status: hard_status_from_str(fields[i_hard], parse_usize(fields[i_n_pseudo])),
            soft_score: parse_f64(fields[i_soft]),
            postfilter_status: if fields[i_status] == "pass" {
                PostfilterStatus::Pass
            } else {
                PostfilterStatus::Fail
            },
        });
    }
    Ok(rows)
}

#[cfg(test)]
fn infer_link_type(seq2: &str, prot2: &str) -> String {
    if seq2.trim().is_empty() && prot2.trim().is_empty() {
        "monolink".to_string()
    } else {
        "crosslink".to_string()
    }
}

#[cfg(test)]
fn job_id_from_composition_and_loss(composition: &str, loss_label: &str) -> String {
    let mut stem = String::new();
    let re = regex_lite_composition(composition);
    for (name, count) in re {
        stem.push_str(&name);
        stem.push('_');
        stem.push_str(&count);
        stem.push('_');
    }
    match loss_label {
        "" | "none" => stem,
        label if label.starts_with('-') => format!("{stem}_{}", &label[1..]),
        label => format!("{stem}_{label}"),
    }
}

#[cfg(test)]
fn regex_lite_composition(composition: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let bytes = composition.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let start = i;
        while i < bytes.len() && bytes[i].is_ascii_alphabetic() {
            i += 1;
        }
        if start == i || i >= bytes.len() || bytes[i] != b'(' {
            break;
        }
        let name = composition[start..i].to_string();
        i += 1;
        let count_start = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        if count_start == i || i >= bytes.len() || bytes[i] != b')' {
            break;
        }
        let count = composition[count_start..i].to_string();
        i += 1;
        out.push((name, count));
    }
    out
}

#[cfg(test)]
pub(crate) fn parse_matched_ions(value: &str) -> Vec<MatchedIon> {
    value
        .split(';')
        .filter_map(|part| {
            let part = part.trim();
            if part.is_empty() {
                return None;
            }
            let (family, rest) = part.split_once('@')?;
            let (mz_text, loss_label) = match rest.split_once('[') {
                Some((mz, loss)) => (mz, loss.trim_end_matches(']').to_string()),
                None => (rest, String::new()),
            };
            let pieces: Vec<&str> = mz_text.split('|').collect();
            let observed_mz = pieces.first()?.parse::<f64>().ok()?;
            let expected_mz = pieces
                .get(1)
                .and_then(|value| value.parse::<f64>().ok())
                .unwrap_or(observed_mz);
            let peak_index = pieces
                .get(2)
                .and_then(|value| value.parse::<usize>().ok())
                .unwrap_or(0);
            let intensity = pieces
                .get(3)
                .and_then(|value| value.parse::<f64>().ok())
                .unwrap_or(0.0);
            let error_ppm = pieces
                .get(4)
                .and_then(|value| value.parse::<f64>().ok())
                .unwrap_or_else(|| {
                    if expected_mz == 0.0 {
                        0.0
                    } else {
                        ((observed_mz - expected_mz) / expected_mz) * 1_000_000.0
                    }
                });
            Some(MatchedIon {
                family: family.to_string(),
                expected_mz,
                observed_mz,
                loss_label,
                peak_index,
                intensity,
                error_ppm,
            })
        })
        .collect()
}

fn format_matched_ions(ions: &[MatchedIon]) -> String {
    ions.iter()
        .map(|ion| {
            let base = format!(
                "{}@{:.4}|{:.4}|{}|{:.4}|{:.3}",
                ion.family,
                ion.observed_mz,
                ion.expected_mz,
                ion.peak_index,
                ion.intensity,
                ion.error_ppm
            );
            if ion.loss_label.is_empty() {
                base
            } else {
                format!("{base}[{}]", ion.loss_label)
            }
        })
        .collect::<Vec<_>>()
        .join(";")
}

pub(crate) fn format_glyco_sites(sites: &[GlycoSite]) -> String {
    sites
        .iter()
        .map(|site| {
            format!(
                "pep{}:{}:{}:{}:{}",
                site.peptide,
                site.peptide_position,
                site.residue,
                opt_bool(site.sequon_present),
                site.plausible
            )
        })
        .collect::<Vec<_>>()
        .join(";")
}

#[cfg(test)]
fn parse_glyco_sites(value: &str) -> Vec<GlycoSite> {
    value
        .split(';')
        .filter_map(|part| {
            let fields: Vec<&str> = part.split(':').collect();
            if fields.len() != 5 {
                return None;
            }
            Some(GlycoSite {
                peptide: fields[0].strip_prefix("pep")?.parse().ok()?,
                peptide_position: fields[1].parse().ok()?,
                residue: fields[2].chars().next()?,
                sequon_present: parse_optional_bool(fields[3]),
                plausible: fields[4].eq_ignore_ascii_case("true"),
            })
        })
        .collect()
}

#[cfg(test)]
fn parse_xquest_matched_ions(value: &str) -> Vec<crate::results::extract::XQuestMatchedIon> {
    value
        .split(';')
        .filter_map(|part| {
            let fields: Vec<&str> = part.split('|').collect();
            if fields.len() < 6 {
                return None;
            }
            Some(crate::results::extract::XQuestMatchedIon {
                label: non_empty(fields[0]),
                ion_type: fields[1].to_string(),
                position: non_empty(fields[2]),
                theoretical_mz: parse_f64(fields[3]),
                observed_mz: parse_f64(fields[4]),
                delta_mz: parse_optional_f64(fields[5]),
                delta_ppm: fields.get(6).and_then(|value| parse_optional_f64(value)),
                intensity: fields.get(7).and_then(|value| parse_optional_f64(value)),
            })
        })
        .collect()
}

fn format_xquest_matched_ions(ions: &[crate::results::extract::XQuestMatchedIon]) -> String {
    ions.iter()
        .map(|ion| {
            [
                ion.label.clone().unwrap_or_default(),
                ion.ion_type.clone(),
                ion.position.clone().unwrap_or_default(),
                format!("{:.4}", ion.theoretical_mz),
                format!("{:.4}", ion.observed_mz),
                opt_f64(ion.delta_mz),
                opt_f64(ion.delta_ppm),
                opt_f64(ion.intensity),
            ]
            .join("|")
        })
        .collect::<Vec<_>>()
        .join(";")
}

#[cfg(test)]
fn hard_status_from_str(value: &str, n_glycan_pseudo: usize) -> HardStatus {
    match value {
        "pass" => HardStatus::Pass,
        "fail_no_xlink" => HardStatus::FailNoXlink,
        "fail_no_glycan" => HardStatus::FailNoGlycan,
        "fail_glycan_limit" => HardStatus::FailGlycanLimit,
        // Compatibility for annotated CSV files written before the outcome was split.
        "fail_glycan_count" if n_glycan_pseudo == 0 => HardStatus::FailNoGlycan,
        "fail_glycan_count" | "fail_multiple_glycans" => HardStatus::FailGlycanLimit,
        "fail_no_diagnostic" => HardStatus::FailNoDiagnostic,
        "fail_precursor_error" => HardStatus::FailPrecursorError,
        "fail_score" => HardStatus::FailScore,
        _ => HardStatus::FailScore,
    }
}

#[cfg(test)]
fn non_empty(value: &str) -> Option<String> {
    if value.is_empty() {
        None
    } else {
        Some(value.to_string())
    }
}

#[cfg(test)]
fn parse_optional_u32(value: &str) -> Option<u32> {
    if value.is_empty() {
        None
    } else {
        value.parse().ok()
    }
}

#[cfg(test)]
fn parse_optional_u8(value: &str) -> Option<u8> {
    if value.is_empty() {
        None
    } else {
        value.parse().ok()
    }
}

#[cfg(test)]
fn parse_optional_f64(value: &str) -> Option<f64> {
    if value.is_empty() {
        None
    } else {
        value.parse().ok()
    }
}

#[cfg(test)]
fn parse_optional_bool(value: &str) -> Option<bool> {
    if value.is_empty() {
        None
    } else {
        Some(value.eq_ignore_ascii_case("true"))
    }
}

#[cfg(test)]
fn parse_u8(value: &str) -> u8 {
    value.parse().unwrap_or(0)
}

#[cfg(test)]
fn parse_usize(value: &str) -> usize {
    value.parse().unwrap_or(0)
}

#[cfg(test)]
fn parse_f64(value: &str) -> f64 {
    value.parse().unwrap_or(0.0)
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
            "glyco_sites",
            "all_sites_plausible",
            "n_glycan_pseudo",
            "sequon_present",
            "charge",
            "charge_plausible",
            "matched_families",
            "matched_ion_count",
            "matched_ions",
            "link_type",
            "seq1",
            "seq2",
            "prot1",
            "prot2",
            "topology",
            "precursor_mz",
            "mr",
            "precursor_error_ppm",
            "xlink_position",
            "xlinker_mass",
            "xquest_version",
            "xlinkions_matched",
            "backboneions_matched",
            "num_matched_ions_alpha",
            "num_matched_ions_beta",
            "num_matched_common_ions_alpha",
            "num_matched_common_ions_beta",
            "num_matched_xlink_ions_alpha",
            "num_matched_xlink_ions_beta",
            "xquest_matched_ions",
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
                format_glyco_sites(&row.glyco_sites),
                row.all_sites_plausible.to_string(),
                row.n_glycan_pseudo.to_string(),
                opt_bool(row.sequon_present),
                hit.charge.to_string(),
                row.charge_plausible.to_string(),
                row.matched_families.join(";"),
                row.matched_ion_count.to_string(),
                format_matched_ions(&row.matched_ions),
                hit.normalized_link_type(),
                hit.seq1.clone(),
                hit.seq2.clone(),
                hit.prot1.clone(),
                hit.prot2.clone(),
                hit.topology.clone(),
                format!("{}", hit.precursor_mz),
                format!("{}", hit.mr),
                format!("{}", hit.precursor_error_ppm),
                hit.xlink_position.clone(),
                opt_f64(hit.xlinker_mass),
                hit.xquest_version.clone().unwrap_or_default(),
                hit.xlinkions_matched.clone().unwrap_or_default(),
                hit.backboneions_matched.clone().unwrap_or_default(),
                opt_u32(hit.num_matched_ions_alpha),
                opt_u32(hit.num_matched_ions_beta),
                opt_u32(hit.num_matched_common_ions_alpha),
                opt_u32(hit.num_matched_common_ions_beta),
                opt_u32(hit.num_matched_xlink_ions_alpha),
                opt_u32(hit.num_matched_xlink_ions_beta),
                format_xquest_matched_ions(&hit.matched_ions),
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
    use crate::jobs::{GlycanVariant, JobManifestEntry, SpectrumKey, build_varmod_plan};
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
                    peak_index: 0,
                    intensity: 5000.0,
                    error_ppm: 0.49,
                }],
            }],
            isotope_pairs: vec![],
            rejected: vec![],
            pruning: vec![],
            stats: Default::default(),
        }
    }

    #[test]
    fn multiple_occurrences_of_searched_glycan_pass_by_default() {
        let settings = Settings::defaults();
        let crosslinker = CrosslinkerProfile::resolve(&settings, Some("dss")).unwrap();
        let source = PathBuf::from("run.mzXML");
        let manifest = manifest_with("HexNAc_1_", source.clone());
        let prefilter = prefilter_with_scan(source, 42);
        let hits = vec![(
            "HexNAc_1_".to_string(),
            XQuestHit {
                seq1: "XATCXAT".into(),
                seq2: "PEPKIDE".into(),
                xlink_position: "2-1".into(),
                charge: 4,
                spectrum_id: "42".into(),
                ..Default::default()
            },
        )];
        let annotated =
            apply_postfilters(hits, &settings, &crosslinker, &prefilter, Some(&manifest));
        assert_eq!(annotated[0].hard_status, HardStatus::Pass);
        assert_eq!(annotated[0].glyco_sites.len(), 2);
        assert_eq!(annotated[0].postfilter_status, PostfilterStatus::Pass);
    }

    #[test]
    fn soft_score_adds_one_point_per_plausible_glycan_site() {
        let settings = Settings::defaults();
        let crosslinker = CrosslinkerProfile::resolve(&settings, Some("dss")).unwrap();
        let source = PathBuf::from("run.mzXML");
        let manifest = manifest_with("HexNAc_1_", source.clone());
        let prefilter = prefilter_with_scan(source, 42);
        let hits = vec![(
            "HexNAc_1_".to_string(),
            XQuestHit {
                seq1: "XATK".into(),
                seq2: "PEPKXPA".into(),
                xlink_position: "4-4".into(),
                charge: 4,
                spectrum_id: "42".into(),
                score: 10.0,
                ..Default::default()
            },
        )];

        let annotated =
            apply_postfilters(hits, &settings, &crosslinker, &prefilter, Some(&manifest));
        let row = &annotated[0];
        assert_eq!(row.hard_status, HardStatus::Pass);
        assert_eq!(row.glyco_sites.len(), 2);
        assert!(!row.all_sites_plausible);
        assert_eq!(row.glyco_sites[0].peptide_position, 1);
        assert_eq!(row.glyco_sites[1].peptide_position, 5);
        assert!((row.soft_score - 11.6).abs() < 1e-6);
    }

    #[test]
    fn implausible_n_glycan_site_does_not_fail_postfilter() {
        let settings = Settings::defaults();
        let crosslinker = CrosslinkerProfile::resolve(&settings, Some("dss")).unwrap();
        let source = PathBuf::from("run.mzXML");
        let manifest = manifest_with("HexNAc_1_", source.clone());
        let prefilter = prefilter_with_scan(source, 42);
        let hits = vec![(
            "HexNAc_1_".to_string(),
            XQuestHit {
                seq1: "XPAK".into(),
                seq2: "PEPKIDE".into(),
                xlink_position: "4-4".into(),
                charge: 4,
                spectrum_id: "42".into(),
                score: 10.0,
                ..Default::default()
            },
        )];

        let annotated =
            apply_postfilters(hits, &settings, &crosslinker, &prefilter, Some(&manifest));
        assert_eq!(annotated[0].hard_status, HardStatus::Pass);
        assert!(!annotated[0].all_sites_plausible);
        assert!((annotated[0].soft_score - 10.6).abs() < 1e-6);
    }

    #[test]
    fn per_peptide_glycan_cap_is_enforced() {
        let mut settings = Settings::defaults();
        settings.max_glycans_per_peptide = 3;
        let crosslinker = CrosslinkerProfile::resolve(&settings, Some("dss")).unwrap();
        let source = PathBuf::from("run.mzXML");
        let manifest = manifest_with("HexNAc_1_", source.clone());
        let prefilter = prefilter_with_scan(source, 42);
        let hits = vec![(
            "HexNAc_1_".to_string(),
            XQuestHit {
                seq1: "XATXATXATXATK".into(),
                seq2: "PEPKIDE".into(),
                xlink_position: "13-4".into(),
                charge: 4,
                spectrum_id: "42".into(),
                score: 10.0,
                ..Default::default()
            },
        )];

        let annotated =
            apply_postfilters(hits, &settings, &crosslinker, &prefilter, Some(&manifest));
        assert_eq!(annotated[0].hard_status, HardStatus::FailGlycanLimit);
        assert_eq!(annotated[0].glyco_sites.len(), 4);
        assert!(annotated[0].all_sites_plausible);
    }

    #[test]
    fn three_glycans_on_each_peptide_pass_at_default_cap() {
        let settings = Settings::defaults();
        let crosslinker = CrosslinkerProfile::resolve(&settings, Some("dss")).unwrap();
        let source = PathBuf::from("run.mzXML");
        let manifest = manifest_with("HexNAc_1_", source.clone());
        let prefilter = prefilter_with_scan(source, 42);
        let hits = vec![(
            "HexNAc_1_".to_string(),
            XQuestHit {
                seq1: "XATXATXATK".into(),
                seq2: "XASXASXASK".into(),
                xlink_position: "10-10".into(),
                charge: 6,
                spectrum_id: "42".into(),
                score: 10.0,
                ..Default::default()
            },
        )];

        let annotated =
            apply_postfilters(hits, &settings, &crosslinker, &prefilter, Some(&manifest));
        assert_eq!(annotated[0].hard_status, HardStatus::Pass);
        assert_eq!(annotated[0].glyco_sites.len(), 6);
        assert_eq!(annotated[0].postfilter_status, PostfilterStatus::Pass);
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
                seq1: "XATDE".into(),
                seq2: "PEPKIDE".into(),
                xlink_position: "3-1".into(),
                charge: 4,
                spectrum_id: "42".into(),
                score: 5.0,
                ..Default::default()
            },
        )];
        let annotated =
            apply_postfilters(hits, &settings, &crosslinker, &prefilter, Some(&manifest));
        let row = &annotated[0];
        assert_eq!(row.hard_status, HardStatus::Pass);
        assert_eq!(row.glycan_composition.as_deref(), Some("HexNAc(1)"));
        assert_eq!(row.glyco_residue, Some('N'));
        assert_eq!(row.glyco_peptide, Some(1));
        assert_eq!(row.scan, Some(42));
        assert_eq!(row.n_glycan_pseudo, 1);
    }

    #[test]
    fn glycosylated_monolink_passes_without_second_peptide() {
        let settings = Settings::defaults();
        let crosslinker = CrosslinkerProfile::resolve(&settings, Some("dss")).unwrap();
        let source = PathBuf::from("run.mzXML");
        let manifest = manifest_with("HexNAc_1_", source.clone());
        let prefilter = prefilter_with_scan(source, 42);
        let hits = vec![(
            "HexNAc_1_".to_string(),
            XQuestHit {
                link_type: "monolink".into(),
                seq1: "XATK".into(),
                seq2: "".into(),
                prot1: "FETUA_BOVIN".into(),
                prot2: "".into(),
                xlink_position: "4".into(),
                xlinker_mass: Some(156.07864),
                charge: 4,
                spectrum_id: "42".into(),
                score: 5.0,
                ..Default::default()
            },
        )];

        let annotated =
            apply_postfilters(hits, &settings, &crosslinker, &prefilter, Some(&manifest));
        let row = &annotated[0];

        assert_eq!(row.hard_status, HardStatus::Pass);
        assert_eq!(row.postfilter_status, PostfilterStatus::Pass);
        assert_eq!(row.hit.link_type, "monolink");
        assert_eq!(row.hit.seq2, "");
        assert_eq!(row.hit.prot2, "");
        assert_eq!(row.glyco_peptide, Some(1));
        assert_eq!(row.n_glycan_pseudo, 1);
    }

    #[test]
    fn annotated_csv_roundtrips_monolink_fields() {
        let settings = Settings::defaults();
        let crosslinker = CrosslinkerProfile::resolve(&settings, Some("dss")).unwrap();
        let source = PathBuf::from("run.mzXML");
        let manifest = manifest_with("HexNAc_1_", source.clone());
        let prefilter = prefilter_with_scan(source, 42);
        let hits = vec![(
            "HexNAc_1_".to_string(),
            XQuestHit {
                link_type: "monolink".into(),
                seq1: "XATK".into(),
                prot1: "FETUA_BOVIN".into(),
                xlink_position: "4".into(),
                xlinker_mass: Some(156.07864),
                charge: 4,
                spectrum_id: "42".into(),
                score: 5.0,
                ..Default::default()
            },
        )];
        let annotated =
            apply_postfilters(hits, &settings, &crosslinker, &prefilter, Some(&manifest));
        let dir = std::env::temp_dir().join(format!("glycoquest_mono_csv_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("glycoquest_xquest.csv");

        write_annotated_csv(&path, &annotated).unwrap();
        let roundtrip = read_annotated_csv(&path).unwrap();

        assert_eq!(roundtrip.len(), 1);
        assert_eq!(roundtrip[0].hit.link_type, "monolink");
        assert_eq!(roundtrip[0].hit.prot2, "");
        assert_eq!(roundtrip[0].hit.seq2, "");
        assert!((roundtrip[0].hit.xlinker_mass.unwrap() - 156.07864).abs() < 1e-6);
        assert_eq!(roundtrip[0].postfilter_status, PostfilterStatus::Pass);
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn non_glycosylated_pair_reports_no_glycan() {
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
        let annotated =
            apply_postfilters(hits, &settings, &crosslinker, &prefilter, Some(&manifest));
        assert_eq!(annotated[0].hard_status, HardStatus::FailNoGlycan);
        assert_eq!(annotated[0].hard_status.as_str(), "fail_no_glycan");
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
