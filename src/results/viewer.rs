//! Interactive viewer bundle (`results/viewer/viewer.json`) for the GlycoQuest CLMS viewer.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::Path;

use serde::Serialize;

use crate::cli::settings::Settings;
use crate::crosslinker::CrosslinkerProfile;
use crate::fasta::FastaDatabase;
use crate::jobs::JobManifest;
use crate::prefilter::MatchedIon;
use crate::prefilter::PrefilterStats;

use super::filter::{AnnotatedHit, PostfilterStatus};
use super::mapping::{map_crosslink, protein_lookup};
use super::report::ReportContext;
use super::spectrum::{ScanSpectrum, load_spectra_for_scans, scans_from_hits};

pub const VIEWER_SCHEMA_VERSION: u32 = 4;

#[derive(Debug, Clone, Serialize)]
pub struct ViewerBundle {
    pub viewer_schema_version: u32,
    pub meta: ViewerMeta,
    pub proteins: Vec<ViewerProtein>,
    pub crosslinks: Vec<ViewerCrosslink>,
    pub qc: ViewerQc,
    pub spectra: HashMap<String, ViewerSpectrum>,
    pub isotope_pairs: HashMap<String, ViewerIsotopePair>,
    pub filtering: HashMap<String, ViewerFiltering>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ViewerMeta {
    pub project: String,
    pub input_label: String,
    pub crosslinker: String,
    pub crosslinker_mw: f64,
    pub xlink_sites: String,
    pub glycan_library: String,
    pub xquest_version: Option<String>,
    pub generated_at: String,
    pub generated_at_iso: String,
    pub generated_at_unix: Option<u64>,
    pub total_hits: usize,
    pub passing_hits: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct ViewerProtein {
    pub id: String,
    pub display_name: String,
    pub sequence: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct ViewerCrosslink {
    pub id: String,
    pub link_type: String,
    pub protein1: String,
    pub pep_pos1: Option<usize>,
    pub pep_seq1: String,
    pub link_pos1: Option<usize>,
    pub abs_pos1: Option<usize>,
    pub protein2: String,
    pub pep_pos2: Option<usize>,
    pub pep_seq2: String,
    pub link_pos2: Option<usize>,
    pub abs_pos2: Option<usize>,
    pub score: f64,
    pub soft_score: f64,
    pub scan: Option<u32>,
    pub retention_time_min: Option<f64>,
    pub source_file: Option<String>,
    pub charge: u8,
    pub precursor_mz: f64,
    pub precursor_error_ppm: f64,
    pub xlinker_mass: Option<f64>,
    pub topology: String,
    pub protein_pair_key: String,
    pub glycan_name: Option<String>,
    pub glycan_composition: Option<String>,
    pub glyco_residue: Option<char>,
    pub glyco_peptide: Option<u8>,
    pub glyco_sites: Vec<ViewerGlycoSite>,
    pub diagnostic_ions: Vec<ViewerDiagnosticIon>,
    pub loss_label: Option<String>,
    pub postfilter_status: String,
    pub mapped: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct ViewerGlycoSite {
    pub peptide: u8,
    pub peptide_position: usize,
    pub residue: char,
    pub sequon_present: Option<bool>,
    pub plausible: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct ViewerDiagnosticIon {
    pub family: String,
    pub expected_mz: f64,
    pub observed_mz: f64,
    pub loss_label: String,
    pub peak_index: usize,
    pub intensity: f64,
    pub error_ppm: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct ViewerQc {
    pub funnel: Vec<NamedCount>,
    pub outcomes: Vec<NamedCount>,
    pub glycan_top: Vec<NamedCount>,
    pub site_dist: Vec<NamedCount>,
    pub score_hist: Histogram,
    pub ppm_hist: Histogram,
}

#[derive(Debug, Clone, Serialize)]
pub struct NamedCount {
    pub label: String,
    pub count: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct Histogram {
    pub bins: usize,
    pub min: f64,
    pub max: f64,
    pub counts: Vec<usize>,
    pub n: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct ViewerSpectrum {
    pub mz: Vec<f32>,
    pub intensity: Vec<f32>,
    pub retention_time_min: f64,
    pub precursor_mz: f64,
    pub charge: u8,
}

#[derive(Debug, Clone, Serialize)]
pub struct ViewerIsotopePair {
    pub id: String,
    pub source_artifact: String,
    pub source_row: Option<usize>,
    pub light_file: Option<String>,
    pub heavy_file: Option<String>,
    pub light_scan: u32,
    pub heavy_scan: u32,
    pub rt_light_min: f64,
    pub rt_heavy_min: f64,
    pub mz_light: f64,
    pub mz_heavy: f64,
    pub light_charge: u8,
    pub heavy_charge: u8,
}

#[derive(Debug, Clone, Serialize)]
pub struct ViewerFiltering {
    pub input_scan: ViewerFilteringInputScan,
    pub diagnostic_prefilter: ViewerFilteringDiagnosticPrefilter,
    pub isotope_pair: Option<ViewerFilteringIsotopePair>,
    pub glycan_pruning: ViewerFilteringGlycanPruning,
    pub xquest_search: ViewerFilteringXquestSearch,
    pub postfilter: ViewerFilteringPostfilter,
}

#[derive(Debug, Clone, Serialize)]
pub struct ViewerFilteringInputScan {
    pub status: String,
    pub source_file: Option<String>,
    pub source_artifact: String,
    pub source_row: Option<usize>,
    pub scan: Option<u32>,
    pub retention_time_min: Option<f64>,
    pub precursor_mz: f64,
    pub charge: u8,
    pub peak_count: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct ViewerFilteringDiagnosticPrefilter {
    pub status: String,
    pub source_artifact: String,
    pub source_row: Option<usize>,
    pub matched_family_count: usize,
    pub matched_families: Vec<String>,
    pub matched_ions: Vec<ViewerDiagnosticIon>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ViewerFilteringIsotopePair {
    pub status: String,
    pub source_artifact: String,
    pub source_row: Option<usize>,
    pub light_scan: u32,
    pub heavy_scan: u32,
    pub rt_light_min: f64,
    pub rt_heavy_min: f64,
    pub mz_light: f64,
    pub mz_heavy: f64,
    pub light_charge: u8,
    pub heavy_charge: u8,
}

#[derive(Debug, Clone, Serialize)]
pub struct ViewerFilteringGlycanPruning {
    pub status: String,
    pub source_artifact: String,
    pub source_row: Option<usize>,
    pub selected_glycan: Option<String>,
    pub selected_composition: Option<String>,
    pub retained_count_for_scan: usize,
    pub required_families: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ViewerFilteringXquestSearch {
    pub status: String,
    pub source_artifact: String,
    pub source_row: Option<usize>,
    pub xquest_version: Option<String>,
    pub score: f64,
    pub rank: u32,
    pub xlinkions_matched: Option<String>,
    pub backboneions_matched: Option<String>,
    pub num_matched_ions_alpha: Option<u32>,
    pub num_matched_ions_beta: Option<u32>,
    pub num_matched_common_ions_alpha: Option<u32>,
    pub num_matched_common_ions_beta: Option<u32>,
    pub num_matched_xlink_ions_alpha: Option<u32>,
    pub num_matched_xlink_ions_beta: Option<u32>,
    pub matched_ions: Vec<ViewerXquestMatchedIon>,
    pub unavailable_reason: Option<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ViewerXquestMatchedIon {
    pub label: String,
    pub ion_type: String,
    pub peptide: Option<String>,
    pub position: Option<String>,
    pub theoretical_mz: f64,
    pub observed_mz: f64,
    pub error_da: Option<f64>,
    pub error_ppm: Option<f64>,
    pub intensity: Option<f64>,
    pub peak_index: Option<usize>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ViewerFilteringPostfilter {
    pub status: String,
    pub source_artifact: String,
    pub source_row: Option<usize>,
    pub hard_status: String,
    pub rules: Vec<ViewerFilteringRule>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ViewerFilteringRule {
    pub name: String,
    pub status: String,
    pub value: String,
    pub threshold: String,
}

/// Build the viewer JSON bundle from consolidated run artifacts.
pub fn build_viewer_bundle(
    hits: &[AnnotatedHit],
    stats: &PrefilterStats,
    manifest: Option<&JobManifest>,
    fasta: &FastaDatabase,
    crosslinker: &CrosslinkerProfile,
    settings: &Settings,
    ctx: &ReportContext,
    spectra_dir: &Path,
) -> ViewerBundle {
    let proteins_map = protein_lookup(fasta);
    let protein_ids = collect_protein_ids(hits, &proteins_map);
    let proteins = build_protein_list(fasta, &protein_ids);

    let mut scans_needed = scans_from_hits(hits, false);
    let isotope_pairs_by_scan = load_isotope_pairs_for_viewer(spectra_dir);
    for scan in scans_needed.iter().copied().collect::<Vec<_>>() {
        if let Some(pair) = isotope_pairs_by_scan.get(&scan) {
            scans_needed.insert(pair.light_scan);
            scans_needed.insert(pair.heavy_scan);
        }
    }
    let spectra_loaded = load_spectra_for_scans(spectra_dir, &scans_needed);
    let artifact_rows = spectra_dir
        .parent()
        .map(load_filtering_artifact_rows)
        .unwrap_or_default();

    let mut crosslinks = Vec::with_capacity(hits.len());
    let mut filtering = HashMap::new();

    for (idx, hit) in hits.iter().enumerate() {
        let plan = manifest
            .and_then(|m| m.by_job_id(&hit.job_id))
            .map(|entry| &entry.varmod_plan);
        let mapping = map_crosslink(
            &hit.hit.seq1,
            &hit.hit.seq2,
            &hit.hit.prot1,
            &hit.hit.prot2,
            &hit.hit.xlink_position,
            &proteins_map,
            plan,
        );
        let id = crosslink_id(hit, idx);
        let link_type = hit.hit.normalized_link_type();
        let is_monolink = hit.hit.is_monolink();
        let mapped = mapping.abs1.is_some() && (is_monolink || mapping.abs2.is_some());

        let scan_spectrum = hit.scan.and_then(|s| spectra_loaded.get(&s));
        let protein_pair_key = protein_pair_key(&mapping.prot1, &mapping.prot2);
        filtering.insert(
            id.clone(),
            build_filtering(
                hit,
                scan_spectrum,
                isotope_pairs_by_scan.get(&hit.scan.unwrap_or_default()),
                &artifact_rows,
                settings,
                idx + 2,
            ),
        );

        crosslinks.push(ViewerCrosslink {
            id,
            link_type,
            protein1: mapping.prot1,
            pep_pos1: mapping.pep_pos1,
            pep_seq1: mapping.pep1,
            link_pos1: mapping.link1,
            abs_pos1: mapping.abs1,
            protein2: mapping.prot2,
            pep_pos2: mapping.pep_pos2,
            pep_seq2: mapping.pep2,
            link_pos2: mapping.link2,
            abs_pos2: mapping.abs2,
            score: hit.hit.score,
            soft_score: hit.soft_score,
            scan: hit.scan,
            retention_time_min: scan_spectrum.map(|s| s.retention_time_min),
            source_file: hit
                .source_file
                .as_ref()
                .map(|path| path.display().to_string()),
            charge: hit.hit.charge,
            precursor_mz: hit.hit.precursor_mz,
            precursor_error_ppm: hit.hit.precursor_error_ppm,
            xlinker_mass: hit.hit.xlinker_mass,
            topology: hit.hit.topology.clone(),
            protein_pair_key,
            glycan_name: hit.glycan_name.clone(),
            glycan_composition: hit.glycan_composition.clone(),
            glyco_residue: hit.glyco_residue,
            glyco_peptide: hit.glyco_peptide,
            glyco_sites: hit
                .glyco_sites
                .iter()
                .map(|site| ViewerGlycoSite {
                    peptide: site.peptide,
                    peptide_position: site.peptide_position,
                    residue: site.residue,
                    sequon_present: site.sequon_present,
                    plausible: site.plausible,
                })
                .collect(),
            diagnostic_ions: hit
                .matched_ions
                .iter()
                .map(|ion| ViewerDiagnosticIon {
                    family: ion.family.clone(),
                    expected_mz: ion.expected_mz,
                    observed_mz: ion.observed_mz,
                    loss_label: ion.loss_label.clone(),
                    peak_index: ion.peak_index,
                    intensity: ion.intensity,
                    error_ppm: ion.error_ppm,
                })
                .collect(),
            loss_label: hit.loss_label.clone(),
            postfilter_status: hit.postfilter_status.as_str().to_string(),
            mapped,
        });
    }

    let passing = hits
        .iter()
        .filter(|h| h.postfilter_status == PostfilterStatus::Pass)
        .count();

    let spectra = spectra_loaded
        .into_iter()
        .map(|(scan, spec)| {
            (
                scan.to_string(),
                ViewerSpectrum {
                    mz: spec.mz,
                    intensity: spec.intensity,
                    retention_time_min: spec.retention_time_min,
                    precursor_mz: spec.precursor_mz,
                    charge: spec.charge,
                },
            )
        })
        .collect();
    let isotope_pairs = isotope_pairs_by_scan
        .into_iter()
        .filter(|(scan, _)| scans_needed.contains(scan))
        .map(|(scan, pair)| (scan.to_string(), pair))
        .collect();

    let generated_at_unix = unix_timestamp();
    let generated_at_iso = iso8601_utc_from_unix(generated_at_unix);

    ViewerBundle {
        viewer_schema_version: VIEWER_SCHEMA_VERSION,
        meta: ViewerMeta {
            project: ctx.project.clone(),
            input_label: ctx.input_label.clone(),
            crosslinker: ctx.crosslinker_name.clone(),
            crosslinker_mw: crosslinker.xlinkermw,
            xlink_sites: ctx.xlink_sites.clone(),
            glycan_library: ctx.glycan_library.clone(),
            xquest_version: hits.iter().find_map(|hit| hit.hit.xquest_version.clone()),
            generated_at: generated_at_iso.clone(),
            generated_at_iso,
            generated_at_unix: Some(generated_at_unix),
            total_hits: hits.len(),
            passing_hits: passing,
        },
        proteins,
        crosslinks,
        qc: build_qc(hits, stats),
        spectra,
        isotope_pairs,
        filtering,
    }
}

#[derive(Debug, Clone, Default)]
struct FilteringArtifactRows {
    filtered_rows: HashMap<String, usize>,
    pruning_rows: HashMap<String, usize>,
    pruning_counts: HashMap<String, usize>,
}

fn build_filtering(
    hit: &AnnotatedHit,
    spectrum: Option<&ScanSpectrum>,
    isotope_pair: Option<&ViewerIsotopePair>,
    artifact_rows: &FilteringArtifactRows,
    settings: &Settings,
    result_row: usize,
) -> ViewerFiltering {
    let scan_key = source_scan_key(hit.source_file.as_deref(), hit.scan);
    let filtered_row = scan_key
        .as_ref()
        .and_then(|key| artifact_rows.filtered_rows.get(key).copied());
    let pruning_row = scan_key.as_ref().and_then(|key| {
        let selected = hit.glycan_composition.as_deref().unwrap_or_default();
        artifact_rows
            .pruning_rows
            .get(&format!("{key}|{selected}"))
            .copied()
    });
    let pruning_count = scan_key
        .as_ref()
        .and_then(|key| artifact_rows.pruning_counts.get(key).copied())
        .unwrap_or(0);

    ViewerFiltering {
        input_scan: ViewerFilteringInputScan {
            status: if hit.scan.is_some() {
                "available".into()
            } else {
                "unavailable".into()
            },
            source_file: hit
                .source_file
                .as_ref()
                .map(|path| path.display().to_string()),
            source_artifact: "spectra/".into(),
            source_row: None,
            scan: hit.scan,
            retention_time_min: spectrum.map(|s| s.retention_time_min),
            precursor_mz: spectrum
                .map(|s| s.precursor_mz)
                .unwrap_or(hit.hit.precursor_mz),
            charge: spectrum.map(|s| s.charge).unwrap_or(hit.hit.charge),
            peak_count: spectrum.map(|s| s.mz.len()).unwrap_or(0),
        },
        diagnostic_prefilter: ViewerFilteringDiagnosticPrefilter {
            status: if hit.matched_ions.is_empty() {
                "unavailable".into()
            } else {
                "matched".into()
            },
            source_artifact: "filtered_spectra.tsv".into(),
            source_row: filtered_row,
            matched_family_count: hit.matched_families.len(),
            matched_families: hit.matched_families.clone(),
            matched_ions: hit.matched_ions.iter().map(viewer_diagnostic_ion).collect(),
        },
        isotope_pair: isotope_pair.map(|pair| ViewerFilteringIsotopePair {
            status: "matched".into(),
            source_artifact: pair.source_artifact.clone(),
            source_row: pair.source_row,
            light_scan: pair.light_scan,
            heavy_scan: pair.heavy_scan,
            rt_light_min: pair.rt_light_min,
            rt_heavy_min: pair.rt_heavy_min,
            mz_light: pair.mz_light,
            mz_heavy: pair.mz_heavy,
            light_charge: pair.light_charge,
            heavy_charge: pair.heavy_charge,
        }),
        glycan_pruning: ViewerFilteringGlycanPruning {
            status: if pruning_count > 0 {
                "retained".into()
            } else {
                "unavailable".into()
            },
            source_artifact: "glycan_pruning.tsv".into(),
            source_row: pruning_row,
            selected_glycan: hit.glycan_name.clone(),
            selected_composition: hit.glycan_composition.clone(),
            retained_count_for_scan: pruning_count,
            required_families: hit.matched_families.clone(),
        },
        xquest_search: build_xquest_filtering(hit, spectrum),
        postfilter: ViewerFilteringPostfilter {
            status: hit.postfilter_status.as_str().to_string(),
            source_artifact: "results/glycoquest_xquest.csv".into(),
            source_row: Some(result_row),
            hard_status: hit.hard_status.as_str().to_string(),
            rules: postfilter_rules(hit, settings),
        },
    }
}

fn viewer_diagnostic_ion(ion: &MatchedIon) -> ViewerDiagnosticIon {
    ViewerDiagnosticIon {
        family: ion.family.clone(),
        expected_mz: ion.expected_mz,
        observed_mz: ion.observed_mz,
        loss_label: ion.loss_label.clone(),
        peak_index: ion.peak_index,
        intensity: ion.intensity,
        error_ppm: ion.error_ppm,
    }
}

fn build_xquest_filtering(
    hit: &AnnotatedHit,
    spectrum: Option<&ScanSpectrum>,
) -> ViewerFilteringXquestSearch {
    let matched_ions: Vec<ViewerXquestMatchedIon> = hit
        .hit
        .matched_ions
        .iter()
        .map(|ion| ViewerXquestMatchedIon {
            label: ion
                .label
                .clone()
                .unwrap_or_else(|| xquest_ion_label(&ion.ion_type, ion.position.as_deref())),
            ion_type: ion.ion_type.clone(),
            peptide: xquest_peptide_from_ion_type(&ion.ion_type),
            position: ion.position.clone(),
            theoretical_mz: ion.theoretical_mz,
            observed_mz: ion.observed_mz,
            error_da: ion.delta_mz.or(Some(ion.observed_mz - ion.theoretical_mz)),
            error_ppm: ion.delta_ppm.or_else(|| {
                (ion.theoretical_mz != 0.0).then(|| {
                    ((ion.observed_mz - ion.theoretical_mz) / ion.theoretical_mz) * 1_000_000.0
                })
            }),
            intensity: ion.intensity,
            peak_index: spectrum.and_then(|s| exact_peak_index(s, ion.observed_mz)),
        })
        .collect();
    let unavailable_reason = if matched_ions.is_empty() {
        Some("xQuest exact matched-ion rows were not emitted".into())
    } else {
        None
    };
    ViewerFilteringXquestSearch {
        status: if matched_ions.is_empty() {
            "unavailable".into()
        } else {
            "matched".into()
        },
        source_artifact: "jobs/<job_id>/results/xquest.xml".replace("<job_id>", &hit.job_id),
        source_row: None,
        xquest_version: hit.hit.xquest_version.clone(),
        score: hit.hit.score,
        rank: hit.hit.search_hit_rank,
        xlinkions_matched: hit.hit.xlinkions_matched.clone(),
        backboneions_matched: hit.hit.backboneions_matched.clone(),
        num_matched_ions_alpha: hit.hit.num_matched_ions_alpha,
        num_matched_ions_beta: hit.hit.num_matched_ions_beta,
        num_matched_common_ions_alpha: hit.hit.num_matched_common_ions_alpha,
        num_matched_common_ions_beta: hit.hit.num_matched_common_ions_beta,
        num_matched_xlink_ions_alpha: hit.hit.num_matched_xlink_ions_alpha,
        num_matched_xlink_ions_beta: hit.hit.num_matched_xlink_ions_beta,
        matched_ions,
        unavailable_reason,
    }
}

fn xquest_ion_label(ion_type: &str, position: Option<&str>) -> String {
    match position {
        Some(pos) if !pos.is_empty() => format!("{ion_type}{pos}"),
        _ => ion_type.to_string(),
    }
}

fn xquest_peptide_from_ion_type(ion_type: &str) -> Option<String> {
    let lower = ion_type.to_ascii_lowercase();
    if lower.contains("beta") || lower.starts_with("p2_") {
        Some("beta".into())
    } else if lower.contains("alpha") || lower.starts_with("p1_") {
        Some("alpha".into())
    } else {
        None
    }
}

fn exact_peak_index(spectrum: &ScanSpectrum, observed_mz: f64) -> Option<usize> {
    spectrum
        .mz
        .iter()
        .position(|mz| ((*mz as f64) - observed_mz).abs() <= 1e-4)
}

fn postfilter_rules(hit: &AnnotatedHit, settings: &Settings) -> Vec<ViewerFilteringRule> {
    let peptide_1_glycans = hit
        .glyco_sites
        .iter()
        .filter(|site| site.peptide == 1)
        .count();
    let peptide_2_glycans = hit
        .glyco_sites
        .iter()
        .filter(|site| site.peptide == 2)
        .count();
    let max = settings.max_glycans_per_peptide as usize;
    let decoded_all_sites = hit.glyco_sites.len() == hit.n_glycan_pseudo;
    let glycan_count_passes = hit.n_glycan_pseudo > 0
        && if decoded_all_sites {
            peptide_1_glycans <= max && peptide_2_glycans <= max
        } else {
            hit.n_glycan_pseudo <= 2 * max
        };
    let plausible_sites = hit.glyco_sites.iter().filter(|site| site.plausible).count();

    vec![
        ViewerFilteringRule {
            name: "xlink_position".into(),
            status: if hit.hit.xlink_position.trim().is_empty()
                && hit.hit.topology.trim().is_empty()
            {
                "fail".into()
            } else {
                "pass".into()
            },
            value: if hit.hit.xlink_position.trim().is_empty() {
                hit.hit.topology.clone()
            } else {
                hit.hit.xlink_position.clone()
            },
            threshold: "required".into(),
        },
        ViewerFilteringRule {
            name: "n_glycan_pseudo".into(),
            status: if glycan_count_passes { "pass" } else { "fail" }.into(),
            value: hit.n_glycan_pseudo.to_string(),
            threshold: format!(">= 1 total; <= {max} per peptide"),
        },
        ViewerFilteringRule {
            name: "plausible_glycan_sites".into(),
            status: "soft".into(),
            value: format!("{plausible_sites}/{}", hit.glyco_sites.len()),
            threshold: "+1.0 soft-score point each".into(),
        },
        ViewerFilteringRule {
            name: "diagnostic_positive".into(),
            status: if hit.matched_ions.is_empty() {
                "fail"
            } else {
                "pass"
            }
            .into(),
            value: (!hit.matched_ions.is_empty()).to_string(),
            threshold: "required".into(),
        },
        ViewerFilteringRule {
            name: "precursor_error_ppm".into(),
            status: if hit.hit.precursor_error_ppm.abs() <= settings.max_precursor_error_ppm {
                "pass"
            } else {
                "fail"
            }
            .into(),
            value: format!("{:.3}", hit.hit.precursor_error_ppm),
            threshold: format!("abs <= {:.3}", settings.max_precursor_error_ppm),
        },
        ViewerFilteringRule {
            name: "score".into(),
            status: if hit.hit.score >= settings.min_score {
                "pass"
            } else {
                "fail"
            }
            .into(),
            value: format!("{:.3}", hit.hit.score),
            threshold: format!(">= {:.3}", settings.min_score),
        },
    ]
}

fn source_scan_key(source_file: Option<&Path>, scan: Option<u32>) -> Option<String> {
    Some(format!("{}|{}", source_file?.display(), scan?))
}

fn load_filtering_artifact_rows(project_dir: &Path) -> FilteringArtifactRows {
    let mut rows = FilteringArtifactRows::default();
    load_filtered_rows(project_dir, &mut rows);
    load_pruning_rows(project_dir, &mut rows);
    rows
}

fn load_filtered_rows(project_dir: &Path, rows: &mut FilteringArtifactRows) {
    let Ok(content) = std::fs::read_to_string(project_dir.join("filtered_spectra.tsv")) else {
        return;
    };
    let mut lines = content.lines();
    let Some(header) = lines.next() else {
        return;
    };
    let cols: Vec<&str> = header.split('\t').collect();
    let (Some(i_source), Some(i_scan)) = (column(&cols, "source_file"), column(&cols, "scan"))
    else {
        return;
    };
    for (idx, line) in lines.enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        let (Some(source), Some(scan)) = (fields.get(i_source), fields.get(i_scan)) else {
            continue;
        };
        rows.filtered_rows
            .insert(format!("{source}|{scan}"), idx + 2);
    }
}

fn load_pruning_rows(project_dir: &Path, rows: &mut FilteringArtifactRows) {
    let Ok(content) = std::fs::read_to_string(project_dir.join("glycan_pruning.tsv")) else {
        return;
    };
    let mut lines = content.lines();
    let Some(header) = lines.next() else {
        return;
    };
    let cols: Vec<&str> = header.split('\t').collect();
    let (Some(i_source), Some(i_scan), Some(i_comp)) = (
        column(&cols, "source_file"),
        column(&cols, "scan"),
        column(&cols, "composition"),
    ) else {
        return;
    };
    for (idx, line) in lines.enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        let (Some(source), Some(scan), Some(comp)) =
            (fields.get(i_source), fields.get(i_scan), fields.get(i_comp))
        else {
            continue;
        };
        let key = format!("{source}|{scan}");
        *rows.pruning_counts.entry(key.clone()).or_insert(0) += 1;
        rows.pruning_rows
            .entry(format!("{key}|{comp}"))
            .or_insert(idx + 2);
    }
}

fn load_isotope_pairs_for_viewer(spectra_dir: &Path) -> HashMap<u32, ViewerIsotopePair> {
    let Some(project_dir) = spectra_dir.parent() else {
        return HashMap::new();
    };
    let path = project_dir.join("isotope_pairs.tsv");
    let Ok(content) = std::fs::read_to_string(path) else {
        return HashMap::new();
    };
    let mut lines = content.lines();
    let Some(header) = lines.next() else {
        return HashMap::new();
    };
    let columns: Vec<&str> = header.split('\t').collect();
    let Some(index) = IsotopePairColumns::from_header(&columns) else {
        return HashMap::new();
    };

    let mut pairs = HashMap::new();
    for (line_idx, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let Some(pair) = parse_isotope_pair(&fields, &index, line_idx + 2) else {
            continue;
        };
        pairs.insert(pair.light_scan, pair.clone());
        pairs.insert(pair.heavy_scan, pair);
    }
    pairs
}

#[derive(Debug, Clone, Copy)]
struct IsotopePairColumns {
    light_file: usize,
    light_scan: usize,
    heavy_file: usize,
    heavy_scan: usize,
    rt_light: usize,
    rt_heavy: usize,
    mz_light: usize,
    mz_heavy: usize,
    light_charge: usize,
    heavy_charge: usize,
}

impl IsotopePairColumns {
    fn from_header(columns: &[&str]) -> Option<Self> {
        Some(Self {
            light_file: column(columns, "light_file")?,
            light_scan: column(columns, "light_scan")?,
            heavy_file: column(columns, "heavy_file")?,
            heavy_scan: column(columns, "heavy_scan")?,
            rt_light: column(columns, "rt_light")?,
            rt_heavy: column(columns, "rt_heavy")?,
            mz_light: column(columns, "mz_light")?,
            mz_heavy: column(columns, "mz_heavy")?,
            light_charge: column(columns, "light_charge")?,
            heavy_charge: column(columns, "heavy_charge")?,
        })
    }
}

fn column(columns: &[&str], name: &str) -> Option<usize> {
    columns.iter().position(|col| *col == name)
}

fn parse_isotope_pair(
    fields: &[&str],
    columns: &IsotopePairColumns,
    source_row: usize,
) -> Option<ViewerIsotopePair> {
    let light_scan = parse_u32(fields.get(columns.light_scan)?)?;
    let heavy_scan = parse_u32(fields.get(columns.heavy_scan)?)?;
    Some(ViewerIsotopePair {
        id: format!("{light_scan}:{heavy_scan}"),
        source_artifact: "isotope_pairs.tsv".into(),
        source_row: Some(source_row),
        light_file: nonempty_string(fields.get(columns.light_file)?),
        heavy_file: nonempty_string(fields.get(columns.heavy_file)?),
        light_scan,
        heavy_scan,
        rt_light_min: parse_f64(fields.get(columns.rt_light)?)?,
        rt_heavy_min: parse_f64(fields.get(columns.rt_heavy)?)?,
        mz_light: parse_f64(fields.get(columns.mz_light)?)?,
        mz_heavy: parse_f64(fields.get(columns.mz_heavy)?)?,
        light_charge: parse_u8(fields.get(columns.light_charge)?)?,
        heavy_charge: parse_u8(fields.get(columns.heavy_charge)?)?,
    })
}

fn nonempty_string(value: &&str) -> Option<String> {
    let trimmed = value.trim();
    (!trimmed.is_empty()).then(|| trimmed.to_string())
}

fn parse_u32(value: &&str) -> Option<u32> {
    value.trim().parse().ok()
}

fn parse_u8(value: &&str) -> Option<u8> {
    value.trim().parse().ok()
}

fn parse_f64(value: &&str) -> Option<f64> {
    value.trim().parse().ok()
}

fn unix_timestamp() -> u64 {
    use std::time::{SystemTime, UNIX_EPOCH};
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0)
}

fn iso8601_utc_from_unix(secs: u64) -> String {
    let days = (secs / 86_400) as i64;
    let seconds_of_day = secs % 86_400;
    let hour = seconds_of_day / 3_600;
    let minute = (seconds_of_day % 3_600) / 60;
    let second = seconds_of_day % 60;
    let (year, month, day) = civil_from_days(days);
    format!("{year:04}-{month:02}-{day:02}T{hour:02}:{minute:02}:{second:02}.000Z")
}

fn civil_from_days(days_since_epoch: i64) -> (i32, u32, u32) {
    let z = days_since_epoch + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = mp + if mp < 10 { 3 } else { -9 };
    let year = y + if month <= 2 { 1 } else { 0 };
    (year as i32, month as u32, day as u32)
}

fn crosslink_id(hit: &AnnotatedHit, idx: usize) -> String {
    let scan = hit
        .scan
        .map(|s| s.to_string())
        .unwrap_or_else(|| "noscan".into());
    format!("xl_{idx}_{scan}_{}", hit.job_id)
}

fn protein_pair_key(left: &str, right: &str) -> String {
    if left.trim().is_empty() {
        return right.to_string();
    }
    if right.trim().is_empty() {
        return left.to_string();
    }
    if left <= right {
        format!("{left}|{right}")
    } else {
        format!("{right}|{left}")
    }
}

fn collect_protein_ids(
    hits: &[AnnotatedHit],
    proteins_map: &HashMap<String, String>,
) -> HashSet<String> {
    let mut ids = HashSet::new();
    for hit in hits {
        let plan = None;
        let m = map_crosslink(
            &hit.hit.seq1,
            &hit.hit.seq2,
            &hit.hit.prot1,
            &hit.hit.prot2,
            &hit.hit.xlink_position,
            proteins_map,
            plan,
        );
        if !m.prot1.is_empty() {
            ids.insert(m.prot1);
        }
        if !m.prot2.is_empty() {
            ids.insert(m.prot2);
        }
    }
    ids
}

fn build_protein_list(fasta: &FastaDatabase, wanted: &HashSet<String>) -> Vec<ViewerProtein> {
    let lookup = protein_lookup(fasta);
    let mut seen = HashSet::new();
    let mut out = Vec::new();

    for id in wanted {
        if !seen.insert(id.clone()) {
            continue;
        }
        if let Some(seq) = lookup.get(id) {
            out.push(ViewerProtein {
                id: id.clone(),
                display_name: id.clone(),
                sequence: seq.clone(),
            });
        }
    }
    out.sort_by(|a, b| a.id.cmp(&b.id));
    out
}

fn build_qc(hits: &[AnnotatedHit], stats: &PrefilterStats) -> ViewerQc {
    ViewerQc {
        funnel: vec![
            NamedCount {
                label: "MS/MS scans".into(),
                count: stats.scans_total as f64,
            },
            NamedCount {
                label: "Diagnostic-ion positive".into(),
                count: stats.diagnostic_positive as f64,
            },
            NamedCount {
                label: "Isotope pairs".into(),
                count: stats.isotope_pairs as f64,
            },
            NamedCount {
                label: "Passed to xQuest".into(),
                count: stats.filtered_scans as f64,
            },
        ],
        outcomes: hard_status_counts(hits),
        glycan_top: glycan_top_counts(hits),
        site_dist: glyco_site_counts(hits),
        score_hist: histogram(hits.iter().map(|h| h.hit.score).collect(), 20),
        ppm_hist: histogram(hits.iter().map(|h| h.hit.precursor_error_ppm).collect(), 20),
    }
}

fn hard_status_counts(hits: &[AnnotatedHit]) -> Vec<NamedCount> {
    let mut counts: BTreeMap<&'static str, usize> = BTreeMap::new();
    for hit in hits {
        *counts.entry(hit.hard_status.as_str()).or_insert(0) += 1;
    }
    counts
        .into_iter()
        .map(|(label, count)| NamedCount {
            label: label.to_string(),
            count: count as f64,
        })
        .collect()
}

fn glycan_top_counts(hits: &[AnnotatedHit]) -> Vec<NamedCount> {
    let mut counts: HashMap<String, usize> = HashMap::new();
    for hit in hits
        .iter()
        .filter(|h| h.postfilter_status == PostfilterStatus::Pass)
    {
        if let Some(comp) = &hit.glycan_composition {
            *counts.entry(comp.clone()).or_insert(0) += 1;
        }
    }
    let mut data: Vec<_> = counts
        .into_iter()
        .map(|(label, count)| NamedCount {
            label,
            count: count as f64,
        })
        .collect();
    data.sort_by(|a, b| {
        b.count
            .partial_cmp(&a.count)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.label.cmp(&b.label))
    });
    data.truncate(15);
    data
}

fn glyco_site_counts(hits: &[AnnotatedHit]) -> Vec<NamedCount> {
    let mut counts: BTreeMap<char, usize> = BTreeMap::new();
    for hit in hits
        .iter()
        .filter(|h| h.postfilter_status == PostfilterStatus::Pass)
    {
        for site in &hit.glyco_sites {
            *counts.entry(site.residue).or_insert(0) += 1;
        }
    }
    counts
        .into_iter()
        .map(|(residue, count)| NamedCount {
            label: glyco_residue_label(residue),
            count: count as f64,
        })
        .collect()
}

fn glyco_residue_label(residue: char) -> String {
    match residue {
        'N' => "N (N-linked)".into(),
        'S' => "S (O-linked)".into(),
        'T' => "T (O-linked)".into(),
        other => other.to_string(),
    }
}

fn histogram(values: Vec<f64>, bins: usize) -> Histogram {
    let finite: Vec<f64> = values.into_iter().filter(|v| v.is_finite()).collect();
    if finite.is_empty() {
        return Histogram {
            bins,
            min: 0.0,
            max: 0.0,
            counts: vec![0; bins.max(1)],
            n: 0,
        };
    }
    let min = finite.iter().copied().fold(f64::INFINITY, f64::min);
    let max = finite.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let span = (max - min).max(1e-9);
    let bins = bins.max(1);
    let mut counts = vec![0usize; bins];
    for v in &finite {
        let mut idx = (((v - min) / span) * bins as f64) as usize;
        if idx >= bins {
            idx = bins - 1;
        }
        counts[idx] += 1;
    }
    Histogram {
        bins,
        min,
        max,
        counts,
        n: finite.len(),
    }
}

/// Write `viewer.json` and copy FASTA snapshot to `viewer_dir`.
pub fn write_viewer_bundle(
    viewer_dir: &Path,
    bundle: &ViewerBundle,
    fasta: &FastaDatabase,
) -> Result<(), String> {
    std::fs::create_dir_all(viewer_dir).map_err(|e| e.to_string())?;
    let json = serde_json::to_string_pretty(bundle).map_err(|e| e.to_string())?;
    std::fs::write(viewer_dir.join("viewer.json"), json).map_err(|e| e.to_string())?;
    write_fasta_snapshot(&viewer_dir.join("database.fasta"), fasta)?;
    Ok(())
}

fn write_fasta_snapshot(path: &Path, fasta: &FastaDatabase) -> Result<(), String> {
    let mut content = String::new();
    for entry in &fasta.entries {
        content.push('>');
        content.push_str(&entry.header);
        content.push('\n');
        content.push_str(&entry.sequence);
        content.push('\n');
    }
    std::fs::write(path, content).map_err(|e| e.to_string())
}

/// Copy pre-built static viewer assets from `source_dir` into `viewer_dir`.
pub fn install_viewer_assets(source_dir: &Path, viewer_dir: &Path) -> Result<(), String> {
    if !source_dir.is_dir() {
        return Err(format!(
            "viewer assets not found at {} (run `npm run build` in viewer/)",
            source_dir.display()
        ));
    }
    std::fs::create_dir_all(viewer_dir).map_err(|e| e.to_string())?;
    for entry in std::fs::read_dir(source_dir).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let name = entry.file_name();
        let dest = viewer_dir.join(&name);
        if entry.path().is_file() {
            std::fs::copy(entry.path(), &dest).map_err(|e| e.to_string())?;
        }
    }
    let serve_script =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("viewer/serve-viewer.sh");
    if serve_script.is_file() {
        let dest = viewer_dir.join("serve-viewer.sh");
        std::fs::copy(&serve_script, &dest).map_err(|e| e.to_string())?;
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            std::fs::set_permissions(&dest, std::fs::Permissions::from_mode(0o755))
                .map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}

pub fn default_viewer_assets_dir() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("viewer/dist")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::{FastaDatabase, FastaEntry};
    use crate::results::extract::XQuestHit;
    use crate::results::filter::HardStatus;
    use std::path::PathBuf;

    fn fasta() -> FastaDatabase {
        FastaDatabase {
            path: PathBuf::from("test.fasta"),
            entries: vec![
                FastaEntry {
                    header: "sp|P00761|TRYP_PIG trypsin".into(),
                    sequence: "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHR".into(),
                },
                FastaEntry {
                    header: "HRP".into(),
                    sequence: "QLTPTFYDNSCPNVSNIVRDTIVNELR".into(),
                },
            ],
        }
    }

    fn pass_hit() -> AnnotatedHit {
        pass_hit_with_scan(Some(101))
    }

    fn pass_hit_with_scan(scan: Option<u32>) -> AnnotatedHit {
        AnnotatedHit {
            hit: XQuestHit {
                seq1: "DTHK".into(),
                seq2: "DTIVNELR".into(),
                prot1: "P00761".into(),
                prot2: "HRP".into(),
                xlink_position: "2-3".into(),
                charge: 4,
                score: 12.5,
                precursor_error_ppm: 2.1,
                topology: "inter".into(),
                xquest_version: Some("xquest 2.1.7".into()),
                xlinkions_matched: Some("1/10".into()),
                backboneions_matched: Some("3/20".into()),
                matched_ions: vec![crate::results::extract::XQuestMatchedIon {
                    ion_type: "p1_b1".into(),
                    position: Some("1".into()),
                    theoretical_mz: 100.1,
                    observed_mz: 100.1,
                    delta_mz: Some(0.0),
                    delta_ppm: Some(0.0),
                    intensity: Some(10.0),
                    label: Some("p1_b1".into()),
                }],
                ..Default::default()
            },
            job_id: "HexNAc_1_".into(),
            source_file: Some(PathBuf::from("run.mzXML")),
            scan,
            glycan_name: Some("HexNAc(1)".into()),
            glycan_composition: Some("HexNAc(1)".into()),
            glycan_mass: Some(203.079),
            loss_label: Some("none".into()),
            glyco_residue: Some('N'),
            glyco_peptide: Some(2),
            glyco_sites: vec![super::super::filter::GlycoSite {
                peptide: 2,
                peptide_position: 1,
                residue: 'N',
                sequon_present: Some(true),
                plausible: true,
            }],
            all_sites_plausible: true,
            n_glycan_pseudo: 1,
            matched_families: vec!["HexNAc".into()],
            matched_ions: vec![crate::prefilter::MatchedIon {
                family: "HexNAc".into(),
                expected_mz: 366.139472,
                observed_mz: 366.1400,
                loss_label: String::new(),
                peak_index: 8,
                intensity: 125_000.0,
                error_ppm: 1.442,
            }],
            matched_ion_count: 1,
            sequon_present: Some(true),
            charge_plausible: true,
            hard_status: HardStatus::Pass,
            soft_score: 13.6,
            postfilter_status: PostfilterStatus::Pass,
        }
    }

    fn monolink_hit() -> AnnotatedHit {
        let mut hit = pass_hit_with_scan(Some(101));
        hit.hit.link_type = "monolink".into();
        hit.hit.seq1 = "DTHK".into();
        hit.hit.seq2 = "".into();
        hit.hit.prot1 = "P00761".into();
        hit.hit.prot2 = "".into();
        hit.hit.xlink_position = "2".into();
        hit.hit.xlinker_mass = Some(156.07864);
        hit.hit.topology = "monolink".into();
        hit.glyco_peptide = Some(1);
        hit.glyco_sites[0].peptide = 1;
        hit
    }

    #[test]
    fn load_isotope_pairs_maps_both_partner_scans() {
        let dir =
            std::env::temp_dir().join(format!("gq_viewer_isotope_pairs_{}", std::process::id()));
        let spectra_dir = dir.join("spectra");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&spectra_dir).unwrap();
        std::fs::write(
            dir.join("isotope_pairs.tsv"),
            "light_file\tlight_scan\theavy_file\theavy_scan\trt_light\trt_heavy\tmz_light\tmz_heavy\tlight_charge\theavy_charge\n\
             run.mzXML\t6715\trun.mzXML\t6626\t19.7142\t19.5497\t938.4606\t941.4794\t4\t4\n",
        )
        .unwrap();

        let pairs = load_isotope_pairs_for_viewer(&spectra_dir);

        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs.get(&6715).unwrap().heavy_scan, 6626);
        assert_eq!(pairs.get(&6626).unwrap().light_scan, 6715);
        assert!((pairs.get(&6715).unwrap().mz_heavy - 941.4794).abs() < 0.0001);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn build_viewer_bundle_generates_unique_ids_for_same_scan_and_job() {
        let ctx = ReportContext {
            project: "test".into(),
            input_label: "run.mzXML".into(),
            crosslinker_name: "dss".into(),
            xlink_sites: "K:K".into(),
            glycan_library: "nglyc309".into(),
        };
        let stats = PrefilterStats::default();
        let crosslinker = crate::crosslinker::CrosslinkerProfile::resolve(
            &crate::cli::settings::Settings::defaults(),
            Some("dss"),
        )
        .unwrap();

        let first = pass_hit_with_scan(Some(23363));
        let mut second = pass_hit_with_scan(Some(23363));
        second.hit.seq1 = "SEIAHR".into();
        second.hit.xlink_position = "2-3".into();

        let bundle = build_viewer_bundle(
            &[first, second],
            &stats,
            None,
            &fasta(),
            &crosslinker,
            &crate::cli::settings::Settings::defaults(),
            &ctx,
            Path::new("/nonexistent"),
        );

        assert_eq!(bundle.crosslinks.len(), 2);
        assert_ne!(bundle.crosslinks[0].id, bundle.crosslinks[1].id);
        for xl in &bundle.crosslinks {
            assert!(bundle.filtering.contains_key(&xl.id));
        }
    }

    #[test]
    fn build_viewer_bundle_has_schema_version() {
        let ctx = ReportContext {
            project: "test".into(),
            input_label: "run.mzXML".into(),
            crosslinker_name: "dss".into(),
            xlink_sites: "K:K".into(),
            glycan_library: "nglyc309".into(),
        };
        let stats = PrefilterStats::default();
        let crosslinker = crate::crosslinker::CrosslinkerProfile::resolve(
            &crate::cli::settings::Settings::defaults(),
            Some("dss"),
        )
        .unwrap();
        let bundle = build_viewer_bundle(
            &[pass_hit()],
            &stats,
            None,
            &fasta(),
            &crosslinker,
            &crate::cli::settings::Settings::defaults(),
            &ctx,
            Path::new("/nonexistent"),
        );
        assert_eq!(bundle.viewer_schema_version, VIEWER_SCHEMA_VERSION);
        assert_eq!(bundle.viewer_schema_version, 4);
        assert!(!bundle.meta.generated_at_iso.is_empty());
        assert!(bundle.meta.generated_at_iso.contains('T'));
        assert!(bundle.meta.generated_at_unix.is_some());
        assert_eq!(bundle.crosslinks.len(), 1);
        assert!(bundle.crosslinks[0].mapped);
        assert_eq!(bundle.crosslinks[0].glyco_sites.len(), 1);
        assert_eq!(
            bundle.crosslinks[0].source_file.as_deref(),
            Some("run.mzXML")
        );
        assert_eq!(bundle.crosslinks[0].protein_pair_key, "HRP|P00761");
        assert_eq!(bundle.proteins.len(), 2);
        let filtering = bundle
            .filtering
            .get(&bundle.crosslinks[0].id)
            .expect("schema v4 must emit Filtering data for each crosslink");
        assert_eq!(filtering.input_scan.scan, Some(101));
        assert_eq!(filtering.diagnostic_prefilter.status, "matched");
        assert_eq!(filtering.postfilter.status, "pass");
        let site_rule = filtering
            .postfilter
            .rules
            .iter()
            .find(|rule| rule.name == "plausible_glycan_sites")
            .unwrap();
        assert_eq!(site_rule.status, "soft");
        assert_eq!(site_rule.value, "1/1");
        assert_eq!(filtering.xquest_search.status, "matched");
    }

    #[test]
    fn qc_site_distribution_counts_every_glycan_site() {
        let ctx = ReportContext {
            project: "test".into(),
            input_label: "run.mzXML".into(),
            crosslinker_name: "dss".into(),
            xlink_sites: "K:K".into(),
            glycan_library: "nglyc309".into(),
        };
        let settings = crate::cli::settings::Settings::defaults();
        let crosslinker =
            crate::crosslinker::CrosslinkerProfile::resolve(&settings, Some("dss")).unwrap();
        let mut hit = pass_hit();
        hit.glyco_sites.push(super::super::filter::GlycoSite {
            peptide: 1,
            peptide_position: 2,
            residue: 'S',
            sequon_present: None,
            plausible: true,
        });
        hit.n_glycan_pseudo = 2;

        let bundle = build_viewer_bundle(
            &[hit],
            &PrefilterStats::default(),
            None,
            &fasta(),
            &crosslinker,
            &settings,
            &ctx,
            Path::new("/nonexistent"),
        );

        assert_eq!(
            bundle.qc.site_dist.iter().map(|row| row.count).sum::<f64>(),
            2.0
        );
    }

    #[test]
    fn build_viewer_bundle_maps_monolink_single_endpoint_and_crosslinker_mass() {
        let ctx = ReportContext {
            project: "test".into(),
            input_label: "run.mzXML".into(),
            crosslinker_name: "dss".into(),
            xlink_sites: "K:K".into(),
            glycan_library: "nglyc309".into(),
        };
        let stats = PrefilterStats::default();
        let crosslinker = crate::crosslinker::CrosslinkerProfile::resolve(
            &crate::cli::settings::Settings::defaults(),
            Some("dss"),
        )
        .unwrap();

        let bundle = build_viewer_bundle(
            &[monolink_hit()],
            &stats,
            None,
            &fasta(),
            &crosslinker,
            &crate::cli::settings::Settings::defaults(),
            &ctx,
            Path::new("/nonexistent"),
        );

        assert_eq!(bundle.meta.crosslinker_mw, crosslinker.xlinkermw);
        assert_eq!(bundle.crosslinks.len(), 1);
        let xl = &bundle.crosslinks[0];
        assert_eq!(xl.link_type, "monolink");
        assert_eq!(xl.protein1, "P00761");
        assert_eq!(xl.protein2, "");
        assert_eq!(xl.abs_pos1, Some(26));
        assert_eq!(xl.abs_pos2, None);
        assert_eq!(xl.mapped, true);
        assert_eq!(xl.protein_pair_key, "P00761");
        assert!((xl.xlinker_mass.unwrap() - 156.07864).abs() < 1e-6);
        assert_eq!(bundle.proteins.len(), 1);
        assert_eq!(bundle.proteins[0].id, "P00761");
        let filtering = bundle.filtering.get(&xl.id).unwrap();
        assert_eq!(filtering.input_scan.scan, Some(101));
        assert_eq!(filtering.postfilter.status, "pass");
    }

    #[test]
    fn build_viewer_bundle_carries_retention_time_to_spectra_and_crosslinks() {
        let ctx = ReportContext {
            project: "test".into(),
            input_label: "run.mzXML".into(),
            crosslinker_name: "dss".into(),
            xlink_sites: "K:K".into(),
            glycan_library: "nglyc309".into(),
        };
        let stats = PrefilterStats::default();
        let crosslinker = crate::crosslinker::CrosslinkerProfile::resolve(
            &crate::cli::settings::Settings::defaults(),
            Some("dss"),
        )
        .unwrap();
        let spectra_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/mzxml");
        let bundle = build_viewer_bundle(
            &[pass_hit_with_scan(Some(42))],
            &stats,
            None,
            &fasta(),
            &crosslinker,
            &crate::cli::settings::Settings::defaults(),
            &ctx,
            &spectra_dir,
        );
        let spectrum = bundle.spectra.get("42").expect("scan 42 should be bundled");
        assert!((spectrum.retention_time_min - 20.0).abs() < 0.001);
        assert!((bundle.crosslinks[0].retention_time_min.unwrap() - 20.0).abs() < 0.001);
    }

    #[test]
    fn write_viewer_bundle_roundtrip() {
        let dir = std::env::temp_dir().join(format!("gq_viewer_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let ctx = ReportContext {
            project: "test".into(),
            input_label: "run.mzXML".into(),
            crosslinker_name: "dss".into(),
            xlink_sites: "K:K".into(),
            glycan_library: "nglyc309".into(),
        };
        let crosslinker = crate::crosslinker::CrosslinkerProfile::resolve(
            &crate::cli::settings::Settings::defaults(),
            Some("dss"),
        )
        .unwrap();
        let bundle = build_viewer_bundle(
            &[pass_hit()],
            &PrefilterStats::default(),
            None,
            &fasta(),
            &crosslinker,
            &crate::cli::settings::Settings::defaults(),
            &ctx,
            Path::new("/nonexistent"),
        );
        write_viewer_bundle(&dir, &bundle, &fasta()).unwrap();
        assert!(dir.join("viewer.json").is_file());
        assert!(dir.join("database.fasta").is_file());
        let json = std::fs::read_to_string(dir.join("viewer.json")).unwrap();
        assert!(json.contains("viewer_schema_version"));
        let _ = std::fs::remove_dir_all(&dir);
    }
}
