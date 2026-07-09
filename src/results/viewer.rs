//! Interactive viewer bundle (`results/viewer/viewer.json`) for the GlycoQuest CLMS viewer.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::Path;

use serde::Serialize;

use crate::crosslinker::CrosslinkerProfile;
use crate::fasta::FastaDatabase;
use crate::jobs::JobManifest;
use crate::prefilter::PrefilterStats;

use super::filter::{AnnotatedHit, PostfilterStatus};
use super::mapping::{map_crosslink, protein_lookup};
use super::report::ReportContext;
use super::spectrum::{annotate_fragments, load_spectra_for_scans, scans_from_hits};

pub const VIEWER_SCHEMA_VERSION: u32 = 2;

#[derive(Debug, Clone, Serialize)]
pub struct ViewerBundle {
    pub viewer_schema_version: u32,
    pub meta: ViewerMeta,
    pub proteins: Vec<ViewerProtein>,
    pub crosslinks: Vec<ViewerCrosslink>,
    pub qc: ViewerQc,
    pub spectra: HashMap<String, ViewerSpectrum>,
    pub fragments: HashMap<String, ViewerFragments>,
    pub mirror_fragments: HashMap<String, ViewerMirrorFragments>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ViewerMeta {
    pub project: String,
    pub input_label: String,
    pub crosslinker: String,
    pub xlink_sites: String,
    pub glycan_library: String,
    pub resume: bool,
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
    pub topology: String,
    pub protein_pair_key: String,
    pub glycan_name: Option<String>,
    pub glycan_composition: Option<String>,
    pub glyco_residue: Option<char>,
    pub glyco_peptide: Option<u8>,
    pub loss_label: Option<String>,
    pub postfilter_status: String,
    pub mapped: bool,
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
pub struct ViewerFragments {
    pub theoretical_mz: Vec<f32>,
    pub labels: Vec<String>,
    pub matched_indices: Vec<usize>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ViewerMirrorFragments {
    pub theoretical_mz: Vec<f32>,
    pub theoretical_intensity: Vec<f32>,
    pub experimental_mz: Vec<f32>,
    pub experimental_intensity: Vec<f32>,
    pub ion_types: Vec<String>,
    pub labels: Vec<String>,
    pub matched_indices_experimental: Vec<usize>,
    pub matched_indices_theoretical: Vec<usize>,
    pub annotation_source: String,
}

/// Build the viewer JSON bundle from consolidated run artifacts.
pub fn build_viewer_bundle(
    hits: &[AnnotatedHit],
    stats: &PrefilterStats,
    manifest: Option<&JobManifest>,
    fasta: &FastaDatabase,
    crosslinker: &CrosslinkerProfile,
    ctx: &ReportContext,
    spectra_dir: &Path,
) -> ViewerBundle {
    let proteins_map = protein_lookup(fasta);
    let protein_ids = collect_protein_ids(hits, &proteins_map);
    let proteins = build_protein_list(fasta, &protein_ids);

    let scans_needed = scans_from_hits(hits, false);
    let spectra_loaded = load_spectra_for_scans(spectra_dir, &scans_needed);

    let mut crosslinks = Vec::with_capacity(hits.len());
    let mut fragments = HashMap::new();
    let mut mirror_fragments = HashMap::new();

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
        let mapped = mapping.abs1.is_some() && mapping.abs2.is_some();

        let scan_spectrum = hit.scan.and_then(|s| spectra_loaded.get(&s));
        let frag = annotate_fragments(&mapping, crosslinker, scan_spectrum);
        fragments.insert(
            id.clone(),
            ViewerFragments {
                theoretical_mz: frag.theoretical_mz.clone(),
                labels: frag.labels.clone(),
                matched_indices: frag.matched_indices.clone(),
            },
        );
        mirror_fragments.insert(
            id.clone(),
            ViewerMirrorFragments {
                theoretical_mz: frag.theoretical_mz,
                theoretical_intensity: frag.theoretical_intensity,
                experimental_mz: frag.experimental_mz,
                experimental_intensity: frag.experimental_intensity,
                ion_types: frag.ion_types,
                labels: frag.labels,
                matched_indices_experimental: frag.matched_indices_experimental,
                matched_indices_theoretical: frag.matched_indices_theoretical,
                annotation_source: "glycoquest_approx".into(),
            },
        );
        let protein_pair_key = protein_pair_key(&mapping.prot1, &mapping.prot2);

        crosslinks.push(ViewerCrosslink {
            id,
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
            topology: hit.hit.topology.clone(),
            protein_pair_key,
            glycan_name: hit.glycan_name.clone(),
            glycan_composition: hit.glycan_composition.clone(),
            glyco_residue: hit.glyco_residue,
            glyco_peptide: hit.glyco_peptide,
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

    let generated_at_unix = unix_timestamp();
    let generated_at_iso = iso8601_utc_from_unix(generated_at_unix);

    ViewerBundle {
        viewer_schema_version: VIEWER_SCHEMA_VERSION,
        meta: ViewerMeta {
            project: ctx.project.clone(),
            input_label: ctx.input_label.clone(),
            crosslinker: ctx.crosslinker_name.clone(),
            xlink_sites: ctx.xlink_sites.clone(),
            glycan_library: ctx.glycan_library.clone(),
            resume: ctx.resume,
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
        fragments,
        mirror_fragments,
    }
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
    format!(
        "xl_{}_{}",
        hit.scan
            .map(|s| s.to_string())
            .unwrap_or_else(|| idx.to_string()),
        hit.job_id
    )
}

fn protein_pair_key(left: &str, right: &str) -> String {
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
        ids.insert(m.prot1);
        ids.insert(m.prot2);
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
        if let Some(residue) = hit.glyco_residue {
            *counts.entry(residue).or_insert(0) += 1;
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
            n_glycan_pseudo: 1,
            matched_families: vec!["HexNAc".into()],
            matched_ion_count: 3,
            sequon_present: Some(true),
            charge_plausible: true,
            hard_status: HardStatus::Pass,
            soft_score: 13.6,
            postfilter_status: PostfilterStatus::Pass,
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
            resume: false,
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
            &ctx,
            Path::new("/nonexistent"),
        );
        assert_eq!(bundle.viewer_schema_version, VIEWER_SCHEMA_VERSION);
        assert_eq!(bundle.viewer_schema_version, 2);
        assert!(!bundle.meta.generated_at_iso.is_empty());
        assert!(bundle.meta.generated_at_iso.contains('T'));
        assert!(bundle.meta.generated_at_unix.is_some());
        assert_eq!(bundle.crosslinks.len(), 1);
        assert!(bundle.crosslinks[0].mapped);
        assert_eq!(
            bundle.crosslinks[0].source_file.as_deref(),
            Some("run.mzXML")
        );
        assert_eq!(bundle.crosslinks[0].protein_pair_key, "HRP|P00761");
        assert_eq!(bundle.proteins.len(), 2);
        assert!(bundle.fragments.contains_key(&bundle.crosslinks[0].id));
        let mirror = bundle
            .mirror_fragments
            .get(&bundle.crosslinks[0].id)
            .expect("mirror fragments should be emitted");
        assert_eq!(mirror.annotation_source, "glycoquest_approx");
        assert_eq!(
            mirror.theoretical_mz.len(),
            mirror.theoretical_intensity.len()
        );
    }

    #[test]
    fn build_viewer_bundle_carries_retention_time_to_spectra_and_crosslinks() {
        let ctx = ReportContext {
            project: "test".into(),
            input_label: "run.mzXML".into(),
            crosslinker_name: "dss".into(),
            xlink_sites: "K:K".into(),
            glycan_library: "nglyc309".into(),
            resume: false,
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
            resume: false,
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
