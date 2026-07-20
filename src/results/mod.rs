// Copyright (c) ETH Zurich, Mateusz Fido

//! Consolidated GlycoQuest xQuest result handling.

mod extract;
mod filter;
mod mapping;
mod report;
mod spectrum;
mod viewer;

pub use extract::{XQuestHit, extract_hits_from_xml, find_result_xmls};
pub use filter::{AnnotatedHit, PostfilterStatus, apply_postfilters, write_annotated_csv};
pub use report::ReportContext;

use std::path::Path;

use crate::cli::settings::Settings;
use crate::crosslinker::CrosslinkerProfile;
use crate::fasta::FastaDatabase;
use crate::jobs::JobManifest;
use crate::prefilter::PrefilterResult;
use crate::xquest::JobRunRecord;

/// Summary written after consolidating xQuest XML hits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ConsolidationSummary {
    pub result_xmls: usize,
    pub hits: usize,
    /// Passing crosslinks written to `results/xiview.csv`.
    pub xiview_rows: usize,
    /// Of those, the number with both peptides mapped to a protein (absolute positions resolved).
    pub xiview_mapped: usize,
    /// Crosslinks written to `results/viewer/viewer.json`.
    pub viewer_crosslinks: usize,
}

/// Write `results/failed_jobs.tsv` listing jobs whose `run.sh` did not succeed.
pub fn write_failed_jobs_tsv(path: &Path, records: &[JobRunRecord]) -> Result<(), String> {
    let mut lines = vec!["job\tlog\treason".to_string()];
    for record in records.iter().filter(|record| !record.success) {
        lines.push(format!(
            "{}\t{}\t{}",
            record.job_name,
            record.log_path.display(),
            record.error.as_deref().unwrap_or("unknown"),
        ));
    }
    std::fs::write(path, lines.join("\n") + "\n").map_err(|err| err.to_string())
}

pub fn consolidate_results(
    layout: &crate::output::ProjectLayout,
    settings: &Settings,
    crosslinker: &CrosslinkerProfile,
    prefilter: &PrefilterResult,
    manifest: Option<&JobManifest>,
    fasta: &FastaDatabase,
    report_ctx: &ReportContext,
) -> Result<ConsolidationSummary, String> {
    let results_dir = layout.results_dir();
    std::fs::create_dir_all(&results_dir).map_err(|err| err.to_string())?;

    let xmls = find_result_xmls(&layout.jobs_dir())?;
    let mut all_hits = Vec::new();
    for (job_id, xml) in &xmls {
        for hit in extract_hits_from_xml(xml)? {
            all_hits.push((job_id.clone(), hit));
        }
    }

    let annotated = apply_postfilters(all_hits, settings, crosslinker, prefilter, manifest);
    let deduped = deduplicate(annotated);
    let hits = deduped.len();
    write_annotated_csv(&results_dir.join("glycoquest_xquest.csv"), &deduped)?;

    let xiview =
        report::write_xiview_csv(&results_dir.join("xiview.csv"), &deduped, manifest, fasta)?;
    report::write_html_report(
        &results_dir.join("report.html"),
        &deduped,
        &prefilter.stats,
        report_ctx,
    )?;

    let viewer_dir = layout.viewer_dir();
    let bundle = viewer::build_viewer_bundle(
        &deduped,
        &prefilter.stats,
        manifest,
        fasta,
        crosslinker,
        settings,
        report_ctx,
        &layout.spectra_dir(),
    );
    viewer::write_viewer_bundle(&viewer_dir, &bundle, fasta)?;
    if let Err(err) =
        viewer::install_viewer_assets(&viewer::default_viewer_assets_dir(), &viewer_dir)
    {
        eprintln!("run: viewer static assets not installed: {err}");
    }

    Ok(ConsolidationSummary {
        result_xmls: xmls.len(),
        hits,
        xiview_rows: xiview.rows,
        xiview_mapped: xiview.mapped,
        viewer_crosslinks: bundle.crosslinks.len(),
    })
}

/// Collapse hits to the best identification per source scan. One MS/MS scan
/// should contribute at most one passing crosslink; when a scan has no passing
/// row, keep the strongest failed row for audit mode. Result is sorted by soft
/// score descending.
fn deduplicate(annotated: Vec<AnnotatedHit>) -> Vec<AnnotatedHit> {
    use std::collections::HashMap;

    let mut best: HashMap<String, AnnotatedHit> = HashMap::new();
    for hit in annotated {
        let key = dedup_key(&hit);
        match best.get(&key) {
            Some(existing) if !is_better_dedup_hit(&hit, existing) => {}
            _ => {
                best.insert(key, hit);
            }
        }
    }

    let mut rows: Vec<AnnotatedHit> = best.into_values().collect();
    rows.sort_by(|a, b| {
        b.soft_score
            .partial_cmp(&a.soft_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    rows
}

fn dedup_key(hit: &AnnotatedHit) -> String {
    match hit.scan {
        Some(scan) => format!(
            "scan|{}|{}",
            hit.source_file
                .as_ref()
                .map(|path| path.display().to_string())
                .unwrap_or_default(),
            scan
        ),
        None => format!(
            "noscan|{}|{}|{}|{}|{}",
            hit.hit.seq1,
            hit.hit.seq2,
            hit.hit.xlink_position,
            hit.glycan_composition.clone().unwrap_or_default(),
            hit.loss_label.clone().unwrap_or_default(),
        ),
    }
}

fn is_better_dedup_hit(candidate: &AnnotatedHit, existing: &AnnotatedHit) -> bool {
    let candidate_pass = candidate.postfilter_status == PostfilterStatus::Pass;
    let existing_pass = existing.postfilter_status == PostfilterStatus::Pass;
    if candidate_pass != existing_pass {
        return candidate_pass;
    }

    candidate
        .soft_score
        .partial_cmp(&existing.soft_score)
        .unwrap_or(std::cmp::Ordering::Equal)
        .then_with(|| {
            candidate
                .hit
                .score
                .partial_cmp(&existing.hit.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .is_gt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::results::extract::XQuestHit;
    use crate::results::filter::{HardStatus, PostfilterStatus};

    fn annotated_hit(
        scan: u32,
        seq1: &str,
        xlink_position: &str,
        score: f64,
        soft_score: f64,
        status: PostfilterStatus,
    ) -> AnnotatedHit {
        AnnotatedHit {
            hit: XQuestHit {
                seq1: seq1.into(),
                seq2: "PEPTIDEK".into(),
                prot1: "P1".into(),
                prot2: "P2".into(),
                xlink_position: xlink_position.into(),
                score,
                ..Default::default()
            },
            job_id: "HexNAc_1_".into(),
            source_file: None,
            scan: Some(scan),
            glycan_name: Some("HexNAc(1)".into()),
            glycan_composition: Some("HexNAc(1)".into()),
            glycan_mass: Some(203.079),
            loss_label: Some("none".into()),
            glyco_residue: Some('N'),
            glyco_peptide: Some(1),
            glyco_sites: vec![crate::results::filter::GlycoSite {
                peptide: 1,
                peptide_position: 1,
                residue: 'N',
                sequon_present: Some(true),
                plausible: true,
            }],
            all_sites_plausible: true,
            n_glycan_pseudo: 1,
            matched_families: vec!["HexNAc".into()],
            matched_ions: vec![],
            matched_ion_count: 1,
            sequon_present: Some(true),
            charge_plausible: true,
            hard_status: if status == PostfilterStatus::Pass {
                HardStatus::Pass
            } else {
                HardStatus::FailScore
            },
            soft_score,
            postfilter_status: status,
        }
    }

    #[test]
    fn deduplicate_keeps_one_best_passing_crosslink_per_scan() {
        let rows = deduplicate(vec![
            annotated_hit(
                32655,
                "AKVFKDVFLEUXIPYSVVR",
                "2,2",
                16.77,
                16.38,
                PostfilterStatus::Pass,
            ),
            annotated_hit(
                32655,
                "AKVFKDVFLEUXIPYSVVR",
                "5,2",
                21.07,
                20.68,
                PostfilterStatus::Pass,
            ),
            annotated_hit(
                32655,
                "AKVFKDVFLEUXIPYSVVR",
                "7,2",
                99.0,
                99.0,
                PostfilterStatus::Fail,
            ),
        ]);

        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].scan, Some(32655));
        assert_eq!(rows[0].hit.xlink_position, "5,2");
        assert_eq!(rows[0].postfilter_status, PostfilterStatus::Pass);
    }

    #[test]
    fn deduplicate_keeps_best_failed_hit_when_scan_has_no_passing_hit() {
        let rows = deduplicate(vec![
            annotated_hit(100, "AAAAK", "2,2", 5.0, 5.0, PostfilterStatus::Fail),
            annotated_hit(100, "BBBBK", "3,2", 7.0, 7.0, PostfilterStatus::Fail),
        ]);

        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].hit.seq1, "BBBBK");
        assert_eq!(rows[0].postfilter_status, PostfilterStatus::Fail);
    }
}
