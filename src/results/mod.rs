//! Consolidated GlycoQuest xQuest result handling.

mod extract;
mod filter;
mod mapping;
mod report;
mod spectrum;
mod viewer;

pub use extract::{extract_hits_from_xml, find_result_xmls, XQuestHit};
pub use filter::{apply_postfilters, write_annotated_csv, AnnotatedHit};
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

    let xiview = report::write_xiview_csv(
        &results_dir.join("xiview.csv"),
        &deduped,
        manifest,
        fasta,
    )?;
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
        report_ctx,
        &layout.spectra_dir(),
    );
    viewer::write_viewer_bundle(&viewer_dir, &bundle, fasta)?;
    if let Err(err) = viewer::install_viewer_assets(&viewer::default_viewer_assets_dir(), &viewer_dir)
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

/// Collapse hits that describe the same scan/peptide-pair/glycan, keeping the
/// highest-scoring copy. xQuest reports the same peptide pair once per glycan
/// job, so without this the CSV is heavily duplicated. Result is sorted by
/// soft score descending.
fn deduplicate(annotated: Vec<AnnotatedHit>) -> Vec<AnnotatedHit> {
    use std::collections::HashMap;

    let mut best: HashMap<String, AnnotatedHit> = HashMap::new();
    for hit in annotated {
        let key = format!(
            "{}|{}|{}|{}|{}|{}",
            hit.scan.map(|s| s.to_string()).unwrap_or_default(),
            hit.hit.seq1,
            hit.hit.seq2,
            hit.hit.xlink_position,
            hit.glycan_composition.clone().unwrap_or_default(),
            hit.loss_label.clone().unwrap_or_default(),
        );
        match best.get(&key) {
            Some(existing) if existing.hit.score >= hit.hit.score => {}
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
