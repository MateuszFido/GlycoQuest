//! Write `plan.json` summarizing a GlycoQuest run.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::crosslinker::CrosslinkerProfile;
use crate::jobs::JobPlan;
use crate::prefilter::PrefilterResult;
use crate::xquest::GeneratedJob;

#[derive(Debug, Clone, PartialEq)]
pub struct RunPlanDocument {
    pub isotope_prefilter_enabled: bool,
    pub crosslinker_name: String,
    pub crosslinker_label: String,
    pub job_count: usize,
    pub total_comparisons: u64,
    pub prefilter_stats: PrefilterResult,
    pub pruned_mzxml_paths: Vec<PathBuf>,
    pub jobs: Vec<GeneratedJob>,
    pub commands: Vec<String>,
}

pub fn write_plan_json(out_dir: &Path, doc: &RunPlanDocument) -> Result<PathBuf, String> {
    let path = out_dir.join("plan.json");
    let mut file = File::create(&path).map_err(|err| err.to_string())?;
    file.write_all(render_json(doc).as_bytes())
        .map_err(|err| err.to_string())?;
    Ok(path)
}

fn render_json(doc: &RunPlanDocument) -> String {
    let stats = &doc.prefilter_stats.stats;
    let jobs: Vec<String> = doc
        .jobs
        .iter()
        .map(|job| {
            format!(
                r#"    {{
      "job_id": {},
      "directory": {},
      "command": {},
      "run_script": {}
    }}"#,
                json_str(&job.job_id),
                json_str(&job.directory.display().to_string()),
                json_str(&job.command),
                json_str(&job.run_script.display().to_string()),
            )
        })
        .collect();

    let pruned: Vec<String> = doc
        .pruned_mzxml_paths
        .iter()
        .map(|p| json_str(&p.display().to_string()))
        .collect();

    format!(
        r#"{{
  "isotope_prefilter_enabled": {},
  "crosslinker_name": {},
  "crosslinker_label": {},
  "job_count": {},
  "total_comparisons": {},
  "prefilter": {{
    "scans_total": {},
    "diagnostic_positive": {},
    "isotope_pairs": {},
    "filtered_scans": {},
    "rejected": {}
  }},
  "spectra": [
{}
  ],
  "jobs": [
{}
  ],
  "commands": [
{}
  ]
}}
"#,
        json_bool(doc.isotope_prefilter_enabled),
        json_str(&doc.crosslinker_name),
        json_str(&doc.crosslinker_label),
        doc.job_count,
        doc.total_comparisons,
        stats.scans_total,
        stats.diagnostic_positive,
        stats.isotope_pairs,
        stats.filtered_scans,
        stats.rejected,
        pruned.join(",\n"),
        jobs.join(",\n"),
        doc.commands
            .iter()
            .map(|c| format!("    {}", json_str(c)))
            .collect::<Vec<_>>()
            .join(",\n"),
    )
}

fn json_str(value: &str) -> String {
    let escaped = value.replace('\\', "\\\\").replace('"', "\\\"");
    format!("\"{escaped}\"")
}

fn json_bool(value: bool) -> &'static str {
    if value {
        "true"
    } else {
        "false"
    }
}

pub fn build_run_plan_document(
    crosslinker: &CrosslinkerProfile,
    prefilter: &PrefilterResult,
    job_plan: &JobPlan,
    pruned_mzxml_paths: Vec<PathBuf>,
    jobs: Vec<GeneratedJob>,
) -> RunPlanDocument {
    let commands = jobs.iter().map(|job| job.command.clone()).collect();
    RunPlanDocument {
        isotope_prefilter_enabled: crosslinker.requires_isotope_pair_prefilter(),
        crosslinker_name: crosslinker.name.clone(),
        crosslinker_label: crosslinker.label.as_str().to_string(),
        job_count: job_plan.jobs.len(),
        total_comparisons: job_plan.total_comparisons,
        prefilter_stats: prefilter.clone(),
        pruned_mzxml_paths,
        jobs,
        commands,
    }
}
