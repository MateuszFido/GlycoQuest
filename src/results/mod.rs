//! Consolidated GlycoQuest xQuest result handling.

mod extract;
mod filter;

pub use extract::{extract_hits_from_xml, find_result_xmls, write_hits_csv, XQuestHit};
pub use filter::{
    apply_postfilters, write_annotated_csv, AnnotatedHit, PostfilterStatus,
};

use std::path::Path;

use crate::cli::settings::Settings;
use crate::crosslinker::CrosslinkerProfile;
use crate::prefilter::PrefilterResult;

pub fn consolidate_results(
    layout: &crate::output::ProjectLayout,
    settings: &Settings,
    crosslinker: &CrosslinkerProfile,
    prefilter: &PrefilterResult,
) -> Result<(), String> {
    let results_dir = layout.results_dir();
    std::fs::create_dir_all(&results_dir).map_err(|err| err.to_string())?;

    let mut all_hits = Vec::new();
    for xml in find_result_xmls(&layout.jobs_dir())? {
        all_hits.extend(extract_hits_from_xml(&xml)?);
    }

    let annotated = apply_postfilters(all_hits, settings, crosslinker, prefilter);
    write_annotated_csv(&results_dir.join("glycoquest_xquest.csv"), &annotated)
}
