//! Write xQuest job directories under `jobs/`.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::cli::settings::Settings;
use crate::crosslinker::CrosslinkerProfile;
use crate::fasta::FastaDatabase;
use crate::jobs::{JobPlan, PlannedJob};
use crate::prefilter::PrefilterResult;
use crate::xquest::defs::write_job_defs;
use crate::xquest::matchlist::{build_matchlist, isotopepairs_path, specxml_filename, write_matchlist};
use crate::xquest::XQuestRuntime;

#[derive(Debug, Clone, PartialEq)]
pub struct GeneratedJob {
    pub job_id: String,
    pub directory: PathBuf,
    pub run_script: PathBuf,
    pub matchlist_path: PathBuf,
    pub command: String,
}

pub fn generate_jobs(
    layout: &crate::output::ProjectLayout,
    plan: &JobPlan,
    prefilter: &PrefilterResult,
    pruned_mzxml_paths: &[PathBuf],
    crosslinker: &CrosslinkerProfile,
    settings: &Settings,
    fasta: &FastaDatabase,
    runtime: &XQuestRuntime,
) -> Result<Vec<GeneratedJob>, String> {
    let jobs_root = layout.jobs_dir();
    fs::create_dir_all(&jobs_root).map_err(|err| err.to_string())?;
    fs::create_dir_all(layout.logs_dir()).map_err(|err| err.to_string())?;
    fs::create_dir_all(layout.tmp_dir()).map_err(|err| err.to_string())?;

    let mut generated = Vec::new();
    for job in &plan.jobs {
        generated.push(write_job(
            &jobs_root,
            job,
            prefilter,
            pruned_mzxml_paths,
            crosslinker,
            settings,
            fasta,
            runtime,
        )?);
    }
    Ok(generated)
}

fn write_job(
    jobs_root: &Path,
    job: &PlannedJob,
    prefilter: &PrefilterResult,
    pruned_mzxml_paths: &[PathBuf],
    crosslinker: &CrosslinkerProfile,
    settings: &Settings,
    fasta: &FastaDatabase,
    runtime: &XQuestRuntime,
) -> Result<GeneratedJob, String> {
    let job_dir = jobs_root.join(&job.job_id);
    fs::create_dir_all(&job_dir).map_err(|err| err.to_string())?;
    let results_dir = job_dir.join("results");
    fs::create_dir_all(&results_dir).map_err(|err| err.to_string())?;

    write_job_defs(
        &job_dir,
        &runtime.root,
        crosslinker,
        settings,
        &job.variant,
        &fasta.path,
    )?;

    let pruned = pruned_mzxml_for_job(job, pruned_mzxml_paths)?;
    let pruned_abs = pruned.canonicalize().map_err(|err| {
        format!(
            "cannot resolve pruned mzXML {}: {err}",
            pruned.display()
        )
    })?;
    let mzxml_link = job_dir.join("input.mzXML");
    symlink_or_copy(&pruned_abs, &mzxml_link)?;

    let rows = build_matchlist(job, prefilter, &PathBuf::from("input.mzXML"), crosslinker)?;
    let matchlist_path = job_dir.join("glycoquest_matched.txt");
    write_matchlist(&matchlist_path, &rows)?;

    let isotopepairs = isotopepairs_path(&matchlist_path);
    let specxml = specxml_filename("results");
    let compare_peaks = runtime.root.join("bin/compare_peaks3.pl");
    let xquest_exe = runtime.executable.canonicalize().unwrap_or_else(|_| runtime.executable.clone());
    let xquest_root = runtime.root.canonicalize().unwrap_or_else(|_| runtime.root.clone());
    let compare_peaks = compare_peaks.canonicalize().unwrap_or(compare_peaks);
    let def_path = job_dir.join("xquest.def");
    let run_script = job_dir.join("run.sh");
    let command = format!(
        "compare_peaks3.pl -match {matchlist} -def {def} -dir . -resultdir results -genxml input.mzXML -cpforce && \
         xquest.pl -def {def} -xquestdir {xquest_root} -list {isotopepairs} -resdir results -dir . -specxml {specxml} -nidx",
        matchlist = matchlist_path.display(),
        def = def_path.display(),
        xquest_root = xquest_root.display(),
        isotopepairs = isotopepairs.display(),
        specxml = specxml,
    );
    let perl5lib = format!(
        "{}/1209/lib/perl5:{}/1209/share/perl5:{}/modules",
        xquest_root.display(),
        xquest_root.display(),
        xquest_root.display()
    );
    let script = format!(
        "#!/bin/sh\nset -euo pipefail\nexport XQUEST_DIR=\"{xquest_root}\"\nexport PERL5LIB=\"{perl5lib}\"\n\
         \"{compare_peaks}\" -match glycoquest_matched.txt -def xquest.def -dir . -resultdir results -genxml input.mzXML -cpforce\n\
         if [ ! -s glycoquest_matched_isotopepairs.txt ]; then echo \"compare_peaks produced no spectra\" >&2; exit 1; fi\n\
         if [ ! -d results/db ]; then NIDX=-nidx; else NIDX=; fi\n\
         \"{xquest_exe}\" -def xquest.def -xquestdir \"$XQUEST_DIR\" -masstab \"$XQUEST_DIR/deffiles/mass_table.def\" -list glycoquest_matched_isotopepairs.txt -resdir results -dir . -specxml {specxml} $NIDX\n",
        xquest_root = xquest_root.display(),
        perl5lib = perl5lib,
        compare_peaks = compare_peaks.display(),
        xquest_exe = xquest_exe.display(),
        specxml = specxml,
    );
    fs::write(&run_script, script).map_err(|err| err.to_string())?;
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        fs::set_permissions(&run_script, fs::Permissions::from_mode(0o755))
            .map_err(|err| err.to_string())?;
    }

    Ok(GeneratedJob {
        job_id: job.job_id.clone(),
        directory: job_dir,
        run_script: run_script.clone(),
        matchlist_path,
        command,
    })
}

fn pruned_mzxml_for_job(job: &PlannedJob, paths: &[PathBuf]) -> Result<PathBuf, String> {
    let source = job
        .spectrum_keys
        .first()
        .ok_or_else(|| format!("job {} has no spectra", job.job_id))?
        .source_file
        .clone();

    paths
        .iter()
        .find(|path| {
            path.file_name()
                .zip(source.file_name())
                .map(|(a, b)| a == b)
                .unwrap_or(false)
        })
        .cloned()
        .ok_or_else(|| format!("no pruned mzXML for job {}", job.job_id))
}

#[cfg(unix)]
fn symlink_or_copy(source: &Path, dest: &Path) -> Result<(), String> {
    if dest.exists() {
        fs::remove_file(dest).ok();
    }
    std::os::unix::fs::symlink(source, dest).or_else(|_| {
        fs::copy(source, dest).map(|_| ()).map_err(|err| err.to_string())
    })
}

#[cfg(not(unix))]
fn symlink_or_copy(source: &Path, dest: &Path) -> Result<(), String> {
    fs::copy(source, dest).map_err(|err| err.to_string())?;
    Ok(())
}
