// Copyright (c) ETH Zurich, Mateusz Fido

//! Write xQuest job directories under `jobs/`.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use rayon::prelude::*;

use crate::cli::settings::Settings;
use crate::crosslinker::CrosslinkerProfile;
use crate::fasta::FastaDatabase;
use crate::jobs::{
    JobManifest, JobManifestEntry, JobPlan, PlannedJob, VarModPlan, build_varmod_plan,
};
use crate::prefilter::PrefilterResult;
use crate::progress::PhaseProgress;
use crate::xquest::XQuestRuntime;
use crate::xquest::defs::write_job_defs;
use crate::xquest::matchlist::{
    FilteredSpectrumIndex, build_matchlist, isotopepairs_path, specxml_filename, write_matchlist,
};

#[derive(Debug, Clone, PartialEq)]
pub struct GeneratedJob {
    pub job_id: String,
    pub directory: PathBuf,
    pub run_script: PathBuf,
    pub matchlist_path: PathBuf,
    pub command: String,
}

/// The generated job folders plus a manifest recording each job's glycan for
/// result annotation.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct GeneratedJobs {
    pub jobs: Vec<GeneratedJob>,
    pub manifest: JobManifest,
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
    progress: Option<&PhaseProgress>,
) -> Result<GeneratedJobs, String> {
    reset_job_workspace(layout)?;
    let jobs_root = layout.jobs_dir();
    fs::create_dir_all(&jobs_root).map_err(|err| err.to_string())?;
    fs::create_dir_all(layout.logs_dir()).map_err(|err| err.to_string())?;
    fs::create_dir_all(layout.tmp_dir()).map_err(|err| err.to_string())?;

    let spectrum_index = FilteredSpectrumIndex::new(prefilter);
    let completed = AtomicUsize::new(0);
    let generated_and_manifest: Result<Vec<_>, String> = plan
        .jobs
        .par_iter()
        .map(|job| {
            let varmod = build_varmod_plan(&job.variant, settings)?;
            let generated = write_job(
                &jobs_root,
                job,
                &varmod,
                prefilter,
                &spectrum_index,
                pruned_mzxml_paths,
                crosslinker,
                settings,
                fasta,
                runtime,
            )?;
            let manifest = JobManifestEntry {
                job_id: job.job_id.clone(),
                variant: job.variant.clone(),
                varmod_plan: varmod,
                source_file: job
                    .spectrum_keys
                    .first()
                    .map(|key| key.source_file.clone())
                    .unwrap_or_default(),
                spectrum_keys: job.spectrum_keys.clone(),
            };
            if let Some(progress) = progress {
                progress.inc(1);
                let count = completed.fetch_add(1, Ordering::Relaxed) + 1;
                progress.set_message(format!("{count} job folders written"));
            }
            Ok((generated, manifest))
        })
        .collect();
    let (generated, entries): (Vec<_>, Vec<_>) = generated_and_manifest?.into_iter().unzip();
    Ok(GeneratedJobs {
        jobs: generated,
        manifest: JobManifest { entries },
    })
}

fn reset_job_workspace(layout: &crate::output::ProjectLayout) -> Result<(), String> {
    let tmp = layout.tmp_dir();
    if tmp.is_dir() {
        fs::remove_dir_all(&tmp).map_err(|err| {
            format!(
                "cannot clear stale xQuest workspace {}: {err}",
                tmp.display()
            )
        })?;
    }
    fs::create_dir_all(&tmp)
        .map_err(|err| format!("cannot create xQuest workspace {}: {err}", tmp.display()))
}

fn write_job(
    jobs_root: &Path,
    job: &PlannedJob,
    varmod: &VarModPlan,
    prefilter: &PrefilterResult,
    spectrum_index: &FilteredSpectrumIndex<'_>,
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
        varmod,
        &fasta.path,
    )?;

    let pruned = pruned_mzxml_for_job(job, pruned_mzxml_paths)?;
    let pruned_abs = pruned
        .canonicalize()
        .map_err(|err| format!("cannot resolve pruned mzXML {}: {err}", pruned.display()))?;
    // Persist the real spectrum filename inside the job directory so it flows
    // through xQuest into the result XML instead of an opaque "input.mzXML".
    let mzxml_name = pruned_abs
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .ok_or_else(|| format!("invalid pruned mzXML path: {}", pruned_abs.display()))?;
    let mzxml_stem = pruned_abs
        .file_stem()
        .map(|stem| stem.to_string_lossy().into_owned())
        .unwrap_or_else(|| mzxml_name.clone());
    let mzxml_link = job_dir.join(&mzxml_name);
    symlink_or_copy(&pruned_abs, &mzxml_link)?;

    let rows = build_matchlist(job, prefilter, spectrum_index, &mzxml_stem, crosslinker)?;
    let matchlist_path = job_dir.join("glycoquest_matched.txt");
    write_matchlist(&matchlist_path, &rows)?;

    let isotopepairs = isotopepairs_path(&matchlist_path);
    let specxml = specxml_filename("results");
    let compare_peaks = runtime.root.join("bin/compare_peaks3.pl");
    let xquest_exe = runtime
        .executable
        .canonicalize()
        .unwrap_or_else(|_| runtime.executable.clone());
    let xquest_root = runtime
        .root
        .canonicalize()
        .unwrap_or_else(|_| runtime.root.clone());
    let compare_peaks = compare_peaks.canonicalize().unwrap_or(compare_peaks);
    let def_path = job_dir.join("xquest.def");
    let run_script = job_dir.join("run.sh");
    let command = format!(
        "compare_peaks3.pl -match {matchlist} -def {def} -dir . -resultdir results -genxml {mzxml} -cpforce && \
         xquest.pl -def {def} -xquestdir {xquest_root} -list {isotopepairs} -resdir results -dir . -specxml {specxml} -nidx",
        matchlist = matchlist_path.display(),
        def = def_path.display(),
        mzxml = mzxml_name,
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
         \"{compare_peaks}\" -match glycoquest_matched.txt -def xquest.def -dir . -resultdir results -genxml \"{mzxml}\" -cpforce\n\
         if [ ! -s glycoquest_matched_isotopepairs.txt ]; then echo \"compare_peaks produced no spectra\" >&2; exit 1; fi\n\
         if [ ! -d results/db ]; then NIDX=-nidx; else NIDX=; fi\n\
         \"{xquest_exe}\" -def xquest.def -xquestdir \"$XQUEST_DIR\" -masstab \"$XQUEST_DIR/deffiles/mass_table.def\" -list glycoquest_matched_isotopepairs.txt -resdir results -dir . -specxml {specxml} $NIDX\n",
        xquest_root = xquest_root.display(),
        perl5lib = perl5lib,
        compare_peaks = compare_peaks.display(),
        mzxml = mzxml_name,
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
        fs::copy(source, dest)
            .map(|_| ())
            .map_err(|err| err.to_string())
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reset_job_workspace_removes_stopped_run_jobs_and_logs() {
        let root = std::env::temp_dir().join(format!(
            "glycoquest_generate_reset_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let layout = crate::output::ProjectLayout::new(root.clone());
        let stale_job = layout.jobs_dir().join("stale-job");
        fs::create_dir_all(&stale_job).unwrap();
        fs::write(stale_job.join("run.sh"), "stale").unwrap();
        fs::create_dir_all(layout.logs_dir()).unwrap();
        fs::write(layout.logs_dir().join("stale.log"), "stale").unwrap();

        reset_job_workspace(&layout).unwrap();

        assert!(layout.tmp_dir().is_dir());
        assert!(!layout.jobs_dir().exists());
        assert!(!layout.logs_dir().exists());
        let _ = fs::remove_dir_all(root);
    }
}

#[cfg(not(unix))]
fn symlink_or_copy(source: &Path, dest: &Path) -> Result<(), String> {
    fs::copy(source, dest).map_err(|err| err.to_string())?;
    Ok(())
}
