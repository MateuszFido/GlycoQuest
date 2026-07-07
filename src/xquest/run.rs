//! Execute xQuest job scripts, continuing past individual failures.
//!
//! Jobs are independent (each runs in its own directory with a private result
//! folder and database copy), so they are executed concurrently on a rayon
//! thread pool sized by `parallelism`.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

use rayon::prelude::*;

/// Outcome of running one xQuest job folder's `run.sh`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JobRunRecord {
    pub job_name: String,
    pub job_dir: PathBuf,
    pub log_path: PathBuf,
    pub success: bool,
    pub error: Option<String>,
}

fn job_has_results(job_dir: &Path) -> bool {
    for candidate in [
        job_dir.join("results/xquest.xml"),
        job_dir.join("result.xml"),
    ] {
        if candidate.is_file() {
            if let Ok(meta) = std::fs::metadata(&candidate) {
                if meta.len() > 100 {
                    return true;
                }
            }
        }
    }
    false
}

/// Run every `run.sh` under `jobs_root`, logging to `logs_dir`.
///
/// Jobs run concurrently on a rayon thread pool. `parallelism` sets the number of
/// concurrent jobs; `0` uses one thread per available CPU core. Individual job
/// failures are recorded and do not stop other jobs. When `skip_completed` is
/// true, jobs that already have a non-empty result XML are skipped. The returned
/// records are sorted by job name for deterministic downstream output.
pub fn execute_jobs(
    jobs_root: &Path,
    logs_dir: &Path,
    skip_completed: bool,
    parallelism: usize,
) -> Result<Vec<JobRunRecord>, String> {
    std::fs::create_dir_all(logs_dir)
        .map_err(|err| format!("cannot create logs directory: {err}"))?;

    let entries = std::fs::read_dir(jobs_root)
        .map_err(|err| format!("cannot read jobs directory: {err}"))?;

    let mut job_scripts: Vec<_> = entries
        .flatten()
        .filter_map(|entry| {
            let run_script = entry.path().join("run.sh");
            run_script.is_file().then_some((entry.path(), run_script))
        })
        .collect();
    job_scripts.sort_by(|(a, _), (b, _)| a.cmp(b));

    let total_jobs = job_scripts.len();
    if total_jobs == 0 {
        return Ok(Vec::new());
    }

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(parallelism)
        .build()
        .map_err(|err| format!("cannot build job thread pool: {err}"))?;

    eprintln!(
        "run: executing {total_jobs} xQuest jobs across {} worker thread(s)",
        pool.current_num_threads()
    );

    let completed = AtomicUsize::new(0);

    let mut records: Vec<JobRunRecord> = pool.install(|| {
        job_scripts
            .par_iter()
            .map(|(job_dir, run_script)| {
                run_single_job(
                    job_dir,
                    run_script,
                    logs_dir,
                    skip_completed,
                    total_jobs,
                    &completed,
                )
            })
            .collect()
    });

    records.sort_by(|a, b| a.job_name.cmp(&b.job_name));
    Ok(records)
}

fn run_single_job(
    job_dir: &Path,
    run_script: &Path,
    logs_dir: &Path,
    skip_completed: bool,
    total_jobs: usize,
    completed: &AtomicUsize,
) -> JobRunRecord {
    let job_name = job_dir
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_default();
    let log_path = logs_dir.join(format!("{job_name}.log"));

    if skip_completed && job_has_results(job_dir) {
        let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
        eprintln!("run: job {done}/{total_jobs} skipped {job_name} (result present)");
        return JobRunRecord {
            job_name,
            job_dir: job_dir.to_path_buf(),
            log_path,
            success: true,
            error: None,
        };
    }

    let record = match Command::new("sh")
        .arg("run.sh")
        .current_dir(job_dir)
        .output()
    {
        Ok(output) => {
            let mut log = output.stdout;
            log.extend_from_slice(&output.stderr);
            let _ = std::fs::write(&log_path, log);

            if output.status.success() {
                JobRunRecord {
                    job_name: job_name.clone(),
                    job_dir: job_dir.to_path_buf(),
                    log_path,
                    success: true,
                    error: None,
                }
            } else {
                let reason = format!("exit status {output}", output = output.status);
                eprintln!(
                    "warning: xQuest job failed ({}) — see {}",
                    job_dir.display(),
                    log_path.display()
                );
                JobRunRecord {
                    job_name: job_name.clone(),
                    job_dir: job_dir.to_path_buf(),
                    log_path,
                    success: false,
                    error: Some(reason),
                }
            }
        }
        Err(err) => {
            eprintln!("warning: failed to launch {}: {err}", run_script.display());
            JobRunRecord {
                job_name: job_name.clone(),
                job_dir: job_dir.to_path_buf(),
                log_path,
                success: false,
                error: Some(err.to_string()),
            }
        }
    };

    let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
    let status = if record.success { "ok" } else { "failed" };
    eprintln!("run: job {done}/{total_jobs} finished {job_name} ({status})");

    record
}

/// Print a summary of job outcomes to stderr.
pub fn log_job_summary(records: &[JobRunRecord]) {
    let total = records.len();
    let ok = records.iter().filter(|record| record.success).count();
    let failed = total.saturating_sub(ok);

    if total == 0 {
        eprintln!("run: no xQuest jobs to execute");
        return;
    }

    eprintln!("run: xQuest jobs finished: {ok} succeeded, {failed} failed (of {total})");

    if failed == total {
        eprintln!("run: all xQuest jobs failed; consolidated hits will be empty");
    } else if failed > 0 {
        eprintln!("run: see results/failed_jobs.tsv for failed job details");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::sync::atomic::{AtomicU64, Ordering};

    static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

    fn temp_dir(name: &str) -> PathBuf {
        let id = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir().join(format!(
            "glycoquest_xquest_run_test_{}_{}_{}",
            std::process::id(),
            name,
            id
        ))
    }

    fn write_run_sh(job_dir: &Path, body: &str) {
        fs::create_dir_all(job_dir).unwrap();
        let script = job_dir.join("run.sh");
        fs::write(&script, body).unwrap();
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            fs::set_permissions(&script, fs::Permissions::from_mode(0o755)).unwrap();
        }
    }

    #[test]
    fn continues_after_job_failure() {
        let root = temp_dir("continue");
        let jobs_root = root.join("jobs");
        let logs_dir = root.join("logs");
        write_run_sh(&jobs_root.join("job_ok"), "#!/bin/sh\nexit 0\n");
        write_run_sh(&jobs_root.join("job_fail"), "#!/bin/sh\nexit 1\n");
        write_run_sh(&jobs_root.join("job_ok2"), "#!/bin/sh\nexit 0\n");

        let records = execute_jobs(&jobs_root, &logs_dir, false, 2).unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records.iter().filter(|r| r.success).count(), 2);
        let failed = records
            .iter()
            .find(|record| record.job_name == "job_fail")
            .expect("job_fail record");
        assert!(!failed.success);
        assert!(failed.error.as_ref().is_some_and(|e| e.contains("exit status")));

        let _ = fs::remove_dir_all(root);
    }

    #[test]
    fn empty_jobs_dir_returns_no_records() {
        let root = temp_dir("empty");
        let jobs_root = root.join("jobs");
        let logs_dir = root.join("logs");
        fs::create_dir_all(&jobs_root).unwrap();

        let records = execute_jobs(&jobs_root, &logs_dir, false, 1).unwrap();
        assert!(records.is_empty());

        let _ = fs::remove_dir_all(root);
    }

    #[test]
    fn skips_jobs_with_existing_results_when_requested() {
        let root = temp_dir("skip");
        let jobs_root = root.join("jobs");
        let logs_dir = root.join("logs");
        let done_dir = jobs_root.join("job_done");
        fs::create_dir_all(done_dir.join("results")).unwrap();
        fs::write(
            done_dir.join("results/xquest.xml"),
            format!("<xquest_results>{}</xquest_results>", "x".repeat(120)),
        )
        .unwrap();
        write_run_sh(&done_dir, "#!/bin/sh\nexit 0\n");
        write_run_sh(&jobs_root.join("job_run"), "#!/bin/sh\nexit 0\n");

        let records = execute_jobs(&jobs_root, &logs_dir, true, 2).unwrap();
        assert_eq!(records.len(), 2);
        assert!(records.iter().all(|record| record.success));
        assert_eq!(records[0].job_name, "job_done");

        let _ = fs::remove_dir_all(root);
    }
}
