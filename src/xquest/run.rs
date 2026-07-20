// Copyright (c) ETH Zurich, Mateusz Fido

//! Execute xQuest job scripts, continuing past individual failures.
//!
//! Jobs are independent (each runs in its own directory with a private result
//! folder and database copy), so they are executed concurrently on a rayon
//! thread pool sized by `parallelism`.

use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Duration;

use rayon::prelude::*;

use crate::progress::PhaseProgress;

const PROGRESS_POLL_INTERVAL: Duration = Duration::from_millis(250);

/// Outcome of running one xQuest job folder's `run.sh`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JobRunRecord {
    pub job_name: String,
    pub job_dir: PathBuf,
    pub log_path: PathBuf,
    pub success: bool,
    pub error: Option<String>,
}

/// Run every `run.sh` under `jobs_root`, logging to `logs_dir`.
///
/// Jobs run concurrently on a rayon thread pool. `parallelism` sets the number of
/// concurrent jobs; `0` uses one thread per available CPU core. Individual job
/// failures are recorded and do not stop other jobs. The returned records are
/// sorted by job name for deterministic downstream output.
#[cfg(test)]
fn execute_jobs(
    jobs_root: &Path,
    logs_dir: &Path,
    parallelism: usize,
) -> Result<Vec<JobRunRecord>, String> {
    execute_jobs_with_progress(jobs_root, logs_dir, parallelism, &HashMap::new(), None)
}

pub(crate) fn execute_jobs_with_progress(
    jobs_root: &Path,
    logs_dir: &Path,
    parallelism: usize,
    job_work: &HashMap<String, u64>,
    progress: Option<&PhaseProgress>,
) -> Result<Vec<JobRunRecord>, String> {
    std::fs::create_dir_all(logs_dir)
        .map_err(|err| format!("cannot create logs directory: {err}"))?;

    let entries =
        std::fs::read_dir(jobs_root).map_err(|err| format!("cannot read jobs directory: {err}"))?;

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

    let execution_message = format!(
        "run: executing {total_jobs} xQuest jobs across {} worker thread(s)",
        pool.current_num_threads()
    );
    if let Some(progress) = progress.filter(|progress| progress.enabled()) {
        progress.println(&execution_message);
    } else {
        eprintln!("{execution_message}");
    }

    let completed = AtomicUsize::new(0);
    let tracker = progress
        .filter(|progress| progress.enabled())
        .map(|progress| {
            Arc::new(JobProgressTracker::new(
                progress.clone(),
                total_jobs,
                if job_work.is_empty() {
                    total_jobs as u64
                } else {
                    job_work.values().copied().sum()
                },
            ))
        });
    let stop_monitor = Arc::new(AtomicBool::new(false));
    let monitor = tracker
        .as_ref()
        .map(|tracker| spawn_progress_monitor(Arc::clone(tracker), Arc::clone(&stop_monitor)));

    let mut records: Vec<JobRunRecord> = pool.install(|| {
        job_scripts
            .par_iter()
            .map(|(job_dir, run_script)| {
                run_single_job(
                    job_dir,
                    run_script,
                    logs_dir,
                    total_jobs,
                    &completed,
                    job_weight(job_dir, job_work),
                    tracker.as_deref(),
                )
            })
            .collect()
    });

    stop_monitor.store(true, Ordering::Release);
    if let Some(monitor) = monitor {
        let _ = monitor.join();
    }
    if let Some(tracker) = tracker {
        tracker.refresh_and_render();
    }

    records.sort_by(|a, b| a.job_name.cmp(&b.job_name));
    Ok(records)
}

fn run_single_job(
    job_dir: &Path,
    run_script: &Path,
    logs_dir: &Path,
    total_jobs: usize,
    completed: &AtomicUsize,
    work: u64,
    tracker: Option<&JobProgressTracker>,
) -> JobRunRecord {
    let job_name = job_dir
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_default();
    let log_path = logs_dir.join(format!("{job_name}.log"));
    let progress_path = job_dir.join("results/results.progress");
    let _ = std::fs::remove_file(&progress_path);
    if let Some(tracker) = tracker {
        tracker.start_job(&job_name, progress_path, work);
    }

    let record = match run_command_to_log(job_dir, &log_path) {
        Ok(status) => {
            if status.success() {
                JobRunRecord {
                    job_name: job_name.clone(),
                    job_dir: job_dir.to_path_buf(),
                    log_path,
                    success: true,
                    error: None,
                }
            } else {
                let reason = format!("exit status {status}");
                log_progress_aware(
                    tracker,
                    format!(
                        "warning: xQuest job failed ({}) — see {}",
                        job_dir.display(),
                        log_path.display()
                    ),
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
            log_progress_aware(
                tracker,
                format!("warning: failed to launch {}: {err}", run_script.display()),
            );
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
    if let Some(tracker) = tracker {
        tracker.finish_job(&job_name, record.success);
    } else {
        eprintln!("run: job {done}/{total_jobs} finished {job_name} ({status})");
    }

    record
}

fn run_command_to_log(job_dir: &Path, log_path: &Path) -> Result<std::process::ExitStatus, String> {
    let stdout = File::create(log_path)
        .map_err(|err| format!("cannot create log {}: {err}", log_path.display()))?;
    let stderr = stdout
        .try_clone()
        .map_err(|err| format!("cannot clone log handle {}: {err}", log_path.display()))?;

    Command::new("sh")
        .arg("run.sh")
        .current_dir(job_dir)
        .stdout(Stdio::from(stdout))
        .stderr(Stdio::from(stderr))
        .status()
        .map_err(|err| err.to_string())
}

fn job_weight(job_dir: &Path, job_work: &HashMap<String, u64>) -> u64 {
    let name = job_dir
        .file_name()
        .map(|name| name.to_string_lossy())
        .unwrap_or_default();
    job_work.get(name.as_ref()).copied().unwrap_or(1).max(1)
}

fn log_progress_aware(tracker: Option<&JobProgressTracker>, message: String) {
    if let Some(tracker) = tracker {
        tracker.phase.println(message);
    } else {
        eprintln!("{message}");
    }
}

#[derive(Debug, Clone)]
struct ActiveJobProgress {
    progress_path: PathBuf,
    work: u64,
    partial_work: u64,
}

#[derive(Debug, Default)]
struct JobProgressState {
    active: HashMap<String, ActiveJobProgress>,
    completed_jobs: usize,
    failed_jobs: usize,
    completed_work: u64,
}

#[derive(Debug)]
struct JobProgressTracker {
    phase: PhaseProgress,
    total_jobs: usize,
    total_work: u64,
    state: Mutex<JobProgressState>,
}

impl JobProgressTracker {
    fn new(phase: PhaseProgress, total_jobs: usize, total_work: u64) -> Self {
        let tracker = Self {
            phase,
            total_jobs,
            total_work: total_work.max(1),
            state: Mutex::new(JobProgressState::default()),
        };
        tracker.render();
        tracker
    }

    fn start_job(&self, name: &str, progress_path: PathBuf, work: u64) {
        self.state.lock().expect("job progress lock").active.insert(
            name.to_string(),
            ActiveJobProgress {
                progress_path,
                work,
                partial_work: 0,
            },
        );
        self.render();
    }

    fn finish_job(&self, name: &str, success: bool) {
        let mut state = self.state.lock().expect("job progress lock");
        if let Some(active) = state.active.remove(name) {
            let completed = if success {
                active.work
            } else {
                active.partial_work
            };
            state.completed_work = state.completed_work.saturating_add(completed);
        }
        state.completed_jobs += 1;
        if !success {
            state.failed_jobs += 1;
        }
        drop(state);
        self.render();
    }

    fn refresh_and_render(&self) {
        let active: Vec<_> = self
            .state
            .lock()
            .expect("job progress lock")
            .active
            .iter()
            .map(|(name, active)| (name.clone(), active.progress_path.clone(), active.work))
            .collect();

        let updates: Vec<_> = active
            .into_iter()
            .filter_map(|(name, path, work)| {
                let text = std::fs::read_to_string(path).ok()?;
                let parsed = parse_xquest_progress(&text)?;
                let partial_work = match parsed {
                    XQuestProgress::Searching { current, total } => {
                        work.saturating_mul(current.saturating_sub(1)) / total.max(1)
                    }
                    XQuestProgress::Finished => work,
                };
                Some((name, partial_work))
            })
            .collect();

        if !updates.is_empty() {
            let mut state = self.state.lock().expect("job progress lock");
            for (name, partial_work) in updates {
                if let Some(active) = state.active.get_mut(&name) {
                    active.partial_work = active.partial_work.max(partial_work.min(active.work));
                }
            }
        }
        self.render();
    }

    fn render(&self) {
        let state = self.state.lock().expect("job progress lock");
        let position = state
            .active
            .values()
            .fold(state.completed_work, |total, active| {
                total.saturating_add(active.partial_work)
            })
            .min(self.total_work);
        self.phase.set_position(position);
        self.phase.set_message(format!(
            "spectrum assignments · {}/{} jobs · {} active · {} failed",
            state.completed_jobs,
            self.total_jobs,
            state.active.len(),
            state.failed_jobs
        ));
    }
}

fn spawn_progress_monitor(
    tracker: Arc<JobProgressTracker>,
    stop: Arc<AtomicBool>,
) -> thread::JoinHandle<()> {
    thread::spawn(move || {
        while !stop.load(Ordering::Acquire) {
            tracker.refresh_and_render();
            thread::sleep(PROGRESS_POLL_INTERVAL);
        }
        tracker.refresh_and_render();
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum XQuestProgress {
    Searching { current: u64, total: u64 },
    Finished,
}

fn parse_xquest_progress(text: &str) -> Option<XQuestProgress> {
    if text.contains("Search finished:") {
        return Some(XQuestProgress::Finished);
    }

    let tokens: Vec<_> = text.split_whitespace().collect();
    for window in tokens.windows(5) {
        if window[0] == "Searching" && window[1] == "spectrum" && window[3] == "of" {
            let current = numeric_token(window[2])?;
            let total = numeric_token(window[4])?;
            if total > 0 {
                return Some(XQuestProgress::Searching { current, total });
            }
        }
    }
    None
}

fn numeric_token(token: &str) -> Option<u64> {
    let digits: String = token.chars().filter(char::is_ascii_digit).collect();
    (!digits.is_empty()).then(|| digits.parse().ok()).flatten()
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

        let records = execute_jobs(&jobs_root, &logs_dir, 2).unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records.iter().filter(|r| r.success).count(), 2);
        let failed = records
            .iter()
            .find(|record| record.job_name == "job_fail")
            .expect("job_fail record");
        assert!(!failed.success);
        assert!(
            failed
                .error
                .as_ref()
                .is_some_and(|e| e.contains("exit status"))
        );

        let _ = fs::remove_dir_all(root);
    }

    #[test]
    fn empty_jobs_dir_returns_no_records() {
        let root = temp_dir("empty");
        let jobs_root = root.join("jobs");
        let logs_dir = root.join("logs");
        fs::create_dir_all(&jobs_root).unwrap();

        let records = execute_jobs(&jobs_root, &logs_dir, 1).unwrap();
        assert!(records.is_empty());

        let _ = fs::remove_dir_all(root);
    }

    #[test]
    fn parses_native_xquest_progress_messages() {
        assert_eq!(
            parse_xquest_progress("\n ### Searching spectrum 12 of 53... ###"),
            Some(XQuestProgress::Searching {
                current: 12,
                total: 53
            })
        );
        assert_eq!(
            parse_xquest_progress("Search finished: 53 spectra were searched."),
            Some(XQuestProgress::Finished)
        );
        assert_eq!(parse_xquest_progress(""), None);
        assert_eq!(parse_xquest_progress("Searching spectrum 1 of 0"), None);
    }

    #[test]
    fn failed_job_keeps_only_observed_partial_progress() {
        let root = temp_dir("failed_progress");
        let progress_path = root.join("results.progress");
        fs::create_dir_all(&root).unwrap();
        fs::write(&progress_path, "Searching spectrum 6 of 20").unwrap();

        let phase = crate::progress::ProgressReporter::new(crate::ProgressMode::Never).determinate(
            3,
            4,
            "xQuest searches",
            100,
        );
        let tracker = JobProgressTracker::new(phase.clone(), 1, 100);
        tracker.start_job("failed", progress_path, 100);
        tracker.refresh_and_render();
        tracker.finish_job("failed", false);

        let (position, message) = phase.test_snapshot();
        assert_eq!(position, 25);
        assert!(message.contains("1 failed"));

        let _ = fs::remove_dir_all(root);
    }

    #[test]
    fn child_output_is_written_to_log() {
        let root = temp_dir("stream_log");
        let jobs_root = root.join("jobs");
        let logs_dir = root.join("logs");
        write_run_sh(
            &jobs_root.join("job_output"),
            "#!/bin/sh\necho standard-output\necho standard-error >&2\n",
        );

        let records = execute_jobs(&jobs_root, &logs_dir, 1).unwrap();
        assert!(records[0].success);
        let log = fs::read_to_string(&records[0].log_path).unwrap();
        assert!(log.contains("standard-output"));
        assert!(log.contains("standard-error"));

        let _ = fs::remove_dir_all(root);
    }
}
