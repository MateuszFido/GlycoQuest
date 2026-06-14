//! GlycoQuest library: CLI parameter types, settings, and the entry-point runner.

mod cli;
mod crosslinker;
mod fasta;
mod glyco;
mod jobs;
mod mzxml;
mod output;
mod prefilter;
mod results;
mod xquest;

pub use cli::{parse_cli, CliParams};
pub use crosslinker::CrosslinkerProfile;
pub use fasta::{FastaDatabase, validate_fasta};
pub use glyco::{
    glycan_data_dir, load_glycan_database, resolve_database, supported_glycan_databases,
    DiagnosticIon, GlycanEntry, GlycanLibrary,
};
pub use cli::input::resolve_input;
pub use cli::settings::{default_settings_path, Settings};
pub use mzxml::{Ms2Scan, parse_scans, write_prefiltered_mzxml};
pub use prefilter::{PrefilterResult, PrunedGlycan, run_prefilter, write_outputs};
pub use xquest::{resolve_runtime, XQuestRuntime};

/// Whether to validate only or execute a full search.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExecutionMode {
    DryRun,
    Run,
}

/// Exit codes
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExitCode {
    Success = 0,
    Validation = 1,
    NoSpectra = 2,
    XquestExecution = 3,
    ResultExtraction = 4,
}

impl From<ExitCode> for i32 {
    fn from(code: ExitCode) -> Self {
        code as i32
    }
}

/// CLI arguments merged with values from `settings.ini` (CLI has precedence)
#[derive(Debug, Clone, PartialEq)]
pub struct RunConfig {
    pub cli: CliParams,
    pub settings: Settings,
    pub crosslinker: CrosslinkerProfile,
}

impl RunConfig {
    pub fn load(cli: CliParams) -> Result<Self, String> {
        let config_path = cli
            .config
            .clone()
            .unwrap_or_else(default_settings_path);
        let mut settings = if config_path.is_file() {
            Settings::load_from_file(&config_path)?
        } else {
            eprintln!(
                "warning: settings file {} not found; using built-in defaults",
                config_path.display()
            );
            Settings::defaults()
        };

        if let Some(ppm) = cli.ppm_tolerance {
            settings.diagnostic_tolerance_ppm = ppm;
        }

        let resolved_out = output::resolve_project_out_dir(&cli.out, &cli.input, &cli.database)?;
        let mut cli = cli;
        cli.out = resolved_out;

        let crosslinker = CrosslinkerProfile::resolve(
            &settings,
            cli.crosslinker.as_deref(),
        )?;

        if crosslinker.label == crate::crosslinker::CrosslinkerLabel::None
            && settings.crosslinker_shift_da.abs() > f64::EPSILON
        {
            eprintln!(
                "warning: crosslinker shift_da={} is ignored for isotope prefilter when label=none",
                settings.crosslinker_shift_da
            );
        }

        Ok(Self {
            cli,
            settings,
            crosslinker,
        })
    }

    pub fn execution_mode(&self) -> ExecutionMode {
        if self.cli.dry_run {
            ExecutionMode::DryRun
        } else {
            ExecutionMode::Run
        }
    }
}

/// Run GlycoQuest from parsed CLI parameters (loads `settings.ini` automatically).
pub fn run(cli: &CliParams) -> i32 {
    match RunConfig::load(cli.clone()) {
        Ok(config) => run_config(&config),
        Err(message) => {
            eprintln!("error: {message}");
            ExitCode::Validation.into()
        }
    }
}

/// Validated inputs ready for prefiltering or xQuest job generation.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidatedInputs {
    pub files: Vec<std::path::PathBuf>,
    pub fasta: FastaDatabase,
    pub glycan_library: GlycanLibrary,
    pub xquest: XQuestRuntime,
    pub layout: output::ProjectLayout,
}

/// Outcome of a single readiness check.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadinessCheck {
    pub label: &'static str,
    pub ok: bool,
    pub error: Option<String>,
}

/// Best-effort assessment used for config summary and readiness reporting.
#[derive(Debug, Clone, PartialEq)]
pub struct ConfigAssessment {
    pub ms_files: Option<Vec<std::path::PathBuf>>,
    pub fasta: Option<FastaDatabase>,
    pub glycan_library: Option<GlycanLibrary>,
    pub xquest: Option<XQuestRuntime>,
    pub checks: Vec<ReadinessCheck>,
}

impl ConfigAssessment {
    pub fn all_passed(&self) -> bool {
        self.checks.iter().all(|check| check.ok)
    }

    pub fn first_error(&self) -> Option<&str> {
        self.checks
            .iter()
            .find(|check| !check.ok)
            .and_then(|check| check.error.as_deref())
    }
}

/// Run GlycoQuest from a fully merged [`RunConfig`].
pub fn run_config(config: &RunConfig) -> i32 {
    let assessment = assess_config(config);
    print_config_summary(config, &assessment);
    print_readiness_report(&assessment);

    let validated = match assessment.try_into_validated(config) {
        Ok(validated) => validated,
        Err(message) => {
            eprintln!("error: {message}");
            return ExitCode::Validation.into();
        }
    };

    match config.execution_mode() {
        ExecutionMode::DryRun => dry_run(config, &validated),
        ExecutionMode::Run => run_search(config, &validated),
    }
}

fn dry_run(config: &RunConfig, validated: &ValidatedInputs) -> i32 {
    match execute_pipeline(config, validated, false) {
        Ok(prefilter) => {
            if prefilter.filtered.is_empty() {
                ExitCode::NoSpectra.into()
            } else {
                ExitCode::Success.into()
            }
        }
        Err(code) => code.into(),
    }
}

fn run_search(config: &RunConfig, validated: &ValidatedInputs) -> i32 {
    let prefilter = match execute_pipeline(config, validated, true) {
        Ok(prefilter) => prefilter,
        Err(code) => return code.into(),
    };

    if prefilter.filtered.is_empty() {
        return ExitCode::NoSpectra.into();
    }

    let layout = &validated.layout;
    let jobs_root = layout.jobs_dir();
    let logs_dir = layout.logs_dir();
    if let Err(message) = std::fs::create_dir_all(&logs_dir) {
        eprintln!("error: cannot create logs directory: {message}");
        return ExitCode::XquestExecution.into();
    }

    let entries = match std::fs::read_dir(&jobs_root) {
        Ok(entries) => entries,
        Err(err) => {
            eprintln!("error: cannot read jobs directory: {err}");
            return ExitCode::XquestExecution.into();
        }
    };

    let mut job_scripts: Vec<_> = entries
        .flatten()
        .filter_map(|entry| {
            let run_script = entry.path().join("run.sh");
            run_script.is_file().then_some((entry, run_script))
        })
        .collect();
    job_scripts.sort_by_key(|(entry, _)| entry.file_name());

    let total_jobs = job_scripts.len();
    for (index, (entry, run_script)) in job_scripts.into_iter().enumerate() {
        eprintln!(
            "run: job {}/{} executing {}",
            index + 1,
            total_jobs,
            run_script.display()
        );
        let log_path = logs_dir.join(format!(
            "{}.log",
            entry.file_name().to_string_lossy()
        ));
        let output = std::process::Command::new("sh")
            .arg("run.sh")
            .current_dir(entry.path())
            .output();
        match output {
            Ok(output) => {
                let mut log = output.stdout;
                log.extend_from_slice(&output.stderr);
                let _ = std::fs::write(&log_path, log);
                if !output.status.success() {
                    eprintln!(
                        "error: xQuest job failed ({}) — see {}",
                        entry.path().display(),
                        log_path.display()
                    );
                    return ExitCode::XquestExecution.into();
                }
            }
            Err(err) => {
                eprintln!("error: failed to launch {}: {err}", run_script.display());
                return ExitCode::XquestExecution.into();
            }
        }
    }

    match results::consolidate_results(
        layout,
        &config.settings,
        &config.crosslinker,
        &prefilter,
    ) {
        Ok(()) => {
            output::cleanup_temp_artifacts(layout);
            ExitCode::Success.into()
        }
        Err(message) => {
            eprintln!("error: {message}");
            ExitCode::ResultExtraction.into()
        }
    }
}

fn execute_pipeline(
    config: &RunConfig,
    validated: &ValidatedInputs,
    _run_jobs: bool,
) -> Result<PrefilterResult, ExitCode> {
    let layout = &validated.layout;
    let prefilter = match run_prefilter(
        &validated.files,
        &validated.glycan_library,
        &config.settings,
        &config.crosslinker,
    ) {
        Ok(result) => result,
        Err(message) => {
            eprintln!("error: {message}");
            return Err(ExitCode::Validation);
        }
    };

    if let Err(message) = write_outputs(layout.root.as_path(), &prefilter) {
        eprintln!("error: {message}");
        return Err(ExitCode::Validation);
    }

    print_prefilter_summary(&prefilter);

    if prefilter.filtered.is_empty() {
        eprintln!("dry-run: no spectra retained after prefilter");
        return Ok(prefilter);
    }

    let pruned_paths = match write_prefiltered_mzxml(layout.spectra_dir().as_path(), &prefilter.filtered) {
        Ok(paths) => paths,
        Err(message) => {
            eprintln!("error: {message}");
            return Err(ExitCode::Validation);
        }
    };

    let job_plan = match jobs::JobPlan::build(
        &prefilter,
        &validated.glycan_library,
        &config.settings,
        &config.crosslinker,
    ) {
        Ok(plan) => plan,
        Err(message) => {
            eprintln!("error: {message}");
            return Err(ExitCode::Validation);
        }
    };

    let generated = match xquest::generate_jobs(
        layout,
        &job_plan,
        &prefilter,
        &pruned_paths,
        &config.crosslinker,
        &config.settings,
        &validated.fasta,
        &validated.xquest,
    ) {
        Ok(jobs) => jobs,
        Err(message) => {
            eprintln!("error: {message}");
            return Err(ExitCode::Validation);
        }
    };

    let plan_doc = output::build_run_plan_document(
        &config.crosslinker,
        &prefilter,
        &job_plan,
        pruned_paths,
        generated,
    );
    if let Err(message) = output::write_plan_json(layout.root.as_path(), &plan_doc) {
        eprintln!("error: {message}");
        return Err(ExitCode::Validation);
    }

    eprintln!(
        "plan: {} jobs, {} estimated comparisons, isotope_prefilter={}",
        plan_doc.job_count,
        plan_doc.total_comparisons,
        plan_doc.isotope_prefilter_enabled
    );

    if config.execution_mode() == ExecutionMode::DryRun {
        eprintln!("dry-run: job folders and plan.json written (xQuest not executed)");
    }

    Ok(prefilter)
}

fn assess_config(config: &RunConfig) -> ConfigAssessment {
    let mut checks = Vec::new();

    let ms_files = match resolve_input(&config.cli.input) {
        Ok(files) => {
            checks.push(ReadinessCheck {
                label: "MS input",
                ok: true,
                error: None,
            });
            Some(files)
        }
        Err(err) => {
            checks.push(ReadinessCheck {
                label: "MS input",
                ok: false,
                error: Some(err),
            });
            None
        }
    };

    let fasta = match validate_fasta(&config.cli.database) {
        Ok(db) => {
            checks.push(ReadinessCheck {
                label: "FASTA database",
                ok: true,
                error: None,
            });
            Some(db)
        }
        Err(err) => {
            checks.push(ReadinessCheck {
                label: "FASTA database",
                ok: false,
                error: Some(err),
            });
            None
        }
    };

    let glycan_library = match load_glycan_database(&config.cli.glycans) {
        Ok(library) => {
            checks.push(ReadinessCheck {
                label: "Glycan library",
                ok: true,
                error: None,
            });
            Some(library)
        }
        Err(err) => {
            checks.push(ReadinessCheck {
                label: "Glycan library",
                ok: false,
                error: Some(err),
            });
            None
        }
    };

    let xquest = match resolve_runtime(&config.cli.xquest_root, &config.settings) {
        Ok(runtime) => {
            checks.push(ReadinessCheck {
                label: "xQuest runtime",
                ok: true,
                error: None,
            });
            Some(runtime)
        }
        Err(err) => {
            checks.push(ReadinessCheck {
                label: "xQuest runtime",
                ok: false,
                error: Some(err),
            });
            None
        }
    };

    match output::ensure_output_dir(&config.cli.out) {
        Ok(()) => {
            checks.push(ReadinessCheck {
                label: "Output directory",
                ok: true,
                error: None,
            });
        }
        Err(err) => {
            checks.push(ReadinessCheck {
                label: "Output directory",
                ok: false,
                error: Some(err),
            });
        }
    }

    ConfigAssessment {
        ms_files,
        fasta,
        glycan_library,
        xquest,
        checks,
    }
}

impl ConfigAssessment {
    fn try_into_validated(self, config: &RunConfig) -> Result<ValidatedInputs, String> {
        if !self.all_passed() {
            return Err(
                self.first_error()
                    .unwrap_or("configuration validation failed")
                    .to_string(),
            );
        }

        Ok(ValidatedInputs {
            files: self.ms_files.expect("checked"),
            fasta: self.fasta.expect("checked"),
            glycan_library: self.glycan_library.expect("checked"),
            xquest: self.xquest.expect("checked"),
            layout: output::ProjectLayout::new(config.cli.out.clone()),
        })
    }
}

fn print_config_summary(config: &RunConfig, assessment: &ConfigAssessment) {
    for line in format_config_summary(config, assessment) {
        eprintln!("{line}");
    }
}

fn print_readiness_report(assessment: &ConfigAssessment) {
    for line in format_readiness_report(assessment) {
        eprintln!("{line}");
    }
}

fn format_config_summary(config: &RunConfig, assessment: &ConfigAssessment) -> Vec<String> {
    let cli = &config.cli;
    let settings = &config.settings;
    let crosslinker = &config.crosslinker;
    let mut lines = vec![
        format!("input: {}", cli.input.display()),
        format!(
            "MS files: {} ({})",
            ms_file_count(assessment),
            ms_file_label(assessment)
        ),
    ];

    match &assessment.ms_files {
        Some(files) => {
            for path in files {
                lines.push(format!("  {}", path.display()));
            }
        }
        None => lines.push("  None".to_string()),
    }

    lines.push(format!("database: {}", cli.database.display()));
    lines.push(format!(
        "FASTA entries: {} ({})",
        fasta_entry_count(assessment),
        optional_label(assessment.fasta.as_ref().map(|_| "loaded"))
    ));
    lines.push(format!("glycans: {}", cli.glycans));
    lines.push(format!(
        "glycan entries: {} ({})",
        glycan_entry_count(assessment),
        optional_label(assessment.glycan_library.as_ref().map(|_| "loaded"))
    ));

    if let Some(library) = &assessment.glycan_library {
        for entry in library.entries.iter().take(3) {
            lines.push(format!(
                "  {}  mass={}  diagnostics={}",
                entry.composition,
                entry.monoisotopic_mass,
                entry.diagnostic_ions.len()
            ));
        }
        if library.entries.len() > 3 {
            lines.push("  …".to_string());
        }
    }

    lines.push(format!("xquest_root: {}", cli.xquest_root.display()));
    lines.push(format!(
        "xquest_bin: {}",
        assessment
            .xquest
            .as_ref()
            .map(|runtime| runtime.executable.display().to_string())
            .unwrap_or_else(|| "None".to_string())
    ));
    if let Some(path) = &settings.xquest_bin {
        lines.push(format!("xquest_bin (settings): {}", path.display()));
    }
    lines.push(format!("out: {}", cli.out.display()));
    lines.push(format!("spectra: {}/spectra", cli.out.display()));
    lines.push(format!("tmp: {}/tmp", cli.out.display()));
    lines.push(format!("crosslinker: {}", crosslinker.name));
    lines.push(format!("crosslinker_label: {}", crosslinker.label.as_str()));
    lines.push(format!(
        "isotope_prefilter: {}",
        if crosslinker.requires_isotope_pair_prefilter() {
            "enabled"
        } else {
            "disabled"
        }
    ));
    lines.push(format!("crosslinker_shift_da: {}", crosslinker.shift_da));
    lines.push(format!(
        "diagnostic_tolerance_ppm: {}",
        settings.diagnostic_tolerance_ppm
    ));

    lines
}

fn format_readiness_report(assessment: &ConfigAssessment) -> Vec<String> {
    let mut lines = vec!["readiness:".to_string()];
    for check in &assessment.checks {
        lines.push(format_readiness_line(check));
    }
    lines.push(format_overall_readiness(assessment.all_passed()));
    lines
}

fn format_readiness_line(check: &ReadinessCheck) -> String {
    if check.ok {
        format!("{}  {}  {}", status_pass(), check.label, "ok")
    } else {
        let detail = check
            .error
            .as_deref()
            .unwrap_or("check failed");
        format!("{}  {}  {}", status_fail(), check.label, detail)
    }
}

fn format_overall_readiness(ok: bool) -> String {
    if ok {
        format!("overall: {}  ready", status_pass())
    } else {
        format!("overall: {}  not ready", status_fail())
    }
}

fn status_pass() -> String {
    format!("{GREEN}✓ PASS{RESET}", GREEN = GREEN, RESET = RESET)
}

fn status_fail() -> String {
    format!("{RED}✗ FAILED{RESET}", RED = RED, RESET = RESET)
}

const GREEN: &str = "\x1b[32m";
const RED: &str = "\x1b[31m";
const RESET: &str = "\x1b[0m";

fn ms_file_count(assessment: &ConfigAssessment) -> usize {
    assessment.ms_files.as_ref().map_or(0, Vec::len)
}

fn ms_file_label(assessment: &ConfigAssessment) -> &'static str {
    if assessment.ms_files.is_some() {
        "resolved"
    } else {
        "None"
    }
}

fn fasta_entry_count(assessment: &ConfigAssessment) -> usize {
    assessment.fasta.as_ref().map_or(0, |db| db.entries.len())
}

fn glycan_entry_count(assessment: &ConfigAssessment) -> usize {
    assessment
        .glycan_library
        .as_ref()
        .map_or(0, |library| library.entries.len())
}

fn optional_label(present: Option<&str>) -> &'static str {
    if present.is_some() {
        "loaded"
    } else {
        "None"
    }
}

fn print_prefilter_summary(result: &PrefilterResult) {
    let stats = &result.stats;
    eprintln!("prefilter: scans={}", stats.scans_total);
    eprintln!("prefilter: diagnostic_positive={}", stats.diagnostic_positive);
    eprintln!("prefilter: isotope_pairs={}", stats.isotope_pairs);
    eprintln!("prefilter: filtered_scans={}", stats.filtered_scans);
    eprintln!("prefilter: rejected={}", stats.rejected);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU64, Ordering};

    static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

    fn temp_dir(name: &str) -> PathBuf {
        let id = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir().join(format!(
            "glycoquest_lib_test_{}_{}_{}",
            std::process::id(),
            name,
            id
        ))
    }

    fn write_xquest(root: &PathBuf) -> PathBuf {
        fs::create_dir_all(root).unwrap();
        let bin = root.join("xquest");
        fs::write(&bin, b"#!/bin/sh\n").unwrap();
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            fs::set_permissions(&bin, fs::Permissions::from_mode(0o755)).unwrap();
        }
        bin
    }

    fn write_fasta(path: &PathBuf) {
        fs::write(path, ">protein\nACDEFGHIK\n").unwrap();
    }

    fn fixture_mzxml(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/mzxml")
            .join(name)
    }

    fn test_config(input: PathBuf, out: PathBuf, dry_run: bool) -> RunConfig {
        let xquest_root = temp_dir("xquest");
        write_xquest(&xquest_root);
        let fasta = out.join("proteins.fasta");
        write_fasta(&fasta);
        RunConfig {
            cli: CliParams {
                input,
                database: fasta,
                glycans: "nglyc309".into(),
                xquest_root,
                out,
                dry_run,
                ..CliParams::default()
            },
            settings: Settings::defaults(),
            crosslinker: CrosslinkerProfile::resolve(&Settings::defaults(), Some("dss")).unwrap(),
        }
    }

    #[test]
    fn rejects_raw_vendor_input() {
        let out = temp_dir("raw_out");
        fs::create_dir_all(&out).unwrap();
        let raw_path = out.join("sample.raw");
        fs::write(&raw_path, b"").unwrap();
        let xquest_root = temp_dir("raw_xquest");
        write_xquest(&xquest_root);
        let fasta = out.join("proteins.fasta");
        write_fasta(&fasta);

        let cli = CliParams {
            input: raw_path,
            database: fasta,
            xquest_root,
            out: out.clone(),
            dry_run: true,
            ..CliParams::default()
        };
        let config = RunConfig {
            cli,
            settings: Settings::defaults(),
            crosslinker: CrosslinkerProfile::resolve(&Settings::defaults(), Some("dss")).unwrap(),
        };
        assert_eq!(run_config(&config), ExitCode::Validation as i32);
    }

    #[test]
    fn dry_run_writes_pipeline_outputs_for_dss_pair_fixture() {
        let out = temp_dir("dry_run_out");
        fs::create_dir_all(&out).unwrap();
        let config = test_config(fixture_mzxml("dss_pair.mzXML"), out.clone(), true);
        assert_eq!(run_config(&config), ExitCode::Success as i32);
        assert!(out.join("filtered_spectra.tsv").is_file());
        assert!(out.join("isotope_pairs.tsv").is_file());
        assert!(out.join("rejected_spectra.tsv").is_file());
        assert!(out.join("glycan_pruning.tsv").is_file());
        assert!(out.join("plan.json").is_file());
        assert!(out.join("spectra").is_dir());
        assert!(out.join("tmp/jobs").is_dir());
        let _ = fs::remove_dir_all(out);
    }

    #[test]
    fn dry_run_dmtmm_writes_jobs_without_isotope_pairs() {
        let out = temp_dir("dmtmm_out");
        fs::create_dir_all(&out).unwrap();
        let xquest_root = temp_dir("dmtmm_xquest");
        write_xquest(&xquest_root);
        let fasta = out.join("proteins.fasta");
        write_fasta(&fasta);
        let config = RunConfig {
            cli: CliParams {
                input: fixture_mzxml("hexnac_positive.mzXML"),
                database: fasta,
                glycans: "nglyc309".into(),
                crosslinker: Some("dmtmm".into()),
                xquest_root,
                out: out.clone(),
                dry_run: true,
                ..CliParams::default()
            },
            settings: Settings::defaults(),
            crosslinker: CrosslinkerProfile::resolve(&Settings::defaults(), Some("dmtmm"))
                .unwrap(),
        };
        assert_eq!(run_config(&config), ExitCode::Success as i32);
        let plan = fs::read_to_string(out.join("plan.json")).unwrap();
        assert!(plan.contains("\"isotope_prefilter_enabled\": false"));
        assert!(plan.contains("\"isotope_pairs\": 0"));
        let jobs_dir = fs::read_dir(out.join("tmp/jobs"))
            .unwrap()
            .next()
            .unwrap()
            .unwrap()
            .path();
        let xquest_def = fs::read_to_string(jobs_dir.join("xquest.def")).unwrap();
        assert!(xquest_def.contains("isotopeshift 0.0000000"));
        assert!(xquest_def.contains("AArequired K:E,K:D"));
        let matchlist = fs::read_to_string(jobs_dir.join("glycoquest_matched.txt")).unwrap();
        assert!(matchlist.contains("\tlight\tlight\t"));
        let _ = fs::remove_dir_all(out);
    }

    #[test]
    fn dry_run_exits_no_spectra_when_all_rejected() {
        let out = temp_dir("no_spectra_out");
        fs::create_dir_all(&out).unwrap();
        let config = test_config(fixture_mzxml("no_diagnostic.mzXML"), out.clone(), true);
        assert_eq!(run_config(&config), ExitCode::NoSpectra as i32);
        assert!(out.join("rejected_spectra.tsv").is_file());
        let _ = fs::remove_dir_all(out);
    }

    #[test]
    fn rejects_missing_fasta() {
        let out = temp_dir("missing_fasta_out");
        fs::create_dir_all(&out).unwrap();
        let xquest_root = temp_dir("missing_fasta_xquest");
        write_xquest(&xquest_root);
        let config = RunConfig {
            cli: CliParams {
                input: fixture_mzxml("dss_pair.mzXML"),
                database: out.join("missing.fasta"),
                xquest_root,
                out: out.clone(),
                dry_run: true,
                ..CliParams::default()
            },
            settings: Settings::defaults(),
            crosslinker: CrosslinkerProfile::resolve(&Settings::defaults(), Some("dss")).unwrap(),
        };
        assert_eq!(run_config(&config), ExitCode::Validation as i32);
        let _ = fs::remove_dir_all(out);
    }

    #[test]
    fn summary_shows_none_when_xquest_missing() {
        let out = temp_dir("summary_xquest_out");
        fs::create_dir_all(&out).unwrap();
        let fasta = out.join("proteins.fasta");
        write_fasta(&fasta);
        let config = RunConfig {
            cli: CliParams {
                input: fixture_mzxml("dss_pair.mzXML"),
                database: fasta,
                xquest_root: temp_dir("summary_missing_xquest"),
                out: out.clone(),
                dry_run: true,
                ..CliParams::default()
            },
            settings: Settings::defaults(),
            crosslinker: CrosslinkerProfile::resolve(&Settings::defaults(), Some("dss")).unwrap(),
        };

        let assessment = assess_config(&config);
        let summary = format_config_summary(&config, &assessment);
        assert!(summary.iter().any(|line| line.contains("xquest_bin: None")));
        assert!(summary.iter().any(|line| line.starts_with("MS files: 1")));
        assert!(summary.iter().any(|line| line.starts_with("FASTA entries: 1")));

        let readiness = format_readiness_report(&assessment);
        assert!(readiness.iter().any(|line| line.contains("xQuest runtime")));
        assert!(readiness.iter().any(|line| line.contains("FAILED")));
        assert!(!assessment.all_passed());

        assert_eq!(run_config(&config), ExitCode::Validation as i32);
        let _ = fs::remove_dir_all(out);
    }

    #[test]
    fn cli_ppm_overrides_settings() {
        let cli = CliParams {
            ppm_tolerance: Some(15.0),
            ..CliParams::default()
        };
        let config = RunConfig::load(cli).unwrap();
        assert_eq!(config.settings.diagnostic_tolerance_ppm, 15.0);
    }

    #[test]
    fn optional_hcg_dry_run_integration() {
        let mzxml = match std::env::var("GLYCOQUEST_HCG_MZXML") {
            Ok(path) => PathBuf::from(path),
            Err(_) => return,
        };
        let fasta = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("data/rcsb_pdb_1HRP_no_contams.fasta");
        if !mzxml.is_file() || !fasta.is_file() {
            return;
        }
        let xquest_root = std::env::var("GLYCOQUEST_XQUEST_ROOT")
            .map(PathBuf::from)
            .unwrap_or_else(|_| {
                PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("xQuest/V2.1.6/xquest")
            });
        if !xquest_root.is_dir() {
            return;
        }
        let out = temp_dir("hcg_integration");
        let _ = fs::create_dir_all(&out);
        let config = RunConfig {
            cli: CliParams {
                input: mzxml,
                database: fasta,
                xquest_root,
                out: out.clone(),
                dry_run: true,
                ..CliParams::default()
            },
            settings: Settings::defaults(),
            crosslinker: CrosslinkerProfile::resolve(&Settings::defaults(), Some("dss")).unwrap(),
        };
        assert_eq!(run_config(&config), ExitCode::Success as i32);
        let plan = fs::read_to_string(out.join("plan.json")).unwrap();
        assert!(plan.contains("\"scans_total\""));
        let _ = fs::remove_dir_all(out);
    }

    #[test]
    fn optional_xquest_integration_skipped_without_env() {
        if std::env::var("GLYCOQUEST_XQUEST_ROOT").is_err() {
            return;
        }
    }
}
