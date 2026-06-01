//! GlycoQuest library: CLI parameter types, settings, and the entry-point runner.

mod cli;
mod glycan;
mod input;
mod settings;

pub use cli::{parse_cli, CliParams};
pub use glycan::{
    glycan_data_dir, load_glycan_database, resolve_database, supported_glycan_databases,
    DiagnosticIon, GlycanEntry, GlycanLibrary,
};
pub use input::resolve_input;
pub use settings::{default_settings_path, Settings};

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
        if let Some(name) = &cli.crosslinker {
            settings.crosslinker_name = name.clone();
        }

        Ok(Self { cli, settings })
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

/// Validated MS input and optional parsed glycan library.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidatedInputs {
    pub files: Vec<std::path::PathBuf>,
    pub glycan_library: Option<GlycanLibrary>,
}

/// Run GlycoQuest from a fully merged [`RunConfig`].
pub fn run_config(config: &RunConfig) -> i32 {
    let validated = match validate_config(config) {
        Ok(validated) => validated,
        Err(message) => {
            eprintln!("error: {message}");
            return ExitCode::Validation.into();
        }
    };

    match config.execution_mode() {
        ExecutionMode::DryRun => {
            eprintln!("dry-run: configuration accepted (no search executed)");
            print_config_summary(config, &validated);
            ExitCode::Success.into()
        }
        ExecutionMode::Run => {
            eprintln!(
                "run: not implemented yet (input={})",
                config.cli.input.display()
            );
            print_config_summary(config, &validated);
            ExitCode::Success.into()
        }
    }
}

fn validate_config(config: &RunConfig) -> Result<ValidatedInputs, String> {
    let files = resolve_input(&config.cli.input)?;
    let glycan_library = Some(load_glycan_database(&config.cli.glycans)?);

    Ok(ValidatedInputs {
        files,
        glycan_library,
    })
}

fn print_config_summary(config: &RunConfig, validated: &ValidatedInputs) {
    let cli = &config.cli;
    let settings = &config.settings;
    eprintln!("input: {}", cli.input.display());
    eprintln!("MS files: {}", validated.files.len());
    for path in &validated.files {
        eprintln!("  {}", path.display());
    }
    eprintln!("database: {}", cli.database.display());
    eprintln!("glycans: {}", cli.glycans);
    if let Some(library) = &validated.glycan_library {
        eprintln!("glycan entries: {} (unique)", library.entries.len());
        for entry in library.entries.iter().take(3) {
            eprintln!(
                "  {}  mass={}  diagnostics={}",
                entry.composition,
                entry.monoisotopic_mass,
                entry.diagnostic_ions.len()
            );
        }
        if library.entries.len() > 3 {
            eprintln!("  …");
        }
    }
    eprintln!("xquest_root: {}", cli.xquest_root.display());
    if let Some(path) = &settings.xquest_bin {
        eprintln!("xquest_bin (settings): {}", path.display());
    }
    eprintln!("out: {}", cli.out.display());
    eprintln!("crosslinker: {}", settings.crosslinker_name);
    eprintln!(
        "diagnostic_tolerance_ppm: {}",
        settings.diagnostic_tolerance_ppm
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU64, Ordering};

    static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

    fn temp_mzxml(name: &str) -> PathBuf {
        let id = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "glycoquest_lib_test_{}_{}_{}.mzXML",
            std::process::id(),
            name,
            id
        ));
        fs::write(&path, b"").unwrap();
        path
    }

    #[test]
    fn rejects_raw_vendor_input() {
        let raw_path = std::env::temp_dir().join(format!(
            "glycoquest_lib_test_{}_sample.raw",
            std::process::id()
        ));
        fs::write(&raw_path, b"").unwrap();

        let cli = CliParams {
            input: raw_path,
            ..CliParams::default()
        };
        let config = RunConfig {
            cli,
            settings: Settings::defaults(),
        };
        assert_eq!(run_config(&config), ExitCode::Validation as i32);
    }

    #[test]
    fn dry_run_accepts_bundled_glycan_database() {
        let cli = CliParams {
            input: temp_mzxml("glycan_dry_run"),
            glycans: "nglyc309".into(),
            dry_run: true,
            ..CliParams::default()
        };
        let config = RunConfig {
            cli,
            settings: Settings::defaults(),
        };
        assert_eq!(run_config(&config), ExitCode::Success as i32);
    }

    #[test]
    fn dry_run_accepts_minimal_cli_params() {
        let cli = CliParams {
            input: temp_mzxml("dry_run"),
            dry_run: true,
            ..CliParams::default()
        };
        let config = RunConfig {
            cli,
            settings: Settings::defaults(),
        };
        assert_eq!(run_config(&config), ExitCode::Success as i32);
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
}
