//! Command-line parsing for the `glycoquest` binary.

pub mod input;
pub mod settings;
use std::path::PathBuf;

use clap::error::ErrorKind;
use clap::Parser;

/// Default bundled glycan database id ([`parse_cli`]).
pub const DEFAULT_GLYCANS: &str = "nglyc309";

/// Parsed GlycoQuest CLI parameters
#[derive(Debug, Clone, PartialEq)]
pub struct CliParams {
    pub input: PathBuf,
    pub database: PathBuf,
    pub glycans: String,
    pub xquest_root: PathBuf,
    pub crosslinker: Option<String>,
    pub ppm_tolerance: Option<f64>,
    /// Number of xQuest jobs to run concurrently (overrides settings.ini; 0 = one per CPU core).
    pub jobs: Option<u32>,
    pub out: PathBuf,
    pub config: Option<PathBuf>,
    /// When true, validate configuration only
    pub dry_run: bool,
    /// Resume xQuest jobs in an existing project output directory, skipping jobs with results.
    pub resume: bool,
}

impl Default for CliParams {
    fn default() -> Self {
        Self {
            input: PathBuf::from("."),
            database: PathBuf::from("proteins.fasta"),
            glycans: DEFAULT_GLYCANS.into(),
            xquest_root: PathBuf::from("."),
            crosslinker: None,
            ppm_tolerance: None,
            jobs: None,
            out: PathBuf::from(crate::output::DEFAULT_OUT_BASE),
            config: None,
            dry_run: false,
            resume: false,
        }
    }
}

#[derive(Parser, Debug)]
#[command(
    name = "glycoquest",
    version,
    about = "Prepare and run xQuest searches for DSS-crosslinked glycopeptide-peptide spectra.",
    arg_required_else_help = false,
    after_help = "Advanced options (xquest_bin, tolerances, modifications, limits, etc.) live in settings.ini.\n\nExamples:\n  glycoquest input.mzXML --database proteins.fasta\n  glycoquest input.mzXML --database proteins.fasta --out job --dry-run"
)]
struct Args {
    /// mzXML file(s) or directory (xQuest-compatible).
    #[arg(value_name = "INPUT")]
    input: Option<PathBuf>,

    /// Protein sequence database (FASTA).
    #[arg(long, required = true, value_name = "FASTA")]
    database: PathBuf,

    /// Glycan library: a bundled database id (e.g. nglyc309, oglyc78) or a path to a
    /// CSV/TSV file (columns: name,composition,monoisotopic_mass,diagnostic_ions,residue_targets).
    #[arg(long, value_name = "DATABASE|FILE", default_value = DEFAULT_GLYCANS)]
    glycans: String,

    /// xQuest installation root (contains `xquest.def` templates).
    #[arg(long, value_name = "DIR", default_value = ".")]
    xquest_root: PathBuf,

    /// Crosslinker chemistry name (e.g. dss). Overrides settings.ini [crosslinker] name.
    #[arg(long, value_name = "NAME")]
    crosslinker: Option<String>,

    /// Diagnostic-ion matching tolerance in ppm. Overrides settings.ini diagnostic_tolerance_ppm.
    #[arg(long, value_name = "PPM")]
    ppm_tolerance: Option<f64>,

    /// Number of xQuest jobs to run concurrently. Overrides settings.ini [execution] job_parallelism.
    /// 0 = one per available CPU core.
    #[arg(long, short = 'j', value_name = "N")]
    jobs: Option<u32>,

    /// Output base directory. Default `out` creates `out/<project>/` from the first input file.
    #[arg(long, value_name = "DIR", default_value = crate::output::DEFAULT_OUT_BASE)]
    out: PathBuf,

    /// Path to settings.ini (default: ./settings.ini).
    #[arg(long, value_name = "FILE", default_value = "settings.ini")]
    config: PathBuf,

    /// Validate inputs and print a summary without running xQuest.
    #[arg(long, conflicts_with = "resume")]
    dry_run: bool,

    /// Resume pending xQuest jobs in an existing --out project directory.
    #[arg(long, conflicts_with = "dry_run")]
    resume: bool,
}

/// Parse `std::env::args()` (skipping the program name) into [`CliParams`].
pub fn parse_cli<I, T>(args: I) -> Result<CliParams, clap::Error>
where
    I: IntoIterator<Item = T>,
    T: Into<std::ffi::OsString> + Clone,
{
    let args = match Args::try_parse_from(args) {
        Err(err) if passes_through_unchanged(&err) => return Err(err),
        Err(err) => return Err(explain_cli_error(err)),
        Ok(args) => args,
    };

    let Some(input) = args.input else {
        return Err(missing_input_error());
    };

    Ok(CliParams {
        input,
        database: args.database,
        glycans: args.glycans,
        xquest_root: args.xquest_root,
        crosslinker: args.crosslinker,
        ppm_tolerance: args.ppm_tolerance,
        jobs: args.jobs,
        out: args.out,
        config: Some(args.config),
        dry_run: args.dry_run,
        resume: args.resume,
    })
}

fn passes_through_unchanged(err: &clap::Error) -> bool {
    matches!(
        err.kind(),
        ErrorKind::DisplayHelp | ErrorKind::DisplayVersion
    )
}

fn missing_input_error() -> clap::Error {
    clap::Error::raw(
        ErrorKind::MissingRequiredArgument,
        "missing required argument: INPUT\n\
         \n\
         Provide one xQuest-compatible mzXML file or a directory of mzXML files \
         as the first (positional) argument.\n\
         \n\
         Example:\n  \
         glycoquest data/run.mzXML --database proteins.fasta\n\
         \n\
         Run `glycoquest --help` for all options.",
    )
}

fn explain_cli_error(err: clap::Error) -> clap::Error {
    match err.kind() {
        ErrorKind::UnknownArgument => clap::Error::raw(
            ErrorKind::UnknownArgument,
            format!(
                "{}\n\nRun `glycoquest --help` for supported options.",
                err.to_string().trim()
            ),
        ),
        ErrorKind::MissingRequiredArgument => clap::Error::raw(
            ErrorKind::MissingRequiredArgument,
            format!(
                "{}\n\nRun `glycoquest --help` for usage.",
                err.to_string().trim()
            ),
        ),
        _ => err,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn missing_input_error_is_explicit() {
        let err = parse_cli(["glycoquest", "--database", "proteins.fasta"]).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("missing required argument: INPUT"));
        assert!(msg.contains("mzXML"));
        assert!(msg.contains("glycoquest data/run.mzXML"));
    }

    #[test]
    fn missing_database_error() {
        let err = parse_cli(["glycoquest", "input.mzXML"]).unwrap_err();
        assert!(err.to_string().contains("--database"));
    }

    #[test]
    fn parses_dry_run_flag() {
        let params = parse_cli([
            "glycoquest",
            "input.mzXML",
            "--database",
            "proteins.fasta",
            "--dry-run",
            "--config",
            "settings.ini",
        ])
        .unwrap();
        assert!(params.dry_run);
    }

    #[test]
    fn parses_glycan_database_flag() {
        let params = parse_cli([
            "glycoquest",
            "input.mzXML",
            "--database",
            "proteins.fasta",
            "--glycans",
            "oglyc78",
            "--config",
            "settings.ini",
        ])
        .unwrap();
        assert_eq!(params.glycans, "oglyc78");
    }

    #[test]
    fn applies_defaults_for_glycans_out_and_xquest_root() {
        let params = parse_cli([
            "glycoquest",
            "input.mzXML",
            "--database",
            "proteins.fasta",
        ])
        .unwrap();
        assert_eq!(params.input, PathBuf::from("input.mzXML"));
        assert_eq!(params.database, PathBuf::from("proteins.fasta"));
        assert_eq!(params.glycans, DEFAULT_GLYCANS);
        assert_eq!(params.out, PathBuf::from(crate::output::DEFAULT_OUT_BASE));
        assert_eq!(params.xquest_root, PathBuf::from("."));
    }

    #[test]
    fn parses_jobs_flag() {
        let params = parse_cli([
            "glycoquest",
            "input.mzXML",
            "--database",
            "proteins.fasta",
            "--jobs",
            "8",
        ])
        .unwrap();
        assert_eq!(params.jobs, Some(8));

        let short = parse_cli([
            "glycoquest",
            "input.mzXML",
            "--database",
            "proteins.fasta",
            "-j",
            "4",
        ])
        .unwrap();
        assert_eq!(short.jobs, Some(4));
    }

    #[test]
    fn jobs_flag_defaults_to_none() {
        let params = parse_cli(["glycoquest", "input.mzXML", "--database", "proteins.fasta"])
            .unwrap();
        assert_eq!(params.jobs, None);
    }

    #[test]
    fn parses_positional_input_and_database() {
        let params = parse_cli([
            "glycoquest",
            "input.mzXML",
            "--database",
            "proteins.fasta",
            "--config",
            "settings.ini",
        ])
        .unwrap();
        assert_eq!(params.input, PathBuf::from("input.mzXML"));
        assert_eq!(params.database, PathBuf::from("proteins.fasta"));
    }
}
