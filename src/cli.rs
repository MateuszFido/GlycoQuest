//! Command-line parsing for the `glycoquest` binary.

use std::path::PathBuf;

use clap::error::ErrorKind;
use clap::Parser;

/// Parsed GlycoQuest CLI parameters 
#[derive(Debug, Clone, Default, PartialEq)]
pub struct CliParams {
    pub input: PathBuf,
    pub database: Option<PathBuf>,
    pub glycans: Option<String>,
    pub xquest_root: Option<PathBuf>,
    pub crosslinker: Option<String>,
    pub ppm_tolerance: Option<f64>,
    pub out: Option<PathBuf>,
    pub config: Option<PathBuf>,
    /// When true, validate configuration only 
    pub dry_run: bool,
}

#[derive(Parser, Debug)]
#[command(
    name = "glycoquest",
    version,
    about = "Prepare and run xQuest searches for DSS-crosslinked glycopeptide-peptide spectra.",
    arg_required_else_help = false,
    after_help = "Advanced options (xquest_bin, tolerances, modifications, limits, etc.) live in settings.ini.\n\nExamples:\n  glycoquest input.mzXML --database proteins.fasta --glycans nglyc309 --xquest-root ./xquest --out job\n  glycoquest input.mzXML --database proteins.fasta --glycans nglyc309 --out job --dry-run"
)]
struct Args {
    /// mzXML file or directory of MS/MS inputs (xQuest-compatible).
    #[arg(value_name = "INPUT")]
    input: Option<PathBuf>,

    /// Protein sequence database (FASTA).
    #[arg(long, value_name = "FASTA")]
    database: Option<PathBuf>,

    /// Bundled glycan database to load (e.g. nglyc309, oglyc78).
    #[arg(long, value_name = "DATABASE")]
    glycans: Option<String>,

    /// xQuest installation root (contains `xquest.def` templates).
    #[arg(long, value_name = "DIR")]
    xquest_root: Option<PathBuf>,

    /// Crosslinker chemistry name (e.g. dss). Overrides settings.ini [crosslinker] name.
    #[arg(long, value_name = "NAME")]
    crosslinker: Option<String>,

    /// Diagnostic-ion matching tolerance in ppm. Overrides settings.ini diagnostic_tolerance_ppm.
    #[arg(long, value_name = "PPM")]
    ppm_tolerance: Option<f64>,

    /// Output directory for generated jobs and results.
    #[arg(long, value_name = "DIR")]
    out: Option<PathBuf>,

    /// Path to settings.ini (default: ./settings.ini).
    #[arg(long, value_name = "FILE", default_value = "settings.ini")]
    config: PathBuf,

    /// Validate inputs and print a summary without running xQuest.
    #[arg(long)]
    dry_run: bool,
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
        out: args.out,
        config: Some(args.config),
        dry_run: args.dry_run,
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
         glycoquest data/run.mzXML \\\n    \
           --database proteins.fasta \\\n    \
           --glycans nglyc309 \\\n    \
           --xquest-root ./xquest \\\n    \
           --out results\n\
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
        let err = parse_cli(["glycoquest"]).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("missing required argument: INPUT"));
        assert!(msg.contains("mzXML"));
        assert!(msg.contains("glycoquest data/run.mzXML"));
    }

    #[test]
    fn parses_dry_run_flag() {
        let params = parse_cli([
            "glycoquest",
            "input.mzXML",
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
            "--glycans",
            "nglyc309",
            "--config",
            "settings.ini",
        ])
        .unwrap();
        assert_eq!(params.glycans.as_deref(), Some("nglyc309"));
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
        assert_eq!(
            params.database.as_deref(),
            Some(std::path::Path::new("proteins.fasta"))
        );
    }
}
