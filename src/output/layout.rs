//! Project output directory layout under a single `out/` base.

use std::path::{Path, PathBuf};

use crate::cli::input::{is_mzxml, resolve_input};

/// Default `--out` value: project runs live under `out/<project>/`.
pub const DEFAULT_OUT_BASE: &str = "out";

/// Resolved paths for one GlycoQuest project run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProjectLayout {
    pub root: PathBuf,
}

impl ProjectLayout {
    pub fn new(root: PathBuf) -> Self {
        Self { root }
    }

    /// Truncated mzXML containing only prefilter-retained MS2 scans (user-inspectable).
    pub fn spectra_dir(&self) -> PathBuf {
        self.root.join("spectra")
    }

    /// Ephemeral xQuest working files (jobs, logs, peptide DB indexes).
    pub fn tmp_dir(&self) -> PathBuf {
        self.root.join("tmp")
    }

    pub fn jobs_dir(&self) -> PathBuf {
        self.tmp_dir().join("jobs")
    }

    pub fn logs_dir(&self) -> PathBuf {
        self.tmp_dir().join("logs")
    }

    pub fn results_dir(&self) -> PathBuf {
        self.root.join("results")
    }

    /// Interactive CLMS viewer bundle (`viewer.json` + static assets).
    pub fn viewer_dir(&self) -> PathBuf {
        self.results_dir().join("viewer")
    }
}

/// Resolve the project output root.
///
/// When `--out` is the default base (`out`), append a slug derived from the first
/// mzXML input file (or the FASTA basename as fallback). Any other `--out` value is
/// treated as an explicit project directory (backward compatible).
pub fn resolve_project_out_dir(
    out_flag: &Path,
    input: &Path,
    database: &Path,
) -> Result<PathBuf, String> {
    if out_flag != Path::new(DEFAULT_OUT_BASE) {
        return Ok(out_flag.to_path_buf());
    }
    let slug = derive_project_slug(input, database)?;
    Ok(PathBuf::from(DEFAULT_OUT_BASE).join(slug))
}

pub fn derive_project_slug(input: &Path, database: &Path) -> Result<String, String> {
    if let Ok(files) = resolve_input(input) {
        if let Some(first) = files.first() {
            if let Some(slug) = slug_from_path(first) {
                return Ok(slug);
            }
        }
    } else if input.is_file() && is_mzxml(input) {
        if let Some(slug) = slug_from_path(input) {
            return Ok(slug);
        }
    }

    database
        .file_stem()
        .and_then(|s| s.to_str())
        .map(sanitize_slug)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| {
            format!(
                "cannot derive project name from input {} or database {}",
                input.display(),
                database.display()
            )
        })
}

fn slug_from_path(path: &Path) -> Option<String> {
    path.file_stem()
        .and_then(|s| s.to_str())
        .map(sanitize_slug)
        .filter(|s| !s.is_empty())
}

fn sanitize_slug(raw: &str) -> String {
    let mut slug = String::with_capacity(raw.len());
    let mut prev_sep = false;
    for ch in raw.chars() {
        if ch.is_ascii_alphanumeric() {
            slug.push(ch.to_ascii_lowercase());
            prev_sep = false;
        } else if !prev_sep {
            slug.push('_');
            prev_sep = true;
        }
    }
    let trimmed = slug.trim_matches('_');
    if trimmed.is_empty() {
        "run".to_string()
    } else {
        trimmed.to_string()
    }
}

/// Remove ephemeral `tmp/` after a successful run.
///
/// Set `GLYCOQUEST_KEEP_TMP=1` to retain the xQuest job folders, logs, and
/// per-job result XML for debugging.
pub fn cleanup_temp_artifacts(layout: &ProjectLayout) {
    if std::env::var_os("GLYCOQUEST_KEEP_TMP").is_some() {
        return;
    }
    if layout.tmp_dir().is_dir() {
        let _ = std::fs::remove_dir_all(layout.tmp_dir());
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
            "glycoquest_layout_test_{}_{}_{}",
            std::process::id(),
            name,
            id
        ))
    }

    #[test]
    fn derives_project_slug_from_mzxml_stem() {
        let dir = temp_dir("slug");
        fs::create_dir_all(&dir).unwrap();
        let mz = dir.join("260607_LU02_disoic_ASF_DSS_1.c.mzXML");
        fs::write(&mz, b"").unwrap();
        let slug = derive_project_slug(&mz, Path::new("target_proteins_asf.fasta")).unwrap();
        assert_eq!(slug, "260607_lu02_disoic_asf_dss_1_c");
    }

    #[test]
    fn default_out_appends_project_slug() {
        let dir = temp_dir("input");
        fs::create_dir_all(&dir).unwrap();
        let mz = dir.join("sample_A.mzXML");
        fs::write(&mz, b"").unwrap();
        let out = resolve_project_out_dir(
            Path::new(DEFAULT_OUT_BASE),
            &dir,
            Path::new("proteins.fasta"),
        )
        .unwrap();
        assert_eq!(out, PathBuf::from("out/sample_a"));
    }

    #[test]
    fn explicit_out_is_used_as_project_root() {
        let out = resolve_project_out_dir(
            Path::new("glycoquest_custom_out"),
            Path::new("run.mzXML"),
            Path::new("proteins.fasta"),
        )
        .unwrap();
        assert_eq!(out, PathBuf::from("glycoquest_custom_out"));
    }
}
