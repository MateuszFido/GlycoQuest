// Copyright (c) ETH Zurich, Mateusz Fido

//! Input path resolution: one xQuest-compatible mzXML file or a directory of them.

use std::path::{Path, PathBuf};

const RAW_VENDOR_EXTENSIONS: &[&str] = &["raw", "wiff", "d", "baf", "tdf"];

/// True when `path` uses a raw vendor extension (Thermo `.raw`, Sciex `.wiff`, etc.).
pub fn is_raw(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| {
            let lower = ext.to_ascii_lowercase();
            RAW_VENDOR_EXTENSIONS.contains(&lower.as_str())
        })
        .unwrap_or(false)
}

/// True when `path` looks like an mzXML file by extension (case-insensitive).
pub fn is_mzxml(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext.eq_ignore_ascii_case("mzxml"))
        .unwrap_or(false)
}

/// Resolve the CLI `INPUT` path to one or more existing mzXML files.
///
/// Accepts a single `.mzXML`/`.mzxml` file or a directory containing such files
/// (non-recursive, only immediate children scanned)
pub fn resolve_input(path: &Path) -> Result<Vec<PathBuf>, String> {
    if is_raw(path) {
        return Err(raw_vendor_error(path));
    }

    if !path.exists() {
        return Err(format!("MS file input not found: {}", path.display()));
    }

    if path.is_file() {
        return resolve_files(&[path]);
    }

    if path.is_dir() {
        return resolve_dir(path);
    }

    Err(format!(
        "MS input is neither a file nor a directory: {}",
        path.display()
    ))
}

/// Validate one or more explicit mzXML file paths.
///
/// Pass a single path as `resolve_files(&[path])` or `resolve_files(std::slice::from_ref(path))`.
pub fn resolve_files(paths: &[&Path]) -> Result<Vec<PathBuf>, String> {
    if paths.is_empty() {
        return Err("MS/MS input: at least one file path is required".to_string());
    }

    let mut files = Vec::with_capacity(paths.len());
    for path in paths {
        if is_raw(path) {
            return Err(raw_vendor_error(path));
        }
        if !path.exists() {
            return Err(format!("MS file input not found: {}", path.display()));
        }
        if !path.is_file() {
            return Err(format!(
                "MS input must be an mzXML file; got non-file: {}",
                path.display()
            ));
        }
        if !is_mzxml(path) {
            return Err(format!(
                "MS input must be an mzXML file or a directory of mzXML files; got file: {}",
                path.display()
            ));
        }
        files.push(path.to_path_buf());
    }

    files.sort();
    Ok(files)
}

fn resolve_dir(dir: &Path) -> Result<Vec<PathBuf>, String> {
    let entries = std::fs::read_dir(dir)
        .map_err(|err| format!("cannot read MS input directory {}: {err}", dir.display()))?;

    let files: Vec<PathBuf> = entries
        .filter_map(|entry| entry.ok())
        .map(|entry| entry.path())
        .filter(|entry_path| entry_path.is_file() && is_mzxml(entry_path))
        .collect();

    if files.is_empty() {
        return Err(format!(
            "no mzXML files found in input directory: {}",
            dir.display()
        ));
    }

    let path_refs: Vec<&Path> = files.iter().map(PathBuf::as_path).collect();
    resolve_files(&path_refs)
}

fn raw_vendor_error(path: &Path) -> String {
    format!(
        "Unsupported raw vendor input: {}. Convert explicitly to mzXML before running GlycoQuest.",
        path.display()
    )
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
            "glycoquest_input_test_{}_{}_{}",
            std::process::id(),
            name,
            id
        ))
    }

    fn write_empty_file(path: &Path) {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).unwrap();
        }
        fs::write(path, b"").unwrap();
    }

    #[test]
    fn accepts_single_mzxml_file() {
        let path = temp_dir("single").join("run.mzXML");
        write_empty_file(&path);

        let resolved = resolve_input(&path).unwrap();
        assert_eq!(resolved, vec![path]);
    }

    #[test]
    fn accepts_mzxml_extension_case_insensitively() {
        let path = temp_dir("case").join("run.mzxml");
        write_empty_file(&path);

        let resolved = resolve_input(&path).unwrap();
        assert_eq!(resolved, vec![path]);
    }

    #[test]
    fn accepts_directory_of_mzxml_files() {
        let dir = temp_dir("dir");
        fs::create_dir_all(&dir).unwrap();
        let b = dir.join("b.mzXML");
        let a = dir.join("a.mzXML");
        write_empty_file(&b);
        write_empty_file(&a);

        let resolved = resolve_input(&dir).unwrap();
        assert_eq!(resolved, vec![a, b]);
    }

    #[test]
    fn accepts_multiple_explicit_files() {
        let dir = temp_dir("multi");
        fs::create_dir_all(&dir).unwrap();
        let a = dir.join("a.mzXML");
        let b = dir.join("b.mzXML");
        write_empty_file(&a);
        write_empty_file(&b);

        let resolved = resolve_files(&[a.as_path(), b.as_path()]).unwrap();
        assert_eq!(resolved, vec![a, b]);
    }

    #[test]
    fn rejects_missing_path() {
        let path = temp_dir("missing").join("nope.mzXML");
        let err = resolve_input(&path).unwrap_err();
        assert!(err.contains("not found"));
    }

    #[test]
    fn rejects_raw_vendor_file() {
        let path = temp_dir("raw").join("sample.raw");
        write_empty_file(&path);

        let err = resolve_input(&path).unwrap_err();
        assert!(err.contains("Unsupported raw vendor input"));
    }

    #[test]
    fn rejects_non_mzxml_file() {
        let path = temp_dir("fasta").join("proteins.fasta");
        write_empty_file(&path);

        let err = resolve_input(&path).unwrap_err();
        assert!(err.contains("must be an mzXML file"));
    }

    #[test]
    fn rejects_empty_directory() {
        let dir = temp_dir("empty");
        fs::create_dir_all(&dir).unwrap();

        let err = resolve_input(&dir).unwrap_err();
        assert!(err.contains("no mzXML files found"));
    }
}
