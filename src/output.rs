// Copyright (c) ETH Zurich, Mateusz Fido

//! Output directory preparation and run plan artifacts.

pub mod layout;

mod plan;

use std::path::Path;

pub use layout::{
    DEFAULT_OUT_BASE, ProjectLayout, cleanup_temp_artifacts, derive_project_slug,
    resolve_project_out_dir,
};
pub use plan::{RunPlanDocument, build_run_plan_document, write_plan_json};

/// Create the output directory if missing, or verify it is a writable directory.
pub fn ensure_output_dir(path: &Path) -> Result<(), String> {
    if path.exists() {
        if !path.is_dir() {
            return Err(format!(
                "output path exists but is not a directory: {}",
                path.display()
            ));
        }
    } else {
        std::fs::create_dir_all(path)
            .map_err(|err| format!("cannot create output directory {}: {err}", path.display()))?;
    }

    let test_file = path.join(".glycoquest_write_test");
    std::fs::write(&test_file, b"").map_err(|err| {
        format!(
            "output directory is not writable: {} ({err})",
            path.display()
        )
    })?;
    let _ = std::fs::remove_file(test_file);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::sync::atomic::{AtomicU64, Ordering};

    static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

    fn temp_out(name: &str) -> std::path::PathBuf {
        let id = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir().join(format!(
            "glycoquest_out_test_{}_{}_{}",
            std::process::id(),
            name,
            id
        ))
    }

    #[test]
    fn creates_missing_directory() {
        let path = temp_out("create");
        assert!(!path.exists());
        ensure_output_dir(&path).unwrap();
        assert!(path.is_dir());
        let _ = fs::remove_dir_all(path);
    }

    #[test]
    fn accepts_existing_writable_directory() {
        let path = temp_out("existing");
        fs::create_dir_all(&path).unwrap();
        ensure_output_dir(&path).unwrap();
        let _ = fs::remove_dir_all(path);
    }

    #[test]
    fn rejects_file_path() {
        let path = temp_out("file");
        fs::write(&path, b"x").unwrap();
        let err = ensure_output_dir(&path).unwrap_err();
        assert!(err.contains("not a directory"));
        let _ = fs::remove_file(path);
    }
}
