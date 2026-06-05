//! xQuest runtime discovery and validation.

use std::path::{Path, PathBuf};

use crate::cli::settings::Settings;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct XQuestRuntime {
    pub executable: PathBuf,
    pub root: PathBuf,
}

/// Resolve the xQuest executable from settings or `--xquest-root`.
pub fn resolve_runtime(xquest_root: &Path, settings: &Settings) -> Result<XQuestRuntime, String> {
    if let Some(bin) = &settings.xquest_bin {
        return validate_executable(bin, xquest_root);
    }

    let candidates = candidate_paths(xquest_root);
    for candidate in candidates {
        if candidate.is_file() {
            return validate_executable(&candidate, xquest_root);
        }
    }

    Err(format!(
        "xQuest executable not found under {} \
         (set xquest_bin in settings.ini or pass --xquest-root); \
         checked: xquest, bin/xquest, xquest/xquest",
        xquest_root.display()
    ))
}

fn candidate_paths(root: &Path) -> Vec<PathBuf> {
    vec![
        root.join("xquest"),
        root.join("bin").join("xquest"),
        root.join("xquest").join("xquest"),
    ]
}

fn validate_executable(path: &Path, root: &Path) -> Result<XQuestRuntime, String> {
    if !path.is_file() {
        return Err(format!("xQuest executable not found: {}", path.display()));
    }

    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let meta = std::fs::metadata(path).map_err(|err| {
            format!("cannot read xQuest executable {}: {err}", path.display())
        })?;
        if meta.permissions().mode() & 0o111 == 0 {
            return Err(format!(
                "xQuest executable is not executable: {}",
                path.display()
            ));
        }
    }

    Ok(XQuestRuntime {
        executable: path.to_path_buf(),
        root: root.to_path_buf(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::sync::atomic::{AtomicU64, Ordering};

    static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

    fn temp_xquest_root(name: &str) -> PathBuf {
        let id = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir().join(format!(
            "glycoquest_xquest_test_{}_{}_{}",
            std::process::id(),
            name,
            id
        ))
    }

    #[test]
    fn resolves_from_settings_xquest_bin() {
        let root = temp_xquest_root("settings");
        fs::create_dir_all(&root).unwrap();
        let bin = root.join("custom_xquest");
        fs::write(&bin, b"#!/bin/sh\n").unwrap();
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            fs::set_permissions(&bin, fs::Permissions::from_mode(0o755)).unwrap();
        }

        let settings = Settings {
            xquest_bin: Some(bin.clone()),
            ..Settings::defaults()
        };
        let runtime = resolve_runtime(&root, &settings).unwrap();
        assert_eq!(runtime.executable, bin);
    }

    #[test]
    fn resolves_xquest_under_root() {
        let root = temp_xquest_root("under_root");
        let bin = root.join("xquest");
        fs::create_dir_all(&root).unwrap();
        fs::write(&bin, b"#!/bin/sh\n").unwrap();
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            fs::set_permissions(&bin, fs::Permissions::from_mode(0o755)).unwrap();
        }

        let settings = Settings::defaults();
        let runtime = resolve_runtime(&root, &settings).unwrap();
        assert_eq!(runtime.executable, bin);
    }

    #[test]
    fn fails_when_executable_missing() {
        let root = temp_xquest_root("missing");
        fs::create_dir_all(&root).unwrap();
        let settings = Settings::defaults();
        let err = resolve_runtime(&root, &settings).unwrap_err();
        assert!(err.contains("xQuest executable not found"));
    }
}
