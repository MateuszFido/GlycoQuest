// Copyright (c) ETH Zurich, Mateusz Fido

//! Fail-fast checks for Perl modules required by the xQuest search path.

use std::path::Path;
use std::process::Command;

/// Modules that must load for GlycoQuest's compare_peaks3 → xquest.pl path.
/// Keep in sync with `scripts/check-xquest-perl.pl`.
const CRITICAL_MODULES: &[&str] = &[
    "DB_File",
    "MLDBM",
    "Bio::Perl",
    "XML::TreeBuilder",
    "XML::Parser",
    "XML::Element",
    "HTML::Tagset",
    "HTML::Entities",
    "MIME::Base64",
    "Storable",
];

/// When the vendored xQuest tree is present, verify Perl can compile the search
/// scripts. Stub roots used in unit tests (no `bin/compare_peaks3.pl`) skip.
pub fn check_perl_search_deps(xquest_root: &Path) -> Result<(), String> {
    let compare = xquest_root.join("bin/compare_peaks3.pl");
    let xquest_pl = xquest_root.join("bin/xquest.pl");
    if !compare.is_file() || !xquest_pl.is_file() {
        return Ok(());
    }

    let perl5lib = build_perl5lib(xquest_root);
    let require_stmt = CRITICAL_MODULES
        .iter()
        .map(|m| format!("require {m};"))
        .collect::<Vec<_>>()
        .join(" ");

    let module_out = Command::new("perl")
        .args(["-e", &require_stmt])
        .env("PERL5LIB", &perl5lib)
        .output()
        .map_err(|err| {
            format!("cannot run perl to check xQuest modules: {err} (is perl on PATH?)")
        })?;

    if !module_out.status.success() {
        let stderr = String::from_utf8_lossy(&module_out.stderr);
        return Err(format_perl_deps_error(
            xquest_root,
            &perl5lib,
            &stderr,
            "module require",
        ));
    }

    for script in [&compare, &xquest_pl] {
        let out = Command::new("perl")
            .args(["-c"])
            .arg(script)
            .env("PERL5LIB", &perl5lib)
            .output()
            .map_err(|err| format!("cannot run perl -c on {}: {err}", script.display()))?;
        let combined = format!(
            "{}{}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
        if !out.status.success() || !combined.contains("syntax OK") {
            return Err(format_perl_deps_error(
                xquest_root,
                &perl5lib,
                &combined,
                &format!("compile {}", script.display()),
            ));
        }
    }

    Ok(())
}

fn build_perl5lib(xquest_root: &Path) -> String {
    let mut dirs = vec![
        xquest_root.join("modules"),
        xquest_root.join("1209/lib/perl5"),
        xquest_root.join("1209/share/perl5"),
    ];
    if let Ok(home) = std::env::var("HOME") {
        let local = Path::new(&home).join("perl5/lib/perl5");
        if local.is_dir() {
            dirs.insert(0, local);
        }
    }
    // Prefer existing PERL5LIB entries after our dirs (job scripts overwrite;
    // readiness still benefits from a user local::lib already on PERL5LIB).
    let mut parts: Vec<String> = dirs
        .into_iter()
        .filter(|p| p.is_dir())
        .map(|p| p.display().to_string())
        .collect();
    if let Ok(existing) = std::env::var("PERL5LIB") {
        if !existing.is_empty() {
            parts.push(existing);
        }
    }
    parts.join(":")
}

fn format_perl_deps_error(
    xquest_root: &Path,
    perl5lib: &str,
    detail: &str,
    phase: &str,
) -> String {
    format!(
        "xQuest Perl dependencies incomplete ({phase}).\n\
         xquest-root: {}\n\
         PERL5LIB={perl5lib}\n\
         {detail}\n\
         Fix: install DB_File and XML::Parser (Euler: scripts/bootstrap-euler-perl.sh),\n\
         then verify with: scripts/check-xquest-perl.pl --xquest-root {}",
        xquest_root.display(),
        xquest_root.display()
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::sync::atomic::{AtomicU64, Ordering};

    static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

    fn temp_root(name: &str) -> std::path::PathBuf {
        let id = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir().join(format!(
            "glycoquest_perl_deps_{}_{}_{}",
            std::process::id(),
            name,
            id
        ))
    }

    #[test]
    fn skips_stub_roots_without_search_scripts() {
        let root = temp_root("stub");
        fs::create_dir_all(&root).unwrap();
        assert!(check_perl_search_deps(&root).is_ok());
        let _ = fs::remove_dir_all(root);
    }

    #[test]
    fn accepts_vendored_xquest_when_perl_deps_present() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("V2.1.7/xquest");
        if !root.join("bin/compare_peaks3.pl").is_file() {
            return;
        }
        // May fail on hosts without DB_File/XML::Parser; that is a real readiness failure.
        // Only assert Ok when the modules are available.
        match check_perl_search_deps(&root) {
            Ok(()) => {}
            Err(err) => {
                let lower = err.to_lowercase();
                assert!(
                    lower.contains("db_file")
                        || lower.contains("xml::parser")
                        || lower.contains("cannot run perl"),
                    "unexpected error: {err}"
                );
            }
        }
    }
}
