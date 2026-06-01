//! Bundled glycan database catalog and data directory resolution.

use std::path::PathBuf;

const GLYCAN_DATA_ENV: &str = "GLYCOQUEST_GLYCAN_DATA_DIR";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GlycanType {
    N,
    O,
}

#[derive(Debug, Clone, Copy)]
pub struct DatabaseEntry {
    pub id: &'static str,
    pub filename: &'static str,
    pub glycan_type: GlycanType,
    pub aliases: &'static [&'static str],
}

const DATABASES: &[DatabaseEntry] = &[
    DatabaseEntry {
        id: "nglyc309",
        filename: "Nglyc309_Byonic.glyc",
        glycan_type: GlycanType::N,
        aliases: &["n-glycan", "nglyc"],
    },
    DatabaseEntry {
        id: "oglyc78",
        filename: "Oglyc78_Byonic.glyc",
        glycan_type: GlycanType::O,
        aliases: &["o-glycan", "oglyc"],
    },
];

pub fn supported_glycan_databases() -> &'static [(&'static str, &'static str)] {
    &[
        ("nglyc309", "N-linked glycans (309 compositions)"),
        ("oglyc78", "O-linked glycans (78 compositions)"),
    ]
}

pub fn supported_ids_message() -> String {
    supported_glycan_databases()
        .iter()
        .map(|(id, _)| *id)
        .collect::<Vec<_>>()
        .join(", ")
}

fn normalize_id(input: &str) -> String {
    input.trim().to_ascii_lowercase()
}

pub fn resolve_database(input: &str) -> Result<&'static DatabaseEntry, String> {
    let trimmed = input.trim();
    if trimmed.is_empty() {
        return Err("glycan database id must not be empty".into());
    }

    let normalized = normalize_id(trimmed);
    DATABASES
        .iter()
        .find(|entry| {
            entry.id == normalized.as_str()
                || entry
                    .aliases
                    .iter()
                    .any(|alias| normalize_id(alias) == normalized)
        })
        .ok_or_else(|| {
            format!(
                "Unknown glycan database: {trimmed}; supported: {}",
                supported_ids_message()
            )
        })
}

pub fn glycan_data_dir() -> PathBuf {
    if let Ok(dir) = std::env::var(GLYCAN_DATA_ENV) {
        if !dir.trim().is_empty() {
            return PathBuf::from(dir);
        }
    }

    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("databases");
    if manifest.is_dir() {
        return manifest;
    }

    let cwd = PathBuf::from("databases");
    if cwd.is_dir() {
        return cwd;
    }

    manifest
}

pub fn residue_targets(entry: &DatabaseEntry) -> Vec<String> {
    match entry.glycan_type {
        GlycanType::N => vec!["N".into()],
        GlycanType::O => vec!["S".into(), "T".into()],
    }
}

pub fn ensure_data_files(entry: &DatabaseEntry) -> Result<(), String> {
    for path in [
        glycan_data_dir().join(entry.filename),
        glycan_data_dir().join("glycan_residues.txt"),
        glycan_data_dir().join("oxonium_ion_list.txt")
    ] {
        if !path.is_file() {
            return Err(format!(
                "bundled glycan data file not found: {} (set {GLYCAN_DATA_ENV} if running outside the repo)",
                path.display()
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolves_catalog_ids_and_aliases() {
        assert_eq!(resolve_database("nglyc309").unwrap().id, "nglyc309");
        assert_eq!(resolve_database("N-glycan").unwrap().id, "nglyc309");
        assert_eq!(resolve_database("oglyc78").unwrap().id, "oglyc78");
    }

    #[test]
    fn rejects_unknown_database() {
        assert!(resolve_database("unknown").unwrap_err().contains("Unknown glycan database"));
    }
}
