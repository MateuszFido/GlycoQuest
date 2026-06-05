//! Bundled glycan database loading and conversion.

mod catalog;
mod composition;
mod diagnostic;

pub use catalog::{glycan_data_dir, resolve_database, supported_glycan_databases};
pub use diagnostic::DiagnosticIon;
pub use composition::{
    composition_mass, contains_family, load_masses, parse_composition, read_compositions,
    Composition,
};

use catalog::{
    ensure_data_files, residue_targets,
    resolve_database as resolve_catalog_entry,
};
use diagnostic::{expand_diagnostic_ions, load_diagnostic_catalog};

#[derive(Debug, Clone, PartialEq)]
pub struct GlycanEntry {
    pub name: String,
    pub composition: String,
    pub monoisotopic_mass: f64,
    pub diagnostic_ions: Vec<DiagnosticIon>,
    pub residue_targets: Vec<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GlycanLibrary {
    pub database_id: String,
    pub entries: Vec<GlycanEntry>,
}

/// Load a bundled glycan database by catalog id (e.g. `nglyc309`, `oglyc78`).
pub fn load_glycan_database(database_id: &str) -> Result<GlycanLibrary, String> {
    let entry = resolve_catalog_entry(database_id)?;
    ensure_data_files(entry)?;

    let masses = load_masses(&glycan_data_dir().join("glycan_residues.txt"))?;
    let diagnostic_catalog =
        load_diagnostic_catalog(&glycan_data_dir().join("diagnostic_ion_catalog.txt"))?;
    let compositions = read_compositions(&glycan_data_dir().join(entry.filename))?;
    let targets = residue_targets(entry);

    let mut entries = Vec::with_capacity(compositions.len());
    for composition_str in compositions {
        let composition = parse_composition(&composition_str)?;
        let monoisotopic_mass = round_mass(composition_mass(&composition, &masses)?);
        let diagnostic_ions = expand_diagnostic_ions(&composition, &diagnostic_catalog);

        if diagnostic_ions.is_empty() {
            return Err(format!(
                "glycan {composition_str} in database {} has no diagnostic ions",
                entry.id
            ));
        }

        entries.push(GlycanEntry {
            name: composition_str.clone(),
            composition: composition_str,
            monoisotopic_mass,
            diagnostic_ions,
            residue_targets: targets.clone(),
        });
    }

    Ok(GlycanLibrary {
        database_id: entry.id.to_string(),
        entries,
    })
}

fn round_mass(value: f64) -> f64 {
    (value * 1_000_000.0).round() / 1_000_000.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loads_bundled_nglyc309_database() {
        let library = load_glycan_database("nglyc309").expect("bundled nglyc309");
        assert_eq!(library.database_id, "nglyc309");
        assert!(!library.entries.is_empty());
        assert!(library.entries.iter().all(|entry| !entry.diagnostic_ions.is_empty()));
        assert!(library
            .entries
            .iter()
            .any(|entry| entry.composition == "HexNAc(1)"));
    }

    #[test]
    fn loads_bundled_oglyc78_database() {
        let library = load_glycan_database("oglyc78").expect("bundled oglyc78");
        assert_eq!(library.database_id, "oglyc78");
        assert_eq!(library.entries.len(), 78);
        assert!(library
            .entries
            .iter()
            .all(|entry| entry.residue_targets == vec!["S".to_string(), "T".to_string()]));
    }

    #[test]
    fn hexnac_mass_matches_residue_table() {
        let library = load_glycan_database("nglyc309").unwrap();
        let entry = library
            .entries
            .iter()
            .find(|entry| entry.composition == "HexNAc(1)")
            .unwrap();
        assert!((entry.monoisotopic_mass - 203.07937).abs() < 0.0001);
    }

    #[test]
    fn rejects_unknown_database_id() {
        let err = load_glycan_database("not-a-database").unwrap_err();
        assert!(err.contains("Unknown glycan database"));
    }

    #[test]
    fn loads_database_from_env_override_dir() {
        let dir = std::env::temp_dir().join(format!(
            "glycoquest_glycan_data_{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("Nglyc309_Byonic.glyc"), "HexNAc(1)\n").unwrap();
        std::fs::copy(
            glycan_data_dir().join("glycan_residues.txt"),
            dir.join("glycan_residues.txt"),
        )
        .unwrap();
        std::fs::copy(
            glycan_data_dir().join("diagnostic_ion_catalog.txt"),
            dir.join("diagnostic_ion_catalog.txt"),
        )
        .unwrap();

        unsafe { std::env::set_var("GLYCOQUEST_GLYCAN_DATA_DIR", &dir) };
        let library = load_glycan_database("nglyc309").unwrap();
        unsafe { std::env::remove_var("GLYCOQUEST_GLYCAN_DATA_DIR") };
        let _ = std::fs::remove_dir_all(dir);

        assert_eq!(library.entries.len(), 1);
        assert_eq!(library.entries[0].composition, "HexNAc(1)");
    }

    #[test]
    /// Specific cherry-picked glycans from the nglyc309 database.
    fn integration_test_nglyc309() {
        let library = load_glycan_database("nglyc309").unwrap();
        assert_eq!(library.database_id, "nglyc309");
        let total_lines = std::fs::read_to_string(glycan_data_dir().join("Nglyc309_Byonic.glyc"))
            .unwrap()
            .lines()
            .count();
        assert_eq!(total_lines, 309);
        assert_eq!(library.entries.len(), 288); // should be deduplicated

        let entry = library
            .entries
            .iter()
            .enumerate()
            .find(|(i, _)| *i == 61)
            .unwrap()
            .1;
        assert_eq!(entry.composition, "HexNAc(4)Hex(4)Fuc(1)NeuGc(1)");
        assert_eq!(entry.monoisotopic_mass, 1913.677);
        assert!(entry.diagnostic_ions.len() > 50);
        let neu_gc_h2o = entry
            .diagnostic_ions
            .iter()
            .find(|ion| ion.family == "NeuGc" && ion.loss_label == "-H2O")
            .expect("NeuGc -H2O variant");
        assert!((neu_gc_h2o.mz - 290.087006).abs() < 0.001);
        let neu_gc_2h2o = entry
            .diagnostic_ions
            .iter()
            .find(|ion| ion.family == "NeuGc" && ion.loss_label == "-2H2O")
            .expect("NeuGc -2H2O variant");
        assert!((neu_gc_2h2o.mz - 272.076406).abs() < 0.001);
    }
}
