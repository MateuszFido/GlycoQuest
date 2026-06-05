//! Spectrum-local glycan candidate pruning by matched diagnostic families.

use crate::glyco::{contains_family, parse_composition, Composition, GlycanLibrary};

#[derive(Debug, Clone, PartialEq)]
pub struct PrunedGlycan {
    pub name: String,
    pub composition: String,
}

/// Keep glycans whose composition contains every matched diagnostic family.
pub fn prune_glycans(
    matched_families: &[String],
    library: &GlycanLibrary,
) -> Result<Vec<PrunedGlycan>, String> {
    if matched_families.is_empty() {
        return Ok(Vec::new());
    }

    let mut candidates = Vec::new();
    for entry in &library.entries {
        let composition = parse_composition(&entry.composition)?;
        if families_supported(&composition, matched_families) {
            candidates.push(PrunedGlycan {
                name: entry.name.clone(),
                composition: entry.composition.clone(),
            });
        }
    }
    Ok(candidates)
}

fn families_supported(glycan: &Composition, families: &[String]) -> bool {
    families
        .iter()
        .all(|family| contains_family(glycan, family))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::glyco::{GlycanEntry, GlycanLibrary};

    fn test_library() -> GlycanLibrary {
        GlycanLibrary {
            database_id: "test".into(),
            entries: vec![
                GlycanEntry {
                    name: "HexNAc(1)".into(),
                    composition: "HexNAc(1)".into(),
                    monoisotopic_mass: 203.079373,
                    diagnostic_ions: vec![],
                    residue_targets: vec!["N".into()],
                },
                GlycanEntry {
                    name: "HexNAc(1)NeuAc(1)".into(),
                    composition: "HexNAc(1)NeuAc(1)".into(),
                    monoisotopic_mass: 494.174789,
                    diagnostic_ions: vec![],
                    residue_targets: vec!["N".into()],
                },
                GlycanEntry {
                    name: "NeuAc(1)".into(),
                    composition: "NeuAc(1)".into(),
                    monoisotopic_mass: 291.095417,
                    diagnostic_ions: vec![],
                    residue_targets: vec!["N".into()],
                },
            ],
        }
    }

    #[test]
    fn hexnac_only_keeps_hexnac_glycans() {
        let library = test_library();
        let pruned = prune_glycans(&["HexNAc".into()], &library).unwrap();
        let names: Vec<_> = pruned.iter().map(|g| g.name.as_str()).collect();
        assert!(names.contains(&"HexNAc(1)"));
        assert!(names.contains(&"HexNAc(1)NeuAc(1)"));
        assert!(!names.contains(&"NeuAc(1)"));
    }

    #[test]
    fn hexnac_and_neuac_keeps_intersection() {
        let library = test_library();
        let pruned = prune_glycans(&["HexNAc".into(), "NeuAc".into()], &library).unwrap();
        assert_eq!(pruned.len(), 1);
        assert_eq!(pruned[0].name, "HexNAc(1)NeuAc(1)");
    }
}
