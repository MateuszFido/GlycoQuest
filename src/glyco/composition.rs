//! Parse Byonic-style glycan composition strings such as `HexNAc(2)Hex(5)Fuc(1)`.

use std::collections::BTreeMap;
use std::collections::HashMap;
use std::path::Path;

/// Parsed residue counts, keyed by the canonical residue name.
pub type Composition = BTreeMap<String, u32>;
/// Mapping of name → monoisotopic mass (Da).
pub type Masses = HashMap<String, f64>;

/// Map glycan residue aliases to canonical class names (mass-equivalent abbreviations).
pub(crate) fn canonical_residue(name: &str) -> &str {
    match name {
        "Fuc" | "F" | "fucose" | "dhex" => "dHex",
        "Man" | "Glc" | "Gal" | "H" | "glc" | "gal" | "hexose" | "glucose" | "galactose"
        | "mannose" => "Hex",
        "GlcA" | "GalA" | "IdoA" | "hexuronic" | "hexuronate" => "HexA",
        "GlcNAc" | "GalNAc" | "N" => "HexNAc",
        "Neu5Ac" | "A" => "NeuAc",
        "Neu5Gc" | "G" => "NeuGc",
        "X" | "xylose" | "xyl" | "rib" | "ribose" | "ara" | "arabinose" | "lyx" | "pentose" => {
            "Pent"
        }
        "Phospho" => "PO3H",
        "Sulfo" => "SO3H",
        other => other,
    }
}

/// Load internal glycan residue masses from a file.
pub fn load_masses(path: &Path) -> Result<Masses, String> {
    let content = std::fs::read_to_string(path)
        .map_err(|err| format!("cannot read glycan residue file {}: {err}", path.display()))?;

    let mut masses = Masses::new();

    for (line_no, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let mut fields = line.split('\t');
        let name = fields
            .next()
            .ok_or_else(|| {
                format!(
                    "malformed residue line {} in {}",
                    line_no + 1,
                    path.display()
                )
            })?
            .trim();
        let mass_str = fields.next().ok_or_else(|| {
            format!(
                "missing mass on residue line {} in {}",
                line_no + 1,
                path.display()
            )
        })?;
        let mass: f64 = mass_str.trim().parse().map_err(|_| {
            format!(
                "invalid mass on residue line {} in {}: {mass_str}",
                line_no + 1,
                path.display()
            )
        })?;
        if mass <= 0.0 {
            return Err(format!(
                "non-positive mass for residue {name} on line {} in {}",
                line_no + 1,
                path.display()
            ));
        }

        masses.insert(name.to_string(), mass);
    }

    if masses.is_empty() {
        return Err(format!("no residues found in {}", path.display()));
    }

    normalize_masses(masses)
}

fn normalize_masses(raw: Masses) -> Result<Masses, String> {
    let mut normalized = Masses::new();
    for (name, mass) in raw {
        let canonical = canonical_residue(&name).to_string();
        if let Some(&existing) = normalized.get(&canonical) {
            if (existing - mass).abs() > 1e-6 {
                return Err(format!(
                    "conflicting masses for canonical residue {canonical}: {existing} vs {mass} (from {name})"
                ));
            }
        } else {
            normalized.insert(canonical, mass);
        }
    }
    Ok(normalized)
}

pub fn composition_mass(composition: &Composition, masses: &Masses) -> Result<f64, String> {
    let mut total = 0.0;
    for (residue, count) in composition {
        let mass = masses
            .get(residue)
            .ok_or_else(|| format!("unknown glycan residue in composition: {residue}"))?;
        total += mass * f64::from(*count);
    }
    Ok(total)
}

/// One composition string per non-empty line; duplicates keep the first occurrence.
pub fn read_compositions(path: &Path) -> Result<Vec<String>, String> {
    let content = std::fs::read_to_string(path)
        .map_err(|err| format!("cannot read glycan library {}: {err}", path.display()))?;

    let mut compositions = Vec::new();
    let mut seen = std::collections::HashSet::new();
    let mut duplicate_count = 0usize;

    for (line_no, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        if !seen.insert(line.to_string()) {
            duplicate_count += 1;
            eprintln!(
                "warning: duplicate glycan composition on line {} of {}: {line}",
                line_no + 1,
                path.display()
            );
            continue;
        }

        compositions.push(line.to_string());
    }

    if compositions.is_empty() {
        return Err(format!(
            "no glycan compositions found in {}",
            path.display()
        ));
    }

    if duplicate_count > 0 {
        eprintln!(
            "warning: skipped {duplicate_count} duplicate composition(s) in {}",
            path.display()
        );
    }

    Ok(compositions)
}

/// Parse a composition string into residue counts.
pub fn parse_composition(line: &str) -> Result<Composition, String> {
    let line = line.trim();
    if line.is_empty() {
        return Err("empty glycan composition".into());
    }

    let mut composition = Composition::new();
    let mut rest = line;

    while !rest.is_empty() {
        let open = rest
            .find('(')
            .ok_or_else(|| format!("invalid glycan composition (missing '('): {line}"))?;
        let close = rest[open..]
            .find(')')
            .ok_or_else(|| format!("invalid glycan composition (missing ')'): {line}"))?;
        let close = open + close;

        let residue = rest[..open].trim();
        if residue.is_empty() {
            return Err(format!(
                "invalid glycan composition (empty residue): {line}"
            ));
        }

        let count_str = rest[open + 1..close].trim();
        let count: u32 = count_str
            .parse()
            .map_err(|_| format!("invalid residue count in composition: {line}"))?;
        if count == 0 {
            return Err(format!(
                "residue count must be positive in composition: {line}"
            ));
        }

        let canonical = canonical_residue(residue).to_string();
        *composition.entry(canonical).or_insert(0) += count;

        rest = rest[close + 1..].trim_start();
    }

    if composition.is_empty() {
        return Err(format!("invalid glycan composition: {line}"));
    }

    Ok(composition)
}

/// True when `glycan` contains at least one of the required family residue.
pub fn contains_family(glycan: &Composition, family: &str) -> bool {
    glycan.get(canonical_residue(family)).copied().unwrap_or(0) > 0
}

/// True when `glycan` has enough counts to supply `needed`.
pub fn can_supply(glycan: &Composition, needed: &Composition) -> bool {
    needed
        .iter()
        .all(|(residue, needed_count)| glycan.get(residue).copied().unwrap_or(0) >= *needed_count)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn skips_hash_header_and_loads_masses() {
        let path = std::env::temp_dir().join(format!(
            "glycoquest_residue_test_{}.txt",
            std::process::id()
        ));
        let mut file = std::fs::File::create(&path).unwrap();
        writeln!(file, "#Name\tMass").unwrap();
        writeln!(file, "HexNAc\t203.07937").unwrap();
        writeln!(file, "Fuc\t146.05791").unwrap();

        let masses = load_masses(&path).unwrap();
        assert_eq!(masses.get("HexNAc"), Some(&203.07937));
        assert_eq!(masses.get("dHex"), Some(&146.05791));
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn parses_single_residue() {
        let comp = parse_composition("HexNAc(1)").unwrap();
        assert_eq!(comp.get("HexNAc"), Some(&1));
    }

    #[test]
    fn parses_multiple_residues() {
        let comp = parse_composition("HexNAc(2)Hex(5)Fuc(1)").unwrap();
        assert_eq!(comp.get("HexNAc"), Some(&2));
        assert_eq!(comp.get("Hex"), Some(&5));
        assert_eq!(comp.get("dHex"), Some(&1));
    }

    #[test]
    fn fuc_and_dhex_are_interchangeable() {
        let from_fuc = parse_composition("HexNAc(1)Fuc(1)").unwrap();
        let from_dhex = parse_composition("HexNAc(1)dHex(1)").unwrap();
        assert_eq!(from_fuc, from_dhex);
    }

    #[test]
    fn hex_aliases_merge() {
        let comp = parse_composition("Man(2)Gal(1)").unwrap();
        assert_eq!(comp.get("Hex"), Some(&3));
    }

    #[test]
    fn modification_aliases_canonicalize() {
        let phospho = parse_composition("Hex(1)Phospho(1)").unwrap();
        assert_eq!(phospho.get("PO3H"), Some(&1));
        let sulfo = parse_composition("HexNAc(1)Sulfo(1)").unwrap();
        assert_eq!(sulfo.get("SO3H"), Some(&1));
    }

    #[test]
    fn rejects_conflicting_masses_for_same_canonical_class() {
        let path = std::env::temp_dir().join(format!(
            "glycoquest_residue_conflict_{}.txt",
            std::process::id()
        ));
        let mut file = std::fs::File::create(&path).unwrap();
        writeln!(file, "Fuc\t146.05791").unwrap();
        writeln!(file, "dHex\t999.0").unwrap();

        let err = load_masses(&path).unwrap_err();
        assert!(err.contains("conflicting masses"));
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn rejects_zero_count() {
        assert!(parse_composition("HexNAc(0)").is_err());
    }

    #[test]
    fn dedupes_duplicate_compositions() {
        let path =
            std::env::temp_dir().join(format!("glycoquest_glyc_test_{}.glyc", std::process::id()));
        let mut file = std::fs::File::create(&path).unwrap();
        writeln!(file, "HexNAc(1)").unwrap();
        writeln!(file, "HexNAc(2)").unwrap();
        writeln!(file, "HexNAc(1)").unwrap();

        let compositions = read_compositions(&path).unwrap();
        assert_eq!(
            compositions,
            vec!["HexNAc(1)".to_string(), "HexNAc(2)".to_string()]
        );
        let _ = std::fs::remove_file(path);
    }
}
