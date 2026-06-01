//! Load internal oxonium rules and derive diagnostic ions and neutral losses.

use std::path::Path;

use super::composition::{self, Composition};
use super::composition::Masses;

pub const PROTON_MASS: f64 = 1.007276;

const CORE_FAMILIES: &[&str] = &["HexNAc", "Hex", "NeuAc", "NeuGc", "Fuc"];

#[derive(Debug, Clone, PartialEq)]
pub struct DiagnosticIon {
    pub family: String,
    pub mz: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct NeutralLoss {
    pub label: String,
    pub delta_da: f64,
}

#[derive(Debug, Clone)]
pub(crate) struct OxoniumRule {
    residue_requirement: String,
    ion_composition: Composition,
    mass_adjustment: f64,
    comment: String,
    diagnostic: bool,
}

pub(crate) fn load_oxonium_rules(path: &Path) -> Result<Vec<OxoniumRule>, String> {
    let content = std::fs::read_to_string(path)
        .map_err(|err| format!("cannot read oxonium ion list {}: {err}", path.display()))?;

    let mut rules = Vec::new();

    for (line_no, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            return Err(format!(
                "malformed oxonium line {} in {}",
                line_no + 1,
                path.display()
            ));
        }

        let residue_requirement = fields[0].trim().to_string();
        let ion_composition = composition::parse_composition(fields[1].trim()).map_err(|err| {
            format!(
                "invalid ion composition on oxonium line {} in {}: {err}",
                line_no + 1,
                path.display()
            )
        })?;
        let mass_adjustment: f64 = fields[2].trim().parse().unwrap_or(0.0);
        let comment = fields[3].trim().to_string();
        let diagnostic = fields
            .get(7)
            .map(|value| value.trim().eq_ignore_ascii_case("true"))
            .unwrap_or(false);

        rules.push(OxoniumRule {
            residue_requirement,
            ion_composition,
            mass_adjustment,
            comment,
            diagnostic,
        });
    }

    if rules.is_empty() {
        return Err(format!("no oxonium rules found in {}", path.display()));
    }

    Ok(rules)
}

pub(crate) fn derive_ions_and_losses(
    glycan: &Composition,
    masses: &Masses,
    rules: &[OxoniumRule],
) -> Result<(Vec<DiagnosticIon>, Vec<NeutralLoss>), String> {
    let mut diagnostics = Vec::new();
    let mut losses = Vec::new();

    for rule in rules {
        if !composition::contains_family(glycan, &rule.residue_requirement) {
            continue;
        }
        if !composition::can_supply(glycan, &rule.ion_composition) {
            continue;
        }

        let ion_mass = match composition::composition_mass(&rule.ion_composition, masses) {
            Ok(mass) => mass,
            Err(_) => continue,
        };
        let adjusted = ion_mass + rule.mass_adjustment;

        if rule.mass_adjustment.abs() > f64::EPSILON {
            let label = neutral_loss_label(&rule.comment, rule.mass_adjustment);
            push_unique_loss(&mut losses, label, rule.mass_adjustment);
        }

        if rule.diagnostic {
            let mz = adjusted + PROTON_MASS;
            push_unique_diagnostic(
                &mut diagnostics,
                rule.residue_requirement.clone(),
                mz,
            );
        }
    }

    for family in CORE_FAMILIES {
        if !composition::contains_family(glycan, family) {
            continue;
        }
        let residue = composition::canonical_residue(family);
        let mass = masses.get(residue).ok_or_else(|| {
            format!("missing core family residue mass for {family}")
        })?;
        let mz = mass + PROTON_MASS;
        push_unique_diagnostic(&mut diagnostics, family.to_string(), mz);
    }

    Ok((diagnostics, losses))
}

fn neutral_loss_label(comment: &str, adjustment: f64) -> String {
    if !comment.is_empty() {
        comment.to_string()
    } else {
        format!("loss@{:.6}", adjustment.abs())
    }
}

fn push_unique_diagnostic(out: &mut Vec<DiagnosticIon>, family: String, mz: f64) {
    let rounded = round_mz(mz);
    if out
        .iter()
        .any(|ion| ion.family == family && (ion.mz - rounded).abs() < 1e-4)
    {
        return;
    }
    out.push(DiagnosticIon { family, mz: rounded });
}

fn push_unique_loss(out: &mut Vec<NeutralLoss>, label: String, delta_da: f64) {
    let delta = round_mz(delta_da);
    if out
        .iter()
        .any(|loss| loss.label == label && (loss.delta_da - delta).abs() < 1e-4)
    {
        return;
    }
    out.push(NeutralLoss { label, delta_da: delta });
}

fn round_mz(value: f64) -> f64 {
    (value * 1_000_000.0).round() / 1_000_000.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn test_masses() -> Masses {
        HashMap::from([
            ("HexNAc".into(), 203.07937),
            ("Hex".into(), 162.05282),
            ("Fuc".into(), 146.05791),
            ("dHex".into(), 146.05791),
            ("NeuAc".into(), 291.09542),
            ("NeuGc".into(), 307.09033),
        ])
    }

    #[test]
    fn hexnac_single_gets_canonical_diagnostic() {
        let glycan = composition::parse_composition("HexNAc(1)").unwrap();
        let rules = vec![OxoniumRule {
            residue_requirement: "HexNAc".into(),
            ion_composition: composition::parse_composition("HexNAc(1)").unwrap(),
            mass_adjustment: 0.0,
            comment: String::new(),
            diagnostic: false,
        }];

        let (diagnostics, _) =
            derive_ions_and_losses(&glycan, &test_masses(), &rules).unwrap();
        assert!(diagnostics.iter().any(|ion| ion.family == "HexNAc"));
        let hexnac = diagnostics
            .iter()
            .find(|ion| ion.family == "HexNAc")
            .unwrap();
        assert!((hexnac.mz - 204.086646).abs() < 0.001);
    }

    #[test]
    fn neuac_glycan_gets_oxonium_diagnostic() {
        let glycan = composition::parse_composition("HexNAc(1)NeuAc(1)").unwrap();
        let rules = vec![OxoniumRule {
            residue_requirement: "NeuAc".into(),
            ion_composition: composition::parse_composition("NeuAc(1)").unwrap(),
            mass_adjustment: 0.0,
            comment: String::new(),
            diagnostic: true,
        }];

        let (diagnostics, _) =
            derive_ions_and_losses(&glycan, &test_masses(), &rules).unwrap();
        assert!(diagnostics.iter().any(|ion| ion.family == "NeuAc"));
    }
}
