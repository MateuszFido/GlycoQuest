//! Load diagnostic ion catalog and expand search targets (m/z × neutral loss delta).

use std::path::Path;

use super::composition::{self, Composition};

#[derive(Debug, Clone, PartialEq)]
pub struct DiagnosticIon {
    pub family: String,
    pub mz: f64,
    pub loss_label: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DiagnosticTemplate {
    pub family: String,
    pub mz: f64,
    pub label: String,
    pub composition: Option<Composition>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct NeutralLoss {
    pub label: String,
    pub delta_da: f64,
}

#[derive(Debug, Clone)]
pub struct DiagnosticCatalog {
    pub templates: Vec<DiagnosticTemplate>,
    pub losses: Vec<NeutralLoss>,
}

pub fn load_diagnostic_catalog(path: &Path) -> Result<DiagnosticCatalog, String> {
    let content = std::fs::read_to_string(path).map_err(|err| {
        format!(
            "cannot read diagnostic ion catalog {}: {err}",
            path.display()
        )
    })?;

    let mut templates = Vec::new();
    let mut losses = Vec::new();
    let mut in_losses = false;

    for (line_no, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') {
            if line.contains("neutral_loss") {
                in_losses = true;
            }
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();

        if in_losses {
            if fields.len() < 1 {
                return Err(format!(
                    "malformed neutral loss line {} in {}",
                    line_no + 1,
                    path.display()
                ));
            }
            let delta_da: f64 = fields[0].trim().parse().map_err(|_| {
                format!(
                    "invalid delta on neutral loss line {} in {}",
                    line_no + 1,
                    path.display()
                )
            })?;
            let label = fields.get(1).map(|s| s.trim().to_string()).unwrap_or_default();
            push_unique_loss(&mut losses, label, delta_da);
            continue;
        }

        if fields.len() < 2 {
            return Err(format!(
                "malformed diagnostic template line {} in {}",
                line_no + 1,
                path.display()
            ));
        }

        let family = composition::canonical_residue(fields[0].trim()).to_string();
        let mz: f64 = fields[1].trim().parse().map_err(|_| {
            format!(
                "invalid m/z on line {} in {}",
                line_no + 1,
                path.display()
            )
        })?;
        let label = fields.get(2).map(|s| s.trim().to_string()).unwrap_or_default();
        let composition = match fields.get(3) {
            Some(comp) if !comp.trim().is_empty() => {
                Some(composition::parse_composition(comp.trim()).map_err(|err| {
                    format!(
                        "invalid composition on line {} in {}: {err}",
                        line_no + 1,
                        path.display()
                    )
                })?)
            }
            _ => None,
        };

        templates.push(DiagnosticTemplate {
            family,
            mz: round_mz(mz),
            label,
            composition,
        });
    }

    if templates.is_empty() {
        return Err(format!(
            "no diagnostic templates found in {}",
            path.display()
        ));
    }
    if losses.is_empty() {
        return Err(format!(
            "no neutral losses found in {}",
            path.display()
        ));
    }

    Ok(DiagnosticCatalog {
        templates,
        losses,
    })
}

pub fn expand_diagnostic_ions(
    glycan: &Composition,
    catalog: &DiagnosticCatalog,
) -> Vec<DiagnosticIon> {
    let mut ions = Vec::new();

    for template in &catalog.templates {
        if !template_eligible(glycan, template) {
            continue;
        }
        for loss in &catalog.losses {
            let mz = round_mz(template.mz + loss.delta_da);
            push_unique_ion(&mut ions, template.family.clone(), mz, loss.label.clone());
        }
    }

    ions.sort_by(|a, b| {
        a.family
            .cmp(&b.family)
            .then(a.mz.partial_cmp(&b.mz).unwrap_or(std::cmp::Ordering::Equal))
    });
    ions
}

fn template_eligible(glycan: &Composition, template: &DiagnosticTemplate) -> bool {
    if !composition::contains_family(glycan, &template.family) {
        return false;
    }
    match &template.composition {
        Some(needed) => composition::can_supply(glycan, needed),
        None => true,
    }
}

fn push_unique_ion(out: &mut Vec<DiagnosticIon>, family: String, mz: f64, loss_label: String) {
    if out.iter().any(|ion| {
        ion.family == family
            && ion.loss_label == loss_label
            && (ion.mz - mz).abs() < 1e-4
    }) {
        return;
    }
    out.push(DiagnosticIon {
        family,
        mz,
        loss_label,
    });
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

    fn fixture_catalog() -> DiagnosticCatalog {
        DiagnosticCatalog {
            templates: vec![DiagnosticTemplate {
                family: "HexNAc".into(),
                mz: 204.086646,
                label: "HexNAc".into(),
                composition: Some(composition::parse_composition("HexNAc(1)").unwrap()),
            }],
            losses: vec![
                NeutralLoss {
                    label: String::new(),
                    delta_da: 0.0,
                },
                NeutralLoss {
                    label: "-H2O".into(),
                    delta_da: -18.0106,
                },
                NeutralLoss {
                    label: "-2H2O".into(),
                    delta_da: -36.0212,
                },
            ],
        }
    }

    #[test]
    fn hexnac_single_expands_to_three_mz() {
        let glycan = composition::parse_composition("HexNAc(1)").unwrap();
        let ions = expand_diagnostic_ions(&glycan, &fixture_catalog());
        assert_eq!(ions.len(), 3);
        assert!(ions.iter().all(|ion| ion.family == "HexNAc"));
        assert!(ions.iter().any(|ion| (ion.mz - 204.086646).abs() < 1e-4));
        assert!(ions.iter().any(|ion| (ion.mz - 186.076046).abs() < 1e-3));
        assert!(ions.iter().any(|ion| (ion.mz - 168.065446).abs() < 1e-3));
    }

    #[test]
    fn neuac_glycan_gets_neuac_template() {
        let glycan = composition::parse_composition("HexNAc(1)NeuAc(1)").unwrap();
        let catalog = DiagnosticCatalog {
            templates: vec![DiagnosticTemplate {
                family: "NeuAc".into(),
                mz: 292.102696,
                label: "NeuAc".into(),
                composition: Some(composition::parse_composition("NeuAc(1)").unwrap()),
            }],
            losses: vec![NeutralLoss {
                label: String::new(),
                delta_da: 0.0,
            }],
        };
        let ions = expand_diagnostic_ions(&glycan, &catalog);
        assert!(ions.iter().any(|ion| ion.family == "NeuAc"));
    }

    #[test]
    fn loads_bundled_catalog() {
        let path = crate::glycan::glycan_data_dir().join("diagnostic_ion_catalog.txt");
        let catalog = load_diagnostic_catalog(&path).expect("bundled catalog");
        assert!(!catalog.templates.is_empty());
        assert_eq!(catalog.losses.len(), 3);
    }
}
