//! Diagnostic-ion matching against expanded glycan library targets.

use std::collections::BTreeSet;

use crate::glyco::{DiagnosticIon, GlycanLibrary};
use crate::mzxml::Ms2Scan;

#[derive(Debug, Clone, PartialEq)]
pub struct MatchedIon {
    pub family: String,
    pub expected_mz: f64,
    pub observed_mz: f64,
    pub loss_label: String,
    pub peak_index: usize,
    pub intensity: f64,
    pub error_ppm: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DiagnosticMatch {
    pub matched_ions: Vec<MatchedIon>,
    pub matched_families: Vec<String>,
    pub passes: bool,
}

pub fn match_diagnostic_ions(
    scan: &Ms2Scan,
    library: &GlycanLibrary,
    tolerance_ppm: f64,
) -> DiagnosticMatch {
    let targets = unique_diagnostic_targets(library);
    let mut matched_ions = Vec::new();
    let mut families = BTreeSet::new();

    for (peak_index, (obs_mz, intensity)) in scan.peaks.iter().enumerate() {
        for ion in &targets {
            if within_ppm(*obs_mz, ion.mz, tolerance_ppm) {
                families.insert(ion.family.clone());
                matched_ions.push(MatchedIon {
                    family: ion.family.clone(),
                    expected_mz: ion.mz,
                    observed_mz: *obs_mz,
                    loss_label: ion.loss_label.clone(),
                    peak_index,
                    intensity: *intensity,
                    error_ppm: ppm_error(*obs_mz, ion.mz),
                });
            }
        }
    }

    matched_ions.sort_by(|a, b| {
        a.family.cmp(&b.family).then(
            a.expected_mz
                .partial_cmp(&b.expected_mz)
                .unwrap_or(std::cmp::Ordering::Equal),
        )
    });

    let matched_families: Vec<String> = families.into_iter().collect();
    let passes = !matched_families.is_empty();

    DiagnosticMatch {
        matched_ions,
        matched_families,
        passes,
    }
}

fn unique_diagnostic_targets(library: &GlycanLibrary) -> Vec<DiagnosticIon> {
    let mut ions = Vec::new();
    for entry in &library.entries {
        for ion in &entry.diagnostic_ions {
            if ions.iter().any(|existing: &DiagnosticIon| {
                existing.family == ion.family
                    && existing.loss_label == ion.loss_label
                    && (existing.mz - ion.mz).abs() < 1e-4
            }) {
                continue;
            }
            ions.push(ion.clone());
        }
    }
    ions
}

fn within_ppm(observed: f64, expected: f64, tolerance_ppm: f64) -> bool {
    if expected <= 0.0 {
        return false;
    }
    ppm_error(observed, expected).abs() <= tolerance_ppm
}

fn ppm_error(observed: f64, expected: f64) -> f64 {
    if expected <= 0.0 {
        return 0.0;
    }
    ((observed - expected) / expected) * 1_000_000.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::glyco::GlycanEntry;
    use crate::mzxml::Ms2Scan;

    fn hexnac_library() -> GlycanLibrary {
        GlycanLibrary {
            database_id: "test".into(),
            entries: vec![GlycanEntry {
                name: "HexNAc(1)".into(),
                composition: "HexNAc(1)".into(),
                monoisotopic_mass: 203.079373,
                diagnostic_ions: vec![
                    DiagnosticIon {
                        family: "HexNAc".into(),
                        mz: 204.086646,
                        loss_label: String::new(),
                    },
                    DiagnosticIon {
                        family: "HexNAc".into(),
                        mz: 186.076046,
                        loss_label: "-H2O".into(),
                    },
                ],
                residue_targets: vec!["N".into()],
            }],
        }
    }

    #[test]
    fn hexnac_positive_scan_passes() {
        let scan = Ms2Scan {
            scan_number: 1,
            retention_time_min: 10.0,
            precursor_mz: 800.0,
            precursor_charge: Some(2),
            peaks: vec![(100.0, 1000.0), (204.0867, 5000.0)],
        };
        let result = match_diagnostic_ions(&scan, &hexnac_library(), 10.0);
        assert!(result.passes);
        assert!(result.matched_families.contains(&"HexNAc".to_string()));
        assert_eq!(result.matched_ions[0].peak_index, 1);
        assert_eq!(result.matched_ions[0].intensity, 5000.0);
        assert!((result.matched_ions[0].error_ppm - 0.264595).abs() < 0.001);
    }

    #[test]
    fn no_diagnostic_scan_fails() {
        let scan = Ms2Scan {
            scan_number: 1,
            retention_time_min: 10.0,
            precursor_mz: 800.0,
            precursor_charge: Some(2),
            peaks: vec![(100.0, 1000.0)],
        };
        let result = match_diagnostic_ions(&scan, &hexnac_library(), 10.0);
        assert!(!result.passes);
    }
}
