// Copyright (c) ETH Zurich, Mateusz Fido

//! DSS light/heavy isotope-pair matching for diagnostic-positive scans.

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use crate::cli::settings::Settings;
use crate::mzxml::Ms2Scan;

#[derive(Debug, Clone, PartialEq)]
pub struct ScanRef {
    pub source_file: PathBuf,
    pub scan: Ms2Scan,
}

#[derive(Debug, Clone, PartialEq)]
pub struct IsotopePair {
    pub light_file: PathBuf,
    pub light_scan: u32,
    pub heavy_file: PathBuf,
    pub heavy_scan: u32,
    pub rt_light_min: f64,
    pub rt_heavy_min: f64,
    pub mz_light: f64,
    pub mz_heavy: f64,
    pub light_charge: u8,
    pub heavy_charge: u8,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct IsotopeMatchResult {
    pub pairs: Vec<IsotopePair>,
    pub unpaired: Vec<ScanRef>,
}

/// Match diagnostic-positive scans to DSS light/heavy partners within each file.
///
/// Every valid pair is retained (non-greedy). Precursor charges may differ; pairing
/// uses monoisotopic precursor mass and the configured crosslinker shift in Da.
pub fn match_isotope_pairs(scans: &[ScanRef], settings: &Settings) -> IsotopeMatchResult {
    let mut indices_by_file: HashMap<&PathBuf, Vec<usize>> = HashMap::new();
    for (index, scan) in scans.iter().enumerate() {
        indices_by_file
            .entry(&scan.source_file)
            .or_default()
            .push(index);
    }

    let mut pairs = Vec::new();
    let mut paired_indices = HashSet::new();

    for indices in indices_by_file.values() {
        for left_pos in 0..indices.len() {
            for right_pos in left_pos + 1..indices.len() {
                let left = &scans[indices[left_pos]];
                let right = &scans[indices[right_pos]];
                if forms_pair(left, right, settings) {
                    paired_indices.insert(indices[left_pos]);
                    paired_indices.insert(indices[right_pos]);
                    pairs.push(build_pair(left, right));
                }
            }
        }
    }

    let unpaired = scans
        .iter()
        .enumerate()
        .filter(|(index, _)| !paired_indices.contains(index))
        .map(|(_, scan)| scan.clone())
        .collect();

    IsotopeMatchResult { pairs, unpaired }
}

fn forms_pair(left: &ScanRef, right: &ScanRef, settings: &Settings) -> bool {
    if left.source_file != right.source_file {
        return false;
    }

    let mass_left = precursor_mass(&left.scan);
    let mass_right = precursor_mass(&right.scan);
    let mass_delta = (mass_left - mass_right).abs();
    if !within_ppm(
        mass_delta,
        settings.crosslinker_shift_da,
        settings.isotope_pair_ms1_tolerance_ppm,
    ) {
        return false;
    }

    let rt_delta = (left.scan.retention_time_min - right.scan.retention_time_min).abs();
    rt_delta <= settings.isotope_pair_rt_tolerance_min
}

fn build_pair(left: &ScanRef, right: &ScanRef) -> IsotopePair {
    let mass_left = precursor_mass(&left.scan);
    let mass_right = precursor_mass(&right.scan);
    let (light, heavy) = if mass_left <= mass_right {
        (left, right)
    } else {
        (right, left)
    };

    IsotopePair {
        light_file: light.source_file.clone(),
        light_scan: light.scan.scan_number,
        heavy_file: heavy.source_file.clone(),
        heavy_scan: heavy.scan.scan_number,
        rt_light_min: light.scan.retention_time_min,
        rt_heavy_min: heavy.scan.retention_time_min,
        mz_light: light.scan.precursor_mz,
        mz_heavy: heavy.scan.precursor_mz,
        light_charge: light.scan.precursor_charge.unwrap_or(2),
        heavy_charge: heavy.scan.precursor_charge.unwrap_or(2),
    }
}

/// Proton monoisotopic mass (Da), used to convert m/z to neutral mass.
const PROTON_MASS: f64 = 1.007276;

fn precursor_mass(scan: &Ms2Scan) -> f64 {
    let charge = f64::from(scan.precursor_charge.unwrap_or(2));
    // Neutral monoisotopic mass so the isotope shift compares correctly even
    // when a light/heavy pair is observed at different precursor charges.
    scan.precursor_mz * charge - charge * PROTON_MASS
}

fn within_ppm(observed: f64, expected: f64, tolerance_ppm: f64) -> bool {
    if expected <= 0.0 {
        return false;
    }
    ((observed - expected).abs() / expected) * 1_000_000.0 <= tolerance_ppm
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mzxml::Ms2Scan;

    fn scan(num: u32, rt: f64, mz: f64, charge: u8) -> Ms2Scan {
        Ms2Scan {
            scan_number: num,
            retention_time_min: rt,
            precursor_mz: mz,
            precursor_charge: Some(charge),
            peaks: vec![(204.0867, 1000.0)],
        }
    }

    fn default_settings() -> Settings {
        Settings::defaults()
    }

    #[test]
    fn finds_dss_pair_at_expected_shift() {
        let file = PathBuf::from("test.mzXML");
        let scans = vec![
            ScanRef {
                source_file: file.clone(),
                scan: scan(1, 20.0, 500.0, 2),
            },
            ScanRef {
                source_file: file.clone(),
                scan: scan(2, 20.1, 506.0376605, 2),
            },
        ];
        let result = match_isotope_pairs(&scans, &default_settings());
        assert_eq!(result.pairs.len(), 1);
        assert_eq!(result.unpaired.len(), 0);
        assert_eq!(result.pairs[0].light_scan, 1);
        assert_eq!(result.pairs[0].heavy_scan, 2);
    }

    #[test]
    fn unpaired_when_shift_wrong() {
        let file = PathBuf::from("test.mzXML");
        let scans = vec![
            ScanRef {
                source_file: file.clone(),
                scan: scan(1, 20.0, 500.0, 2),
            },
            ScanRef {
                source_file: file.clone(),
                scan: scan(2, 20.1, 510.0, 2),
            },
        ];
        let result = match_isotope_pairs(&scans, &default_settings());
        assert!(result.pairs.is_empty());
        assert_eq!(result.unpaired.len(), 2);
    }

    #[test]
    fn pairs_within_file_only() {
        let file_a = PathBuf::from("a.mzXML");
        let file_b = PathBuf::from("b.mzXML");
        let scans = vec![
            ScanRef {
                source_file: file_a.clone(),
                scan: scan(1, 20.0, 500.0, 2),
            },
            ScanRef {
                source_file: file_b,
                scan: scan(2, 20.1, 506.0376605, 2),
            },
        ];
        let result = match_isotope_pairs(&scans, &default_settings());
        assert!(result.pairs.is_empty());
        assert_eq!(result.unpaired.len(), 2);
    }

    fn mz_for(neutral_mass: f64, charge: u8) -> f64 {
        neutral_mass / f64::from(charge) + super::PROTON_MASS
    }

    #[test]
    fn pairs_across_different_precursor_charges_by_mass() {
        let file = PathBuf::from("test.mzXML");
        let light_mass = 2000.0;
        let heavy_mass = light_mass + 12.075321;
        let scans = vec![
            ScanRef {
                source_file: file.clone(),
                scan: scan(1, 20.0, mz_for(light_mass, 4), 4),
            },
            ScanRef {
                source_file: file.clone(),
                scan: scan(2, 20.1, mz_for(heavy_mass, 3), 3),
            },
        ];
        let result = match_isotope_pairs(&scans, &default_settings());
        assert_eq!(result.pairs.len(), 1);
        assert_eq!(result.pairs[0].light_charge, 4);
        assert_eq!(result.pairs[0].heavy_charge, 3);
    }

    #[test]
    fn retains_all_valid_pairs_non_greedy() {
        let file = PathBuf::from("test.mzXML");
        let base_mass = 2000.0;
        let shift = 12.075321;
        let scans = vec![
            ScanRef {
                source_file: file.clone(),
                scan: scan(1, 20.0, base_mass / 2.0, 2),
            },
            ScanRef {
                source_file: file.clone(),
                scan: scan(2, 20.05, (base_mass + shift) / 2.0, 2),
            },
            ScanRef {
                source_file: file.clone(),
                scan: scan(3, 20.1, (base_mass + shift) / 2.0, 2),
            },
        ];
        let result = match_isotope_pairs(&scans, &default_settings());
        assert_eq!(result.pairs.len(), 2);
        assert_eq!(result.unpaired.len(), 0);
    }
}
