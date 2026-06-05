//! DSS light/heavy isotope-pair matching for diagnostic-positive scans.

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
    pub charge: u8,
}

#[derive(Debug, Clone, PartialEq)]
pub enum IsotopePairOutcome {
    Paired(IsotopePair),
    Unpaired(ScanRef),
}

/// Match diagnostic-positive scans to DSS light/heavy partners across all files.
pub fn match_isotope_pairs(
    scans: &[ScanRef],
    settings: &Settings,
) -> Vec<IsotopePairOutcome> {
    let mut outcomes = Vec::with_capacity(scans.len());
    let mut used = vec![false; scans.len()];

    for i in 0..scans.len() {
        if used[i] {
            continue;
        }
        let left = &scans[i];
        let charge = left.scan.precursor_charge.unwrap_or(2);

        if let Some(j) = find_partner(i, scans, settings, charge) {
            used[i] = true;
            used[j] = true;
            let right = &scans[j];
            let (light, heavy, light_idx, heavy_idx) = if left.scan.precursor_mz <= right.scan.precursor_mz {
                (left, right, i, j)
            } else {
                (right, left, j, i)
            };
            let _ = (light_idx, heavy_idx);
            outcomes.push(IsotopePairOutcome::Paired(IsotopePair {
                light_file: light.source_file.clone(),
                light_scan: light.scan.scan_number,
                heavy_file: heavy.source_file.clone(),
                heavy_scan: heavy.scan.scan_number,
                rt_light_min: light.scan.retention_time_min,
                rt_heavy_min: heavy.scan.retention_time_min,
                mz_light: light.scan.precursor_mz,
                mz_heavy: heavy.scan.precursor_mz,
                charge,
            }));
        } else {
            used[i] = true;
            outcomes.push(IsotopePairOutcome::Unpaired(left.clone()));
        }
    }

    outcomes
}

fn find_partner(
    index: usize,
    scans: &[ScanRef],
    settings: &Settings,
    charge: u8,
) -> Option<usize> {
    let query = &scans[index];
    let expected_shift = settings.crosslinker_shift_da / f64::from(charge);

    for (j, candidate) in scans.iter().enumerate() {
        if j == index {
            continue;
        }
        let candidate_charge = candidate.scan.precursor_charge.unwrap_or(2);
        if candidate_charge != charge {
            continue;
        }

        let mz_delta = (query.scan.precursor_mz - candidate.scan.precursor_mz).abs();
        if !within_ppm(mz_delta, expected_shift, settings.isotope_pair_ms1_tolerance_ppm) {
            continue;
        }

        let rt_delta = (query.scan.retention_time_min - candidate.scan.retention_time_min).abs();
        if rt_delta > settings.isotope_pair_rt_tolerance_min {
            continue;
        }

        return Some(j);
    }

    None
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
        let outcomes = match_isotope_pairs(&scans, &default_settings());
        assert_eq!(outcomes.len(), 1);
        match &outcomes[0] {
            IsotopePairOutcome::Paired(pair) => {
                assert_eq!(pair.light_scan, 1);
                assert_eq!(pair.heavy_scan, 2);
            }
            IsotopePairOutcome::Unpaired(_) => panic!("expected pair"),
        }
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
        let outcomes = match_isotope_pairs(&scans, &default_settings());
        assert!(outcomes.iter().all(|o| matches!(o, IsotopePairOutcome::Unpaired(_))));
    }
}
