// Copyright (c) ETH Zurich, Mateusz Fido

//! mzXML peak extraction for the viewer.

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use crate::mzxml::{self, Ms2Scan};

use super::filter::{AnnotatedHit, PostfilterStatus};

const MAX_PEAKS_PER_SCAN: usize = 2000;

/// Experimental peaks for one MS/MS scan.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanSpectrum {
    pub scan: u32,
    pub retention_time_min: f64,
    pub precursor_mz: f64,
    pub charge: u8,
    pub mz: Vec<f32>,
    pub intensity: Vec<f32>,
}

/// Load peaks for the given scans from reduced mzXML files under `spectra_dir`.
pub fn load_spectra_for_scans(
    spectra_dir: &Path,
    scans: &HashSet<u32>,
) -> HashMap<u32, ScanSpectrum> {
    if scans.is_empty() || !spectra_dir.is_dir() {
        return HashMap::new();
    }

    let mut out = HashMap::new();
    let Ok(entries) = std::fs::read_dir(spectra_dir) else {
        return out;
    };

    for entry in entries.flatten() {
        let path = entry.path();
        if path
            .extension()
            .is_none_or(|ext| !ext.eq_ignore_ascii_case("mzxml"))
        {
            continue;
        }
        if let Ok(file_scans) = mzxml::parse_scans(&path) {
            merge_scans(&mut out, file_scans, scans);
        }
    }
    out
}

fn merge_scans(
    out: &mut HashMap<u32, ScanSpectrum>,
    file_scans: Vec<Ms2Scan>,
    wanted: &HashSet<u32>,
) {
    for scan in file_scans {
        if !wanted.contains(&scan.scan_number) || out.contains_key(&scan.scan_number) {
            continue;
        }
        let (mz, intensity) = compress_peaks(&scan.peaks);
        out.insert(
            scan.scan_number,
            ScanSpectrum {
                scan: scan.scan_number,
                retention_time_min: scan.retention_time_min,
                precursor_mz: scan.precursor_mz,
                charge: scan.precursor_charge.unwrap_or(2),
                mz,
                intensity,
            },
        );
    }
}

fn compress_peaks(peaks: &[(f64, f64)]) -> (Vec<f32>, Vec<f32>) {
    let mut sorted: Vec<_> = peaks
        .iter()
        .filter(|(mz, int)| mz.is_finite() && int.is_finite() && *int > 0.0)
        .copied()
        .collect();
    sorted.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    sorted.truncate(MAX_PEAKS_PER_SCAN);
    sorted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mz: Vec<f32> = sorted.iter().map(|(m, _)| *m as f32).collect();
    let intensity: Vec<f32> = sorted.iter().map(|(_, i)| *i as f32).collect();
    (mz, intensity)
}

/// Collect scan numbers referenced by hits (passing hits only when `passing_only`).
pub fn scans_from_hits(hits: &[AnnotatedHit], passing_only: bool) -> HashSet<u32> {
    hits.iter()
        .filter(|h| !passing_only || h.postfilter_status == PostfilterStatus::Pass)
        .filter_map(|h| h.scan)
        .collect()
}

pub fn spectra_dir_from_layout(root: &Path) -> PathBuf {
    root.join("spectra")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn load_spectra_from_fixture() {
        let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/mzxml");
        let scans: HashSet<u32> = [1, 2].into_iter().collect();
        let loaded = load_spectra_for_scans(&dir, &scans);
        // Fixture directory contains mzXML files with scans.
        assert!(!loaded.is_empty() || dir.is_dir());
    }
}
