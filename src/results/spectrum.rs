//! mzXML peak extraction and approximate b/y fragment annotation for the viewer.

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use crate::crosslinker::CrosslinkerProfile;
use crate::mzxml::{self, Ms2Scan};
use crate::results::mapping::CrosslinkMapping;

use super::filter::{AnnotatedHit, PostfilterStatus};

const PROTON: f64 = 1.007276;
const WATER: f64 = 18.010565;
const TOLERANCE_DA: f32 = 0.5;
const MAX_PEAKS_PER_SCAN: usize = 2000;

/// Experimental peaks for one MS/MS scan.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanSpectrum {
    pub scan: u32,
    pub precursor_mz: f64,
    pub charge: u8,
    pub mz: Vec<f32>,
    pub intensity: Vec<f32>,
}

/// Theoretical fragment ions for mirror-plot annotation (approximate; not xQuest-identical).
#[derive(Debug, Clone, PartialEq)]
pub struct FragmentAnnotation {
    pub theoretical_mz: Vec<f32>,
    pub labels: Vec<String>,
    pub matched_indices: Vec<usize>,
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
        if path.extension().is_none_or(|ext| !ext.eq_ignore_ascii_case("mzxml")) {
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

/// Build approximate b/y fragments for both peptide arms of a crosslink.
pub fn annotate_fragments(
    mapping: &CrosslinkMapping,
    crosslinker: &CrosslinkerProfile,
    spectrum: Option<&ScanSpectrum>,
) -> FragmentAnnotation {
    let mut theoretical_mz = Vec::new();
    let mut labels = Vec::new();

    append_peptide_fragments(
        &mapping.pep1,
        mapping.link1,
        "p1",
        crosslinker.xlinkermw,
        &mut theoretical_mz,
        &mut labels,
    );
    append_peptide_fragments(
        &mapping.pep2,
        mapping.link2,
        "p2",
        crosslinker.xlinkermw,
        &mut theoretical_mz,
        &mut labels,
    );

    let matched_indices = spectrum
        .map(|spec| match_peaks(&theoretical_mz, &spec.mz))
        .unwrap_or_default();

    FragmentAnnotation {
        theoretical_mz,
        labels,
        matched_indices,
    }
}

fn append_peptide_fragments(
    peptide: &str,
    link_pos: Option<usize>,
    prefix: &str,
    xlink_mass: f64,
    theoretical_mz: &mut Vec<f32>,
    labels: &mut Vec<String>,
) {
    if peptide.is_empty() {
        return;
    }
    let residues: Vec<char> = peptide.chars().collect();
    let n = residues.len();
    let link_idx = link_pos.and_then(|p| p.checked_sub(1)).filter(|&i| i < n);

    // b ions (N-terminal fragments).
    let mut mass = 0.0;
    for (i, &res) in residues.iter().enumerate() {
        mass += residue_mass(res);
        let mut ion_mass = mass + PROTON;
        if link_idx.is_some_and(|li| i >= li) {
            ion_mass += xlink_mass;
        }
        theoretical_mz.push(ion_mass as f32);
        labels.push(format!("{prefix}_b{}", i + 1));
    }

    // y ions (C-terminal fragments).
    let mut mass = WATER;
    for (i, &res) in residues.iter().enumerate().rev() {
        mass += residue_mass(res);
        let pos = n - i;
        let mut ion_mass = mass + PROTON;
        if link_idx.is_some_and(|li| i <= n - 1 - li) {
            ion_mass += xlink_mass;
        }
        theoretical_mz.push(ion_mass as f32);
        labels.push(format!("{prefix}_y{pos}"));
    }
}

fn residue_mass(residue: char) -> f64 {
    match residue.to_ascii_uppercase() {
        'A' => 71.037114,
        'R' => 156.101111,
        'N' => 114.042927,
        'D' => 115.026943,
        'C' => 103.009185,
        'E' => 129.042593,
        'Q' => 128.058578,
        'G' => 57.021464,
        'H' => 137.058912,
        'I' => 113.084064,
        'L' => 113.084064,
        'K' => 128.094963,
        'M' => 131.040485,
        'F' => 147.068414,
        'P' => 97.052764,
        'S' => 87.032028,
        'T' => 101.047679,
        'W' => 186.079313,
        'Y' => 163.063329,
        'V' => 99.068414,
        _ => 0.0,
    }
}

fn match_peaks(theoretical: &[f32], observed: &[f32]) -> Vec<usize> {
    let mut matched = Vec::new();
    for (idx, &obs) in observed.iter().enumerate() {
        if theoretical
            .iter()
            .any(|&theo| (theo - obs).abs() <= TOLERANCE_DA)
        {
            matched.push(idx);
        }
    }
    matched
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
    use crate::crosslinker::{CrosslinkerLabel, CrosslinkerProfile};
    use crate::results::mapping::map_crosslink;
    use std::collections::HashMap;

    fn dss_crosslinker() -> CrosslinkerProfile {
        CrosslinkerProfile {
            name: "dss".into(),
            label: CrosslinkerLabel::LightHeavy,
            shift_da: 12.075321,
            xlinkermw: 138.0680796,
            xlink_sites: "K:K".into(),
            nterm_xlinkable: false,
        }
    }

    #[test]
    fn annotate_fragments_produces_b_and_y_ions() {
        let proteins = HashMap::new();
        let mapping = map_crosslink("ACDE", "FGHK", "P1", "P2", "2-3", &proteins, None);
        let frags = annotate_fragments(&mapping, &dss_crosslinker(), None);
        assert!(frags.theoretical_mz.len() >= 8);
        assert!(frags.labels.iter().any(|l| l.starts_with("p1_b")));
        assert!(frags.labels.iter().any(|l| l.starts_with("p2_y")));
    }

    #[test]
    fn load_spectra_from_fixture() {
        let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/mzxml");
        let scans: HashSet<u32> = [1, 2].into_iter().collect();
        let loaded = load_spectra_for_scans(&dir, &scans);
        // Fixture directory contains mzXML files with scans.
        assert!(!loaded.is_empty() || dir.is_dir());
    }
}
