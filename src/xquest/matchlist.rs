// Copyright (c) ETH Zurich, Mateusz Fido

//! xQuest compare_peaks matchlist generation.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::crosslinker::CrosslinkerProfile;
use crate::jobs::{PlannedJob, SpectrumKey};
use crate::prefilter::{FilteredSpectrum, PrefilterResult};

#[derive(Debug, Clone, PartialEq)]
pub enum MatchlistRow {
    Paired {
        id: String,
        precursor_mz: f64,
        charge: u8,
        light_label: String,
        heavy_label: String,
        light_scan: u32,
        heavy_scan: u32,
        rt_light_sec: f64,
        rt_heavy_sec: f64,
        mz_light: f64,
        mz_heavy: f64,
    },
    Single {
        id: String,
        precursor_mz: f64,
        charge: u8,
        label: String,
        scan: u32,
        rt_sec: f64,
    },
}

/// Constant-time lookup shared by every job generated from one prefilter run.
pub struct FilteredSpectrumIndex<'a> {
    by_key: HashMap<SpectrumKey, &'a FilteredSpectrum>,
}

impl<'a> FilteredSpectrumIndex<'a> {
    pub fn new(prefilter: &'a PrefilterResult) -> Self {
        let mut by_key = HashMap::with_capacity(prefilter.filtered.len());
        for spectrum in &prefilter.filtered {
            by_key
                .entry(SpectrumKey {
                    source_file: spectrum.source_file.clone(),
                    scan_number: spectrum.scan_number,
                })
                .or_insert(spectrum);
        }
        Self { by_key }
    }

    fn get(&self, key: &SpectrumKey) -> Option<&'a FilteredSpectrum> {
        self.by_key.get(key).copied()
    }
}

/// Build a per-spectrum label for matchlist columns 4/5. compare_peaks3.pl uses
/// `basename(col4)_basename(col5)` both as its redundancy key and as the spectrum
/// name that flows into the result XML. The scan is placed first so the label is
/// unique per pair (fixing xQuest's over-aggressive dedup) and the leading integer
/// stays trivially parseable regardless of the source filename, while the real
/// mzXML stem is preserved for readable output.
fn spectrum_label(scan: u32, stem: &str) -> String {
    format!("{scan}.{stem}")
}

pub fn build_matchlist(
    job: &PlannedJob,
    prefilter: &PrefilterResult,
    spectrum_index: &FilteredSpectrumIndex<'_>,
    stem: &str,
    crosslinker: &CrosslinkerProfile,
) -> Result<Vec<MatchlistRow>, String> {
    let mut rows = Vec::new();

    if crosslinker.requires_isotope_pair_prefilter() {
        let job_scans: HashSet<&SpectrumKey> = job.spectrum_keys.iter().collect();
        let mut seen_pairs = HashSet::new();

        for pair in &prefilter.isotope_pairs {
            let light_key = SpectrumKey {
                source_file: pair.light_file.clone(),
                scan_number: pair.light_scan,
            };
            let heavy_key = SpectrumKey {
                source_file: pair.heavy_file.clone(),
                scan_number: pair.heavy_scan,
            };

            if !job_scans.contains(&light_key) || !job_scans.contains(&heavy_key) {
                continue;
            }

            if !seen_pairs.insert((pair.light_scan, pair.heavy_scan)) {
                continue;
            }

            let light_spec = spectrum_index
                .get(&light_key)
                .ok_or_else(|| missing_spec(&light_key))?;

            rows.push(MatchlistRow::Paired {
                id: format!("{},{}", pair.light_scan, pair.heavy_scan),
                precursor_mz: light_spec.precursor_mz,
                charge: pair.light_charge,
                light_label: spectrum_label(pair.light_scan, stem),
                heavy_label: spectrum_label(pair.heavy_scan, stem),
                light_scan: pair.light_scan,
                heavy_scan: pair.heavy_scan,
                rt_light_sec: pair.rt_light_min * 60.0,
                rt_heavy_sec: pair.rt_heavy_min * 60.0,
                mz_light: pair.mz_light,
                mz_heavy: pair.mz_heavy,
            });
        }
    } else {
        for key in job.spectrum_keys.iter() {
            let spec = spectrum_index.get(key).ok_or_else(|| missing_spec(key))?;
            rows.push(MatchlistRow::Single {
                id: spec.scan_number.to_string(),
                precursor_mz: spec.precursor_mz,
                charge: spec.precursor_charge.unwrap_or(2),
                label: spectrum_label(spec.scan_number, stem),
                scan: spec.scan_number,
                rt_sec: spec.retention_time_min * 60.0,
            });
        }
    }

    Ok(rows)
}

fn missing_spec(key: &SpectrumKey) -> String {
    format!(
        "filtered spectrum missing for {} scan {}",
        key.source_file.display(),
        key.scan_number
    )
}

/// Write compare_peaks3 `-match` input (tab-separated rows).
pub fn write_matchlist(path: &Path, rows: &[MatchlistRow]) -> Result<(), String> {
    let file = File::create(path).map_err(|err| err.to_string())?;
    let mut w = BufWriter::new(file);

    for row in rows {
        match row {
            MatchlistRow::Paired {
                id,
                precursor_mz,
                charge,
                light_label,
                heavy_label,
                light_scan,
                heavy_scan,
                rt_light_sec,
                rt_heavy_sec,
                mz_light,
                mz_heavy,
            } => {
                writeln!(
                    w,
                    "{id}\t{precursor_mz:.6}\t{charge}\t{light_label}\t{heavy_label}\tlight\theavy\t{light_scan}:{heavy_scan}\t{rt_light_sec:.4}:{rt_heavy_sec:.4}\t{mz_light:.8}:{mz_heavy:.8}"
                )
                .map_err(|err| err.to_string())?;
            }
            MatchlistRow::Single {
                id,
                precursor_mz,
                charge,
                label,
                scan,
                rt_sec,
            } => {
                writeln!(
                    w,
                    "{id}\t{precursor_mz:.6}\t{charge}\t{label}\t{label}\tlight\tlight\t{scan}:{scan}\t{rt_sec:.4}:{rt_sec:.4}\t{precursor_mz:.8}:{precursor_mz:.8}"
                )
                .map_err(|err| err.to_string())?;
            }
        }
    }

    Ok(())
}

pub fn isotopepairs_path(matchlist_path: &Path) -> PathBuf {
    let stem = matchlist_path
        .file_stem()
        .map(|stem| stem.to_string_lossy().into_owned())
        .unwrap_or_else(|| matchlist_path.display().to_string());
    let parent = matchlist_path.parent().unwrap_or_else(|| Path::new("."));
    parent.join(format!("{stem}_isotopepairs.txt"))
}

pub fn specxml_filename(result_dir: &str) -> String {
    format!("{result_dir}.spec.xml")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_compare_peaks_matchlist_row() {
        let rows = vec![MatchlistRow::Paired {
            id: "1,2".into(),
            precursor_mz: 500.0,
            charge: 2,
            light_label: spectrum_label(1, "sample.c"),
            heavy_label: spectrum_label(2, "sample.c"),
            light_scan: 1,
            heavy_scan: 2,
            rt_light_sec: 1200.0,
            rt_heavy_sec: 1206.0,
            mz_light: 500.0,
            mz_heavy: 506.0376605,
        }];
        let path =
            std::env::temp_dir().join(format!("glycoquest_matchlist_{}", std::process::id()));
        write_matchlist(&path, &rows).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        // Columns 4/5 must be distinct per pair so compare_peaks3.pl does not
        // treat every spectral pair as redundant.
        assert!(content.contains("1,2\t500.000000\t2\t1.sample.c\t2.sample.c\t"));
        assert!(content.contains("\tlight\theavy\t1:2\t"));
        let _ = std::fs::remove_file(path);
    }
}
