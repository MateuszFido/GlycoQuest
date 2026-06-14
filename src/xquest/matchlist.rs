//! xQuest compare_peaks matchlist generation.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::crosslinker::CrosslinkerProfile;
use crate::jobs::{filtered_for_key, isotope_pair_for_scan, PlannedJob, SpectrumKey};
use crate::prefilter::PrefilterResult;

#[derive(Debug, Clone, PartialEq)]
pub enum MatchlistRow {
    Paired {
        id: String,
        precursor_mz: f64,
        charge: u8,
        mzxml_file: String,
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
        mzxml_file: String,
        scan: u32,
        rt_sec: f64,
    },
}

pub fn build_matchlist(
    job: &PlannedJob,
    prefilter: &PrefilterResult,
    pruned_mzxml: &Path,
    crosslinker: &CrosslinkerProfile,
) -> Result<Vec<MatchlistRow>, String> {
    let mzxml_file = pruned_mzxml.display().to_string();
    let mut rows = Vec::new();

    if crosslinker.requires_isotope_pair_prefilter() {
        for key in &job.spectrum_keys {
            let spec = filtered_for_key(prefilter, key).ok_or_else(|| missing_spec(key))?;
            let pair = isotope_pair_for_scan(
                &prefilter.isotope_pairs,
                &key.source_file,
                key.scan_number,
            )
            .ok_or_else(|| {
                format!(
                    "missing isotope pair for {} scan {}",
                    key.source_file.display(),
                    key.scan_number
                )
            })?;

            let (light_scan, heavy_scan, rt_light_sec, rt_heavy_sec, mz_light, mz_heavy) =
                if pair.light_file == key.source_file && pair.light_scan == key.scan_number {
                    (
                        pair.light_scan,
                        pair.heavy_scan,
                        pair.rt_light_min * 60.0,
                        pair.rt_heavy_min * 60.0,
                        pair.mz_light,
                        pair.mz_heavy,
                    )
                } else {
                    (
                        pair.light_scan,
                        pair.heavy_scan,
                        pair.rt_light_min * 60.0,
                        pair.rt_heavy_min * 60.0,
                        pair.mz_light,
                        pair.mz_heavy,
                    )
                };

            if rows.iter().any(|row| matchlist_light_scan(row) == Some(light_scan)) {
                continue;
            }

            rows.push(MatchlistRow::Paired {
                id: format!("{light_scan},{heavy_scan}"),
                precursor_mz: spec.precursor_mz,
                charge: spec.precursor_charge.unwrap_or(2),
                mzxml_file: mzxml_file.clone(),
                light_scan,
                heavy_scan,
                rt_light_sec,
                rt_heavy_sec,
                mz_light,
                mz_heavy,
            });
        }
    } else {
        for key in &job.spectrum_keys {
            let spec = filtered_for_key(prefilter, key).ok_or_else(|| missing_spec(key))?;
            rows.push(MatchlistRow::Single {
                id: spec.scan_number.to_string(),
                precursor_mz: spec.precursor_mz,
                charge: spec.precursor_charge.unwrap_or(2),
                mzxml_file: mzxml_file.clone(),
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

fn matchlist_light_scan(row: &MatchlistRow) -> Option<u32> {
    match row {
        MatchlistRow::Paired { light_scan, .. } => Some(*light_scan),
        MatchlistRow::Single { scan, .. } => Some(*scan),
    }
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
                mzxml_file,
                light_scan,
                heavy_scan,
                rt_light_sec,
                rt_heavy_sec,
                mz_light,
                mz_heavy,
            } => {
                writeln!(
                    w,
                    "{id}\t{precursor_mz:.6}\t{charge}\t{mzxml_file}\t{mzxml_file}\tlight\theavy\t{light_scan}:{heavy_scan}\t{rt_light_sec:.4}:{rt_heavy_sec:.4}\t{mz_light:.8}:{mz_heavy:.8}"
                )
                .map_err(|err| err.to_string())?;
            }
            MatchlistRow::Single {
                id,
                precursor_mz,
                charge,
                mzxml_file,
                scan,
                rt_sec,
            } => {
                writeln!(
                    w,
                    "{id}\t{precursor_mz:.6}\t{charge}\t{mzxml_file}\t{mzxml_file}\tlight\tlight\t{scan}:{scan}\t{rt_sec:.4}:{rt_sec:.4}\t{precursor_mz:.8}:{precursor_mz:.8}"
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
            mzxml_file: "./input.mzXML".into(),
            light_scan: 1,
            heavy_scan: 2,
            rt_light_sec: 1200.0,
            rt_heavy_sec: 1206.0,
            mz_light: 500.0,
            mz_heavy: 506.0376605,
        }];
        let path = std::env::temp_dir().join(format!(
            "glycoquest_matchlist_{}",
            std::process::id()
        ));
        write_matchlist(&path, &rows).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("1,2\t500.000000\t2\t./input.mzXML"));
        assert!(content.contains("\tlight\theavy\t1:2\t"));
        let _ = std::fs::remove_file(path);
    }
}
