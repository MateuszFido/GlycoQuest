//! Write reduced mzXML files containing only retained scans.

use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use base64::Engine;
use base64::engine::general_purpose::STANDARD;

use crate::mzxml::{Ms2Scan, parse_scans};
use crate::prefilter::FilteredSpectrum;

/// Write `spectra/` containing one reduced mzXML per source file (prefilter-retained scans only).
pub fn write_prefiltered_mzxml(
    spectra_dir: &Path,
    filtered: &[FilteredSpectrum],
) -> Result<Vec<PathBuf>, String> {
    let target_dir = spectra_dir;
    fs::create_dir_all(target_dir).map_err(|err| {
        format!(
            "cannot create spectra directory {}: {err}",
            target_dir.display()
        )
    })?;

    let mut by_file: std::collections::HashMap<PathBuf, HashSet<u32>> =
        std::collections::HashMap::new();
    for row in filtered {
        by_file
            .entry(row.source_file.clone())
            .or_default()
            .insert(row.scan_number);
    }

    let mut written = Vec::new();
    for (source, scans) in by_file {
        let out_path = target_dir.join(
            source
                .file_name()
                .ok_or_else(|| format!("invalid source path: {}", source.display()))?,
        );
        write_subset_mzxml(&source, &out_path, &scans)?;
        written.push(out_path);
    }

    Ok(written)
}

fn write_subset_mzxml(source: &Path, dest: &Path, keep: &HashSet<u32>) -> Result<(), String> {
    let scans: Vec<Ms2Scan> = parse_scans(source)?
        .into_iter()
        .filter(|scan| keep.contains(&scan.scan_number))
        .collect();

    if scans.is_empty() {
        return Err(format!(
            "no matching scans to write for {}",
            source.display()
        ));
    }

    let mut w = BufWriter::new(
        File::create(dest).map_err(|err| format!("cannot write {}: {err}", dest.display()))?,
    );
    writeln!(
        w,
        r#"<?xml version="1.0" encoding="UTF-8"?>
<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.1" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://sashimi.sourceforge.net/schema_revision/mzXML_3.1 http://sashimi.sourceforge.net/schema_revision/mzXML_3.1/mzXML_3.1.xsd" version="3.1">
<msRun>"#
    )
    .map_err(|err| err.to_string())?;

    for scan in scans {
        write_scan(&mut w, &scan)?;
    }

    writeln!(w, "</msRun>\n</mzXML>").map_err(|err| err.to_string())?;
    Ok(())
}

fn write_scan(w: &mut impl Write, scan: &Ms2Scan) -> Result<(), String> {
    let rt_seconds = scan.retention_time_min * 60.0;
    let peaks_b64 = encode_peaks_xquest_base64(&scan.peaks);
    let charge_attr = scan
        .precursor_charge
        .map(|charge| format!(r#" precursorCharge="{charge}""#))
        .unwrap_or_default();

    writeln!(
        w,
        r#"<scan num="{num}" msLevel="2" retentionTime="PT{rt_seconds:.6}S" peaksCount="{}">
<precursorMz{charge_attr}>{precursor:.6}</precursorMz>
<peaks precision="32" byteOrder="network" contentType="m/z-int">{peaks_b64}</peaks>
</scan>"#,
        scan.peaks.len(),
        num = scan.scan_number,
        rt_seconds = rt_seconds,
        charge_attr = charge_attr,
        precursor = scan.precursor_mz,
        peaks_b64 = peaks_b64,
    )
    .map_err(|err| err.to_string())
}

fn encode_peaks_xquest_base64(peaks: &[(f64, f64)]) -> String {
    let mut raw = Vec::with_capacity(peaks.len() * 8);
    for (mz, intensity) in peaks {
        raw.extend_from_slice(&(*mz as f32).to_be_bytes());
        raw.extend_from_slice(&(*intensity as f32).to_be_bytes());
    }
    STANDARD.encode(raw)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prefilter::FilteredSpectrum;
    use std::path::PathBuf;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/mzxml")
            .join(name)
    }

    #[test]
    fn writes_subset_for_dss_pair() {
        let source = fixture("dss_pair.mzXML");
        let out = std::env::temp_dir().join(format!("glycoquest_pruned_{}", std::process::id()));
        let _ = fs::remove_dir_all(&out);
        fs::create_dir_all(&out).unwrap();

        let filtered = vec![
            FilteredSpectrum {
                source_file: source.clone(),
                scan_number: 1,
                retention_time_min: 20.0,
                precursor_mz: 500.0,
                precursor_charge: Some(2),
                matched_families: vec!["HexNAc".into()],
                matched_ions: vec![],
            },
            FilteredSpectrum {
                source_file: source.clone(),
                scan_number: 2,
                retention_time_min: 20.1,
                precursor_mz: 506.0376605,
                precursor_charge: Some(2),
                matched_families: vec!["HexNAc".into()],
                matched_ions: vec![],
            },
        ];

        let paths = write_prefiltered_mzxml(&out.join("spectra"), &filtered).unwrap();
        assert_eq!(paths.len(), 1);
        let text = fs::read_to_string(&paths[0]).unwrap();
        assert!(text.contains("precision=\"32\""));
        assert!(text.matches("<scan").count() == 2);
        let _ = fs::remove_dir_all(out);
    }
}
