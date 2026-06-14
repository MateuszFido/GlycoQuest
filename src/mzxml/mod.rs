//! Narrow mzXML reader for MS2 scans used by GlycoQuest prefilters.

mod write;

use std::io::Read;
use std::path::Path;

pub use write::write_prefiltered_mzxml;

#[derive(Debug, Clone, PartialEq)]
pub struct Ms2Scan {
    pub scan_number: u32,
    pub retention_time_min: f64,
    pub precursor_mz: f64,
    pub precursor_charge: Option<u8>,
    pub peaks: Vec<(f64, f64)>,
}

/// Parse all MS2 scans from an mzXML file.
pub fn parse_scans(path: &Path) -> Result<Vec<Ms2Scan>, String> {
    let content = std::fs::read_to_string(path).map_err(|err| {
        format!("cannot read mzXML file {}: {err}", path.display())
    })?;

    if !content.contains("<scan") {
        return Err(format!("no scans found in mzXML file: {}", path.display()));
    }

    let mut scans = Vec::new();
    let mut pos = 0usize;

    while let Some(start) = content[pos..].find("<scan") {
        let abs_start = pos + start;
        let tag_end = content[abs_start..]
            .find('>')
            .ok_or_else(|| format!("malformed scan tag in {}", path.display()))?;
        let scan_open_end = abs_start + tag_end + 1;

        let close_tag = "</scan>";
        let close_rel = content[scan_open_end..]
            .find(close_tag)
            .ok_or_else(|| format!("unclosed scan in {}", path.display()))?;
        let scan_body = &content[scan_open_end..scan_open_end + close_rel];
        let scan_open_tag = &content[abs_start..scan_open_end];

        if let Some(scan) = parse_scan_block(scan_open_tag, scan_body, path)? {
            scans.push(scan);
        }

        pos = scan_open_end + close_rel + close_tag.len();
    }

    Ok(scans)
}

fn parse_scan_block(
    scan_open_tag: &str,
    scan_body: &str,
    path: &Path,
) -> Result<Option<Ms2Scan>, String> {
    let ms_level = extract_tag_value(scan_body, "msLevel")
        .or_else(|| extract_attribute(scan_open_tag, "msLevel"))
        .unwrap_or_else(|| "2".to_string());
    if ms_level != "2" {
        return Ok(None);
    }

    let scan_number = extract_attribute(scan_open_tag, "num")
        .or_else(|| extract_tag_value(scan_body, "num"))
        .and_then(|v| v.parse().ok())
        .ok_or_else(|| format!("scan missing num attribute in {}", path.display()))?;

    let retention_time_min = parse_retention_time(scan_body, scan_open_tag)?;

    let (precursor_mz, precursor_charge) = parse_precursor(scan_body)?;

    let peaks = parse_peaks(scan_body, path, scan_number)?;

    Ok(Some(Ms2Scan {
        scan_number,
        retention_time_min,
        precursor_mz,
        precursor_charge,
        peaks,
    }))
}

fn parse_precursor(scan_body: &str) -> Result<(f64, Option<u8>), String> {
    if let Some(start) = scan_body.find("<precursorMz") {
        let tag_end = scan_body[start..]
            .find('>')
            .ok_or_else(|| "malformed precursorMz tag".to_string())?
            + start;
        let open_tag = &scan_body[start..=tag_end];
        let value_start = tag_end + 1;
        let value_end = scan_body[value_start..]
            .find("</precursorMz>")
            .ok_or_else(|| "unclosed precursorMz tag".to_string())?
            + value_start;
        let value = scan_body[value_start..value_end].trim();
        let precursor_mz = value
            .split_whitespace()
            .next()
            .and_then(|v| v.parse().ok())
            .ok_or_else(|| format!("invalid precursorMz value: {value}"))?;
        let charge = extract_attribute(open_tag, "precursorCharge")
            .or_else(|| extract_tag_value(scan_body, "precursorCharge"))
            .and_then(|v| v.parse().ok());
        return Ok((precursor_mz, charge));
    }

    let precursor_mz = extract_tag_value(scan_body, "precursorMz")
        .or_else(|| extract_tag_value(scan_body, "precursorMZ"))
        .and_then(|v| v.split_whitespace().next()?.parse().ok())
        .ok_or_else(|| "scan missing precursorMz".to_string())?;
    let precursor_charge = extract_tag_value(scan_body, "precursorCharge")
        .or_else(|| extract_tag_value(scan_body, "precursor charge"))
        .and_then(|v| v.parse().ok());
    Ok((precursor_mz, precursor_charge))
}

fn parse_retention_time(scan_body: &str, scan_open_tag: &str) -> Result<f64, String> {
    if let Some(rt) = extract_attribute(scan_open_tag, "retTime") {
        if let Ok(minutes) = rt.parse::<f64>() {
            return Ok(minutes);
        }
    }

    if let Some(rt) = extract_attribute(scan_open_tag, "retentionTime") {
        return parse_retention_time_value(&rt);
    }

    if let Some(rt) = extract_tag_value(scan_body, "retentionTime") {
        return parse_retention_time_value(&rt);
    }

    Err("scan missing retention time".into())
}

fn parse_retention_time_value(raw: &str) -> Result<f64, String> {
    let trimmed = raw.trim();
    if let Some(seconds) = trimmed.strip_prefix("PT").and_then(|s| s.strip_suffix('S')) {
        let sec: f64 = seconds.parse().map_err(|_| format!("invalid retention time: {raw}"))?;
        return Ok(sec / 60.0);
    }
    trimmed
        .parse::<f64>()
        .map_err(|_| format!("invalid retention time: {raw}"))
}

fn parse_peaks(scan_body: &str, path: &Path, scan_number: u32) -> Result<Vec<(f64, f64)>, String> {
    let start = scan_body
        .find("<peaks")
        .ok_or_else(|| format!("scan {scan_number} missing peaks in {}", path.display()))?;
    let tag_end = scan_body[start..]
        .find('>')
        .ok_or_else(|| format!("malformed peaks tag in {}", path.display()))?;
    let peaks_tag = &scan_body[start..=start + tag_end];
    let content_start = start + tag_end + 1;
    let content_end = scan_body[content_start..]
        .find("</peaks>")
        .ok_or_else(|| format!("unclosed peaks in {}", path.display()))?;
    let peaks_text = scan_body[content_start..content_start + content_end].trim();

    if peaks_text.is_empty() {
        return Ok(Vec::new());
    }

    let compression = extract_attribute(peaks_tag, "compressionType");
    let precision = extract_attribute(peaks_tag, "precision");

    if compression.as_deref() == Some("zlib") {
        return decode_compressed_peaks(peaks_text, peaks_tag, scan_number, path);
    }

    if precision.as_deref().is_some_and(|p| p == "32" || p == "64") {
        return decode_base64_peaks(peaks_text, peaks_tag, scan_number, path);
    }

    parse_plain_peaks(peaks_text, scan_number, path)
}

fn parse_plain_peaks(
    peaks_text: &str,
    scan_number: u32,
    path: &Path,
) -> Result<Vec<(f64, f64)>, String> {
    let values: Result<Vec<f64>, String> = peaks_text
        .split_whitespace()
        .map(|v| {
            v.parse::<f64>()
                .map_err(|_| format!("invalid peak value in scan {scan_number}: {v}"))
        })
        .collect();
    let values = values?;

    if values.len() % 2 != 0 {
        return Err(format!(
            "odd number of peak values in scan {scan_number} in {}",
            path.display()
        ));
    }

    let mut peaks = Vec::with_capacity(values.len() / 2);
    for chunk in values.chunks(2) {
        peaks.push((chunk[0], chunk[1]));
    }
    Ok(peaks)
}

fn decode_base64_bytes(peaks_text: &str, scan_number: u32, path: &Path) -> Result<Vec<u8>, String> {
    use base64::Engine;
    use base64::engine::general_purpose::STANDARD;

    STANDARD
        .decode(peaks_text.replace('\n', "").replace('\r', ""))
        .map_err(|err| {
            format!(
                "invalid base64 peaks in scan {scan_number} in {}: {err}",
                path.display()
            )
        })
}

fn decode_binary_peak_values(
    raw: &[u8],
    peaks_tag: &str,
    scan_number: u32,
    path: &Path,
) -> Result<Vec<f64>, String> {
    let precision = extract_attribute(peaks_tag, "precision").unwrap_or_else(|| "32".into());
    let byte_order = extract_attribute(peaks_tag, "byteOrder").unwrap_or_else(|| "network".into());
    let little_endian = byte_order.eq_ignore_ascii_case("little");

    match precision.as_str() {
        "64" => decode_f64_peak_values(raw, little_endian, scan_number, path),
        "32" => decode_f32_peak_values(raw, little_endian, scan_number, path),
        other => Err(format!(
            "unsupported peaks precision {other} in scan {scan_number} in {}",
            path.display()
        )),
    }
}

fn values_to_peak_pairs(values: Vec<f64>, scan_number: u32, path: &Path) -> Result<Vec<(f64, f64)>, String> {
    if values.len() % 2 != 0 {
        return Err(format!(
            "odd number of peak values in scan {scan_number} in {}",
            path.display()
        ));
    }

    let mut peaks = Vec::with_capacity(values.len() / 2);
    for chunk in values.chunks(2) {
        peaks.push((chunk[0], chunk[1]));
    }
    Ok(peaks)
}

fn decode_base64_peaks(
    peaks_text: &str,
    peaks_tag: &str,
    scan_number: u32,
    path: &Path,
) -> Result<Vec<(f64, f64)>, String> {
    let raw = decode_base64_bytes(peaks_text, scan_number, path)?;
    let values = decode_binary_peak_values(&raw, peaks_tag, scan_number, path)?;
    values_to_peak_pairs(values, scan_number, path)
}

fn decode_compressed_peaks(
    peaks_text: &str,
    peaks_tag: &str,
    scan_number: u32,
    path: &Path,
) -> Result<Vec<(f64, f64)>, String> {
    use flate2::read::ZlibDecoder;

    let compressed = decode_base64_bytes(peaks_text, scan_number, path)?;
    let mut decoder = ZlibDecoder::new(compressed.as_slice());
    let mut raw = Vec::new();
    decoder
        .read_to_end(&mut raw)
        .map_err(|err| format!("zlib inflate failed for scan {scan_number}: {err}"))?;

    let values = decode_binary_peak_values(&raw, peaks_tag, scan_number, path)?;
    values_to_peak_pairs(values, scan_number, path)
}

fn decode_f32_peak_values(
    raw: &[u8],
    little_endian: bool,
    scan_number: u32,
    path: &Path,
) -> Result<Vec<f64>, String> {
    if raw.len() % 4 != 0 {
        return Err(format!(
            "compressed peaks length {} is not aligned to 4-byte values in scan {scan_number} in {}",
            raw.len(),
            path.display()
        ));
    }
    raw.chunks_exact(4)
        .map(|chunk| {
            let bytes: [u8; 4] = chunk.try_into().expect("chunk length");
            Ok(if little_endian {
                f32::from_le_bytes(bytes) as f64
            } else {
                f32::from_be_bytes(bytes) as f64
            })
        })
        .collect()
}

fn decode_f64_peak_values(
    raw: &[u8],
    little_endian: bool,
    scan_number: u32,
    path: &Path,
) -> Result<Vec<f64>, String> {
    if raw.len() % 8 != 0 {
        return Err(format!(
            "compressed peaks length {} is not aligned to 8-byte values in scan {scan_number} in {}",
            raw.len(),
            path.display()
        ));
    }
    raw.chunks_exact(8)
        .map(|chunk| {
            let bytes: [u8; 8] = chunk.try_into().expect("chunk length");
            Ok(if little_endian {
                f64::from_le_bytes(bytes)
            } else {
                f64::from_be_bytes(bytes)
            })
        })
        .collect()
}

fn extract_attribute(tag: &str, name: &str) -> Option<String> {
    let pattern = format!("{name}=\"");
    let start = tag.find(&pattern)? + pattern.len();
    let end = tag[start..].find('"')? + start;
    Some(tag[start..end].to_string())
}

fn extract_tag_value(body: &str, tag: &str) -> Option<String> {
    let open = format!("<{tag}>");
    let close = format!("</{tag}>");
    let start = body.find(&open)? + open.len();
    let end = body[start..].find(&close)? + start;
    Some(body[start..end].trim().to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/mzxml")
            .join(name)
    }

    #[test]
    fn parses_hexnac_fixture() {
        let scans = parse_scans(&fixture("hexnac_positive.mzXML")).unwrap();
        assert_eq!(scans.len(), 1);
        assert_eq!(scans[0].scan_number, 1);
        assert!((scans[0].precursor_mz - 800.0).abs() < 0.01);
        assert!(scans[0].peaks.iter().any(|(mz, _)| (*mz - 204.0867).abs() < 0.01));
    }

    #[test]
    fn parses_dss_pair_fixture() {
        let scans = parse_scans(&fixture("dss_pair.mzXML")).unwrap();
        assert_eq!(scans.len(), 2);
        assert_eq!(scans[0].precursor_charge, Some(2));
        assert_eq!(scans[1].precursor_charge, Some(2));
    }

    #[test]
    fn parses_no_diagnostic_fixture() {
        let scans = parse_scans(&fixture("no_diagnostic.mzXML")).unwrap();
        assert_eq!(scans.len(), 1);
        assert!(scans[0].peaks.iter().all(|(mz, _)| (*mz - 204.0).abs() > 1.0));
    }

    #[test]
    fn parses_msconvert_retention_time_attribute_and_zlib_peaks() {
        let scans = parse_scans(&fixture("msconvert_style.mzXML")).unwrap();
        assert_eq!(scans.len(), 1);
        assert_eq!(scans[0].scan_number, 42);
        assert!((scans[0].retention_time_min - 20.0).abs() < 0.001);
        assert_eq!(scans[0].precursor_charge, Some(2));
        assert!(scans[0].peaks.iter().any(|(mz, _)| (*mz - 204.0867).abs() < 0.01));
        assert!(scans[0].peaks.iter().any(|(mz, _)| (*mz - 300.0).abs() < 0.01));
    }

    #[test]
    fn parses_asf_uncompressed_base64_peaks() {
        let path = PathBuf::from(
            "/home/user/Nextcloud/postdoc/data/Dinko/ASF_TRFE/260607_LU02_disoic_ASF_DSS_1.c.mzXML",
        );
        if !path.is_file() {
            return;
        }
        let scans = parse_scans(&path).expect("ASF mzXML should parse");
        assert!(!scans.is_empty());
        assert!(scans.iter().any(|s| s.scan_number == 56));
        let scan56 = scans.iter().find(|s| s.scan_number == 56).unwrap();
        assert!(!scan56.peaks.is_empty());
    }

    #[test]
    fn optional_hcg_integration_prefilter() {
        let path = match std::env::var("GLYCOQUEST_HCG_MZXML") {
            Ok(path) => PathBuf::from(path),
            Err(_) => return,
        };
        if !path.is_file() {
            return;
        }
        let scans = parse_scans(&path).expect("hCG mzXML should parse");
        assert!(!scans.is_empty(), "expected MS2 scans in hCG file");
        assert!(
            scans.iter().all(|scan| scan.retention_time_min.is_finite()),
            "every scan should have retention time"
        );
        assert!(
            scans.iter().all(|scan| !scan.peaks.is_empty()),
            "every MS2 scan should have peaks"
        );
    }
}
