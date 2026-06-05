//! Narrow mzXML reader for MS2 scans used by GlycoQuest prefilters.

use std::path::Path;

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

    let precursor_mz = extract_tag_value(scan_body, "precursorMz")
        .or_else(|| extract_tag_value(scan_body, "precursorMZ"))
        .and_then(|v| v.split_whitespace().next()?.parse().ok())
        .ok_or_else(|| {
            format!(
                "scan {scan_number} missing precursorMz in {}",
                path.display()
            )
        })?;

    let precursor_charge = extract_tag_value(scan_body, "precursorCharge")
        .or_else(|| extract_tag_value(scan_body, "precursor charge"))
        .and_then(|v| v.parse().ok());

    let peaks = parse_peaks(scan_body, path, scan_number)?;

    Ok(Some(Ms2Scan {
        scan_number,
        retention_time_min,
        precursor_mz,
        precursor_charge,
        peaks,
    }))
}

fn parse_retention_time(scan_body: &str, scan_open_tag: &str) -> Result<f64, String> {
    if let Some(rt) = extract_attribute(scan_open_tag, "retTime") {
        if let Ok(minutes) = rt.parse::<f64>() {
            return Ok(minutes);
        }
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
    let content_start = start + tag_end + 1;
    let content_end = scan_body[content_start..]
        .find("</peaks>")
        .ok_or_else(|| format!("unclosed peaks in {}", path.display()))?;
    let peaks_text = scan_body[content_start..content_start + content_end].trim();

    if peaks_text.is_empty() {
        return Ok(Vec::new());
    }

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
}
