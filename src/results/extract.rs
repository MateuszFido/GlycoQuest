//! Flatten xQuest XML search results to CSV rows.

use std::fs;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, PartialEq, Default)]
pub struct XQuestHit {
    pub spectrum_id: String,
    pub score: f64,
    pub seq1: String,
    pub seq2: String,
    pub precursor_error_ppm: f64,
    pub xlink_position: String,
}

pub fn extract_hits_from_xml(path: &Path) -> Result<Vec<XQuestHit>, String> {
    let content = fs::read_to_string(path).map_err(|err| err.to_string())?;
    let mut hits = Vec::new();
    let mut pos = 0usize;

    while let Some(start) = content[pos..].find("<search_hit").map(|i| pos + i) {
        let tag_end = content[start..]
            .find('>')
            .ok_or_else(|| format!("malformed search_hit in {}", path.display()))?
            + start
            + 1;
        let open_tag = &content[start..tag_end];
        let block = if open_tag.ends_with("/>") {
            open_tag.to_string()
        } else {
            let close = content[tag_end..]
                .find("</search_hit>")
                .ok_or_else(|| format!("unclosed search_hit in {}", path.display()))?;
            content[start..tag_end + close].to_string()
        };

        hits.push(XQuestHit {
            spectrum_id: attr(&block, "spectrumid").unwrap_or_default(),
            score: attr(&block, "score")
                .and_then(|v| v.parse().ok())
                .unwrap_or(0.0),
            seq1: attr(&block, "seq1").unwrap_or_default(),
            seq2: attr(&block, "seq2").unwrap_or_default(),
            precursor_error_ppm: attr(&block, "precursorerrorppm")
                .or_else(|| attr(&block, "precursorerror"))
                .and_then(|v| v.parse().ok())
                .unwrap_or(0.0),
            xlink_position: attr(&block, "xlinkposition").unwrap_or_default(),
        });

        pos = tag_end;
    }

    Ok(hits)
}

fn attr(tag: &str, name: &str) -> Option<String> {
    for needle in [
        format!(" {name}=\""),
        format!("\n{name}=\""),
        format!("\t{name}=\""),
        format!("<search_hit {name}=\""),
    ] {
        if let Some(start) = tag.find(&needle) {
            let value_start = start + needle.len();
            let end = tag[value_start..].find('"')? + value_start;
            return Some(tag[value_start..end].to_string());
        }
    }
    None
}

pub fn write_hits_csv(path: &Path, hits: &[XQuestHit]) -> Result<(), String> {
    let mut lines = vec![
        "spectrum_id\tscore\tseq1\tseq2\tprecursor_error_ppm\txlink_position\tpostfilter_status"
            .to_string(),
    ];
    for hit in hits {
        lines.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            hit.spectrum_id,
            hit.score,
            hit.seq1,
            hit.seq2,
            hit.precursor_error_ppm,
            hit.xlink_position,
            "raw",
        ));
    }
    fs::write(path, lines.join("\n") + "\n").map_err(|err| err.to_string())
}

pub fn find_result_xmls(jobs_dir: &Path) -> Result<Vec<PathBuf>, String> {
    let mut paths = Vec::new();
    if !jobs_dir.is_dir() {
        return Ok(paths);
    }
    for entry in fs::read_dir(jobs_dir).map_err(|err| err.to_string())? {
        let entry = entry.map_err(|err| err.to_string())?;
        if !entry.file_type().map_err(|err| err.to_string())?.is_dir() {
            continue;
        }
        let job_dir = entry.path();
        for candidate in [
            job_dir.join("results/xquest.xml"),
            job_dir.join("result.xml"),
        ] {
            if candidate.is_file() {
                paths.push(candidate);
                break;
            }
        }
    }
    paths.sort();
    Ok(paths)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_minimal_search_hit() {
        let xml = r#"<search_hit spectrumid="1" score="12.3" seq1="PEPTXIDE" seq2="PEPTIDE" precursorerrorppm="4.5" xlinkposition="3-7"/>"#;
        let dir = std::env::temp_dir().join(format!("glycoquest_xml_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("result.xml"), xml).unwrap();
        let hits = extract_hits_from_xml(&dir.join("result.xml")).unwrap();
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].score, 12.3);
        assert_eq!(hits[0].seq1, "PEPTXIDE");
        let _ = std::fs::remove_dir_all(dir);
    }
}
