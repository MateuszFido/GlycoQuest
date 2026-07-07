//! Flatten xQuest XML search results to CSV rows.

use std::fs;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, PartialEq, Default)]
pub struct XQuestHit {
    pub spectrum_id: String,
    pub search_hit_rank: u32,
    pub score: f64,
    pub seq1: String,
    pub seq2: String,
    pub prot1: String,
    pub prot2: String,
    pub topology: String,
    pub charge: u8,
    pub precursor_mz: f64,
    pub mr: f64,
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

        let parent_spectrum = parent_spectrum(&content, start);
        hits.push(parse_search_hit(&block, parent_spectrum.as_deref()));

        pos = tag_end;
    }

    Ok(hits)
}

fn parse_search_hit(block: &str, parent_spectrum: Option<&str>) -> XQuestHit {
    // The <spectrum_search spectrum="..."> attribute is the scan-bearing spectrum
    // name (e.g. "2173.sample.c_1896.sample.c"); the search_hit's own `id` is the
    // crosslink identity (peptide-peptide-topology), which does not carry a scan.
    // Prefer the parent spectrum so downstream scan lookup/dedup work correctly.
    let spectrum_id = parent_spectrum
        .map(str::to_string)
        .or_else(|| attr(block, "spectrumid"))
        .or_else(|| attr(block, "id"))
        .unwrap_or_default();

    XQuestHit {
        spectrum_id,
        search_hit_rank: attr(block, "search_hit_rank")
            .and_then(|v| v.parse().ok())
            .unwrap_or(0),
        score: attr(block, "score").and_then(|v| v.parse().ok()).unwrap_or(0.0),
        seq1: attr(block, "seq1").unwrap_or_default(),
        seq2: attr(block, "seq2").unwrap_or_default(),
        prot1: attr(block, "prot1").unwrap_or_default(),
        prot2: attr(block, "prot2").unwrap_or_default(),
        topology: attr(block, "topology").unwrap_or_default(),
        charge: attr(block, "charge")
            .and_then(|v| v.parse().ok())
            .unwrap_or(0),
        precursor_mz: attr(block, "mz").and_then(|v| v.parse().ok()).unwrap_or(0.0),
        mr: attr(block, "Mr").and_then(|v| v.parse().ok()).unwrap_or(0.0),
        precursor_error_ppm: attr(block, "precursorerrorppm")
            .or_else(|| attr(block, "precursorerror"))
            .or_else(|| attr(block, "error"))
            .or_else(|| attr(block, "error_rel"))
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.0),
        xlink_position: attr(block, "xlinkposition").unwrap_or_default(),
    }
}

fn parent_spectrum(content: &str, hit_start: usize) -> Option<String> {
    let before = &content[..hit_start];
    let tag_start = before.rfind("<spectrum_search")?;
    let tag_end = before[tag_start..].find('>')? + tag_start;
    attr(&before[tag_start..tag_end], "spectrum")
}

fn attr(tag: &str, name: &str) -> Option<String> {
    for needle in [
        format!(" {name}=\""),
        format!("\n{name}=\""),
        format!("\t{name}=\""),
        format!("<search_hit {name}=\""),
        format!("<spectrum_search {name}=\""),
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
        "spectrum_id\trank\tscore\tseq1\tseq2\tprot1\tprot2\ttopology\tcharge\tprecursor_mz\tmr\tprecursor_error_ppm\txlink_position\tpostfilter_status"
            .to_string(),
    ];
    for hit in hits {
        lines.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\traw",
            hit.spectrum_id,
            hit.search_hit_rank,
            hit.score,
            hit.seq1,
            hit.seq2,
            hit.prot1,
            hit.prot2,
            hit.topology,
            hit.charge,
            hit.precursor_mz,
            hit.mr,
            hit.precursor_error_ppm,
            hit.xlink_position,
        ));
    }
    fs::write(path, lines.join("\n") + "\n").map_err(|err| err.to_string())
}

/// Find each job's result XML, paired with the job id (the job directory name)
/// so hits can be tied back to the job manifest.
pub fn find_result_xmls(jobs_dir: &Path) -> Result<Vec<(String, PathBuf)>, String> {
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
        let job_id = job_dir
            .file_name()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_default();
        for candidate in [
            job_dir.join("results/xquest.xml"),
            job_dir.join("result.xml"),
        ] {
            if candidate.is_file() {
                paths.push((job_id.clone(), candidate));
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
        assert_eq!(hits[0].spectrum_id, "1");
        assert_eq!(hits[0].score, 12.3);
        assert_eq!(hits[0].seq1, "PEPTXIDE");
        let _ = std::fs::remove_dir_all(dir);
    }
}
