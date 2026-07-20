// Copyright (c) ETH Zurich, Mateusz Fido

//! Flatten xQuest XML search results to CSV rows.

use std::fs;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, PartialEq)]
pub struct XQuestHit {
    pub xquest_version: Option<String>,
    pub spectrum_id: String,
    pub search_hit_rank: u32,
    pub link_type: String,
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
    pub xlinker_mass: Option<f64>,
    pub xlinkions_matched: Option<String>,
    pub backboneions_matched: Option<String>,
    pub num_matched_ions_alpha: Option<u32>,
    pub num_matched_ions_beta: Option<u32>,
    pub num_matched_common_ions_alpha: Option<u32>,
    pub num_matched_common_ions_beta: Option<u32>,
    pub num_matched_xlink_ions_alpha: Option<u32>,
    pub num_matched_xlink_ions_beta: Option<u32>,
    pub matched_ions: Vec<XQuestMatchedIon>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct XQuestMatchedIon {
    pub ion_type: String,
    pub position: Option<String>,
    pub theoretical_mz: f64,
    pub observed_mz: f64,
    pub delta_mz: Option<f64>,
    pub delta_ppm: Option<f64>,
    pub intensity: Option<f64>,
    pub label: Option<String>,
}

impl Default for XQuestHit {
    fn default() -> Self {
        Self {
            xquest_version: None,
            spectrum_id: String::new(),
            search_hit_rank: 0,
            link_type: "crosslink".to_string(),
            score: 0.0,
            seq1: String::new(),
            seq2: String::new(),
            prot1: String::new(),
            prot2: String::new(),
            topology: String::new(),
            charge: 0,
            precursor_mz: 0.0,
            mr: 0.0,
            precursor_error_ppm: 0.0,
            xlink_position: String::new(),
            xlinker_mass: None,
            xlinkions_matched: None,
            backboneions_matched: None,
            num_matched_ions_alpha: None,
            num_matched_ions_beta: None,
            num_matched_common_ions_alpha: None,
            num_matched_common_ions_beta: None,
            num_matched_xlink_ions_alpha: None,
            num_matched_xlink_ions_beta: None,
            matched_ions: Vec::new(),
        }
    }
}

impl XQuestHit {
    pub fn normalized_link_type(&self) -> String {
        let raw = self.link_type.trim();
        if raw.is_empty() {
            return if self.seq2.trim().is_empty() && self.prot2.trim().is_empty() {
                "monolink".to_string()
            } else {
                "crosslink".to_string()
            };
        }
        match raw.to_ascii_lowercase().as_str() {
            "mono" | "monolink" => "monolink".to_string(),
            "xlink" | "crosslink" => "crosslink".to_string(),
            other => other.to_string(),
        }
    }

    pub fn is_monolink(&self) -> bool {
        self.normalized_link_type() == "monolink"
    }
}

pub fn extract_hits_from_xml(path: &Path) -> Result<Vec<XQuestHit>, String> {
    let content = fs::read_to_string(path).map_err(|err| err.to_string())?;
    let mut hits = Vec::new();
    let mut pos = 0usize;
    let xquest_version = root_xquest_version(&content);

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
        hits.push(parse_search_hit(
            &block,
            parent_spectrum.as_deref(),
            xquest_version.clone(),
        ));

        pos = tag_end;
    }

    Ok(hits)
}

fn parse_search_hit(
    block: &str,
    parent_spectrum: Option<&str>,
    xquest_version: Option<String>,
) -> XQuestHit {
    // The <spectrum_search spectrum="..."> attribute is the scan-bearing spectrum
    // name (e.g. "2173.sample.c_1896.sample.c"); the search_hit's own `id` is the
    // crosslink identity (peptide-peptide-topology), which does not carry a scan.
    // Prefer the parent spectrum so downstream scan lookup/dedup work correctly.
    let spectrum_id = parent_spectrum
        .map(str::to_string)
        .or_else(|| attr(block, "spectrumid"))
        .or_else(|| attr(block, "id"))
        .unwrap_or_default();
    let seq1 = attr(block, "seq1").unwrap_or_default();
    let seq2 = attr(block, "seq2").unwrap_or_default();
    let prot1 = attr(block, "prot1").unwrap_or_default();
    let prot2 = attr(block, "prot2").unwrap_or_default();

    XQuestHit {
        xquest_version,
        spectrum_id,
        search_hit_rank: attr(block, "search_hit_rank")
            .and_then(|v| v.parse().ok())
            .unwrap_or(0),
        link_type: normalize_link_type(attr(block, "type"), &seq2, &prot2),
        score: attr(block, "score")
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.0),
        seq1,
        seq2,
        prot1,
        prot2,
        topology: attr(block, "topology").unwrap_or_default(),
        charge: attr(block, "charge")
            .and_then(|v| v.parse().ok())
            .unwrap_or(0),
        precursor_mz: attr(block, "mz")
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.0),
        mr: attr(block, "Mr")
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.0),
        precursor_error_ppm: attr(block, "precursorerrorppm")
            .or_else(|| attr(block, "precursorerror"))
            .or_else(|| attr(block, "error"))
            .or_else(|| attr(block, "error_rel"))
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.0),
        xlink_position: attr(block, "xlinkposition").unwrap_or_default(),
        xlinker_mass: attr(block, "xlinkermass")
            .or_else(|| attr(block, "xlinker_mass"))
            .and_then(|v| v.parse().ok()),
        xlinkions_matched: attr(block, "xlinkions_matched"),
        backboneions_matched: attr(block, "backboneions_matched"),
        num_matched_ions_alpha: attr_u32(block, "num_of_matched_ions_alpha"),
        num_matched_ions_beta: attr_u32(block, "num_of_matched_ions_beta"),
        num_matched_common_ions_alpha: attr_u32(block, "num_of_matched_common_ions_alpha"),
        num_matched_common_ions_beta: attr_u32(block, "num_of_matched_common_ions_beta"),
        num_matched_xlink_ions_alpha: attr_u32(block, "num_of_matched_xlink_ions_alpha"),
        num_matched_xlink_ions_beta: attr_u32(block, "num_of_matched_xlink_ions_beta"),
        matched_ions: parse_matched_ions(block),
    }
}

fn root_xquest_version(content: &str) -> Option<String> {
    let start = content.find("<xquest_results")?;
    let end = content[start..].find('>')? + start;
    attr(&content[start..end], "xquest_version")
}

fn attr_u32(tag: &str, name: &str) -> Option<u32> {
    attr(tag, name).and_then(|value| value.parse().ok())
}

fn attr_f64(tag: &str, name: &str) -> Option<f64> {
    attr(tag, name).and_then(|value| value.parse().ok())
}

fn parse_matched_ions(block: &str) -> Vec<XQuestMatchedIon> {
    let mut rows = Vec::new();
    let mut pos = 0usize;
    while let Some(start) = block[pos..].find("<matched_ion").map(|i| pos + i) {
        if block[start..].starts_with("<matched_ions") {
            pos = start + "<matched_ions".len();
            continue;
        }
        let Some(end) = block[start..].find('>').map(|i| start + i + 1) else {
            break;
        };
        let tag = &block[start..end];
        let Some(theoretical_mz) =
            attr_f64(tag, "theoretical_mz").or_else(|| attr_f64(tag, "ion_mz"))
        else {
            pos = end;
            continue;
        };
        let Some(observed_mz) = attr_f64(tag, "observed_mz").or_else(|| attr_f64(tag, "peak_mz"))
        else {
            pos = end;
            continue;
        };
        rows.push(XQuestMatchedIon {
            ion_type: attr(tag, "ion_type").unwrap_or_else(|| "unknown".to_string()),
            position: attr(tag, "position"),
            theoretical_mz,
            observed_mz,
            delta_mz: attr_f64(tag, "delta_mz"),
            delta_ppm: attr_f64(tag, "delta_ppm"),
            intensity: attr_f64(tag, "intensity"),
            label: attr(tag, "label"),
        });
        pos = end;
    }
    rows
}

fn normalize_link_type(value: Option<String>, seq2: &str, prot2: &str) -> String {
    let raw = value.unwrap_or_default();
    match raw.trim().to_ascii_lowercase().as_str() {
        "" if seq2.trim().is_empty() && prot2.trim().is_empty() => "monolink".to_string(),
        "" | "xlink" | "crosslink" => "crosslink".to_string(),
        "mono" | "monolink" => "monolink".to_string(),
        other => other.to_string(),
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
        "spectrum_id\trank\tlink_type\tscore\tseq1\tseq2\tprot1\tprot2\ttopology\tcharge\tprecursor_mz\tmr\tprecursor_error_ppm\txlink_position\txlinker_mass\tpostfilter_status"
            .to_string(),
    ];
    for hit in hits {
        lines.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\traw",
            hit.spectrum_id,
            hit.search_hit_rank,
            hit.normalized_link_type(),
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
            hit.xlinker_mass
                .map(|value| value.to_string())
                .unwrap_or_default(),
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

    #[test]
    fn parses_monolink_type_and_crosslinker_mass() {
        let xml = r#"<spectrum_search spectrum="42.run.c"><search_hit search_hit_rank="1" type="monolink" score="9.5" seq1="PEPXK" seq2="" prot1="FETUA_BOVIN" prot2="" topology="" xlinkposition="5" xlinkermass="156.07864" charge="4" mz="800.1" Mr="3199.2" error_rel="2.2"/></spectrum_search>"#;
        let dir = std::env::temp_dir().join(format!("glycoquest_mono_xml_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("result.xml"), xml).unwrap();

        let hits = extract_hits_from_xml(&dir.join("result.xml")).unwrap();

        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].spectrum_id, "42.run.c");
        assert_eq!(hits[0].link_type, "monolink");
        assert_eq!(hits[0].seq2, "");
        assert_eq!(hits[0].prot2, "");
        assert_eq!(hits[0].xlink_position, "5");
        assert!((hits[0].xlinker_mass.unwrap() - 156.07864).abs() < 1e-6);
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn parses_xquest_v217_aggregate_and_exact_matched_ions() {
        let xml = r#"<xquest_results xquest_version="xquest 2.1.7">
<spectrum_search spectrum="12705.run.c">
<search_hit search_hit_rank="1" type="xlink" score="12.94" seq1="VHTECCHGDLLECADDRADLAKYICENQDSISSK" seq2="MVLSPADKTXVKAAWGK" prot1="ALBU" prot2="HBA" topology="inter" xlinkposition="22,8" charge="6" mz="1470.96960" Mr="8820.1" error_rel="0.9" xlinkions_matched="2/80" backboneions_matched="7/150" num_of_matched_common_ions_alpha="4" num_of_matched_common_ions_beta="3" num_of_matched_xlink_ions_alpha="1" num_of_matched_xlink_ions_beta="1">
  <matched_ions source="xquest_v2.1.7">
    <matched_ion ion_type="common_b" position="2" theoretical_mz="244.1770" observed_mz="244.1780" delta_mz="-0.0010" delta_ppm="-4.1" intensity="12345.0"/>
  </matched_ions>
</search_hit>
</spectrum_search>
</xquest_results>"#;
        let dir =
            std::env::temp_dir().join(format!("glycoquest_xquest_v217_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("result.xml"), xml).unwrap();

        let hits = extract_hits_from_xml(&dir.join("result.xml")).unwrap();

        assert_eq!(hits.len(), 1);
        let hit = &hits[0];
        assert_eq!(hit.xquest_version.as_deref(), Some("xquest 2.1.7"));
        assert_eq!(hit.xlinkions_matched.as_deref(), Some("2/80"));
        assert_eq!(hit.backboneions_matched.as_deref(), Some("7/150"));
        assert_eq!(hit.num_matched_common_ions_alpha, Some(4));
        assert_eq!(hit.num_matched_xlink_ions_beta, Some(1));
        assert_eq!(hit.matched_ions.len(), 1);
        assert_eq!(hit.matched_ions[0].ion_type, "common_b");
        assert_eq!(hit.matched_ions[0].position.as_deref(), Some("2"));
        assert!((hit.matched_ions[0].observed_mz - 244.1780).abs() < 1e-6);
        assert!((hit.matched_ions[0].intensity.unwrap() - 12345.0).abs() < 1e-6);

        let _ = std::fs::remove_dir_all(dir);
    }
}
