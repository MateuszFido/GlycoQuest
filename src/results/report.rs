// Copyright (c) ETH Zurich, Mateusz Fido

//! CLMS-CSV-compatible export and a self-contained HTML QC/glycan report.
//!
//! Two artifacts are written next to `glycoquest_xquest.csv`:
//!
//! * `results/xiview.csv` — passing crosslinks in the CLMS-CSV column layout,
//!   with glycan annotations appended as trailing columns.
//! * `results/report.html` — a dependency-free HTML report (inline CSS + SVG
//!   charts, no network or JavaScript) summarizing the prefilter funnel, hit
//!   quality, and glycan distribution.

use std::fmt::Write as _;
use std::path::Path;

use crate::crosslinker::CrosslinkerProfile;
use crate::fasta::FastaDatabase;
use crate::jobs::JobManifest;
use crate::prefilter::PrefilterStats;

use super::filter::{AnnotatedHit, PostfilterStatus, format_glyco_sites};
use super::mapping::{
    abs_position, first_protein, locate_peptide, parse_link_positions, protein_lookup,
    resolve_peptide,
};

/// Metadata rendered in the HTML report header.
#[derive(Debug, Clone)]
pub struct ReportContext {
    pub project: String,
    pub input_label: String,
    pub crosslinker_name: String,
    pub xlink_sites: String,
    pub glycan_library: String,
}

impl ReportContext {
    pub fn new(
        project: String,
        input_label: String,
        crosslinker: &CrosslinkerProfile,
        glycan_library: String,
    ) -> Self {
        Self {
            project,
            input_label,
            crosslinker_name: crosslinker.name.clone(),
            xlink_sites: crosslinker.xlink_sites.clone(),
            glycan_library,
        }
    }
}

/// Counts returned by [`write_xiview_csv`] for logging.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct XiviewSummary {
    /// Passing crosslinks written to the CSV.
    pub rows: usize,
    /// Rows for which both peptides mapped to a protein sequence (absolute positions resolved).
    pub mapped: usize,
}

// -----------------------------------------------------------------------------
// CLMS-CSV export
// -----------------------------------------------------------------------------

/// Write passing crosslinks to `path` in CLMS-CSV layout.
///
/// Peptide sequences are resolved back from xQuest pseudo-residues (via the job
/// manifest) so they can be located in the FASTA to yield 1-based peptide start
/// positions (`PepPos*`) and absolute cross-link residue positions (`AbsPos*`).
/// Glycan annotations are appended as trailing columns; conforming importers ignore
/// unrecognized columns.
pub fn write_xiview_csv(
    path: &Path,
    hits: &[AnnotatedHit],
    manifest: Option<&JobManifest>,
    fasta: &FastaDatabase,
) -> Result<XiviewSummary, String> {
    let proteins = protein_lookup(fasta);

    let header = [
        "Protein1",
        "PepPos1",
        "PepSeq1",
        "LinkPos1",
        "AbsPos1",
        "Protein2",
        "PepPos2",
        "PepSeq2",
        "LinkPos2",
        "AbsPos2",
        "Score",
        "Id",
        "Scan",
        "Glycan",
        "GlycoResidue",
        "GlycoSites",
        "AllSitesPlausible",
        "LossLabel",
        "LinkType",
        "XlinkerMass",
    ]
    .join(",");
    let mut lines = vec![header];
    let mut summary = XiviewSummary::default();

    for hit in hits
        .iter()
        .filter(|h| h.postfilter_status == PostfilterStatus::Pass)
    {
        let plan = manifest
            .and_then(|m| m.by_job_id(&hit.job_id))
            .map(|entry| &entry.varmod_plan);

        let pep1 = resolve_peptide(&hit.hit.seq1, plan);
        let pep2 = resolve_peptide(&hit.hit.seq2, plan);
        let (link1, link2) = parse_link_positions(&hit.hit.xlink_position);

        let prot1 = first_protein(&hit.hit.prot1);
        let prot2 = first_protein(&hit.hit.prot2);
        let pep_pos1 = locate_peptide(&proteins, prot1, &pep1);
        let pep_pos2 = locate_peptide(&proteins, prot2, &pep2);

        let abs1 = abs_position(pep_pos1, link1);
        let abs2 = abs_position(pep_pos2, link2);
        let is_monolink = hit.hit.is_monolink();
        if abs1.is_some() && (is_monolink || abs2.is_some()) {
            summary.mapped += 1;
        }

        let id = format!(
            "{}:{}",
            hit.scan.map(|s| s.to_string()).unwrap_or_default(),
            hit.job_id
        );

        let row = [
            csv_field(prot1),
            opt_num(pep_pos1),
            csv_field(&pep1),
            opt_num(link1),
            opt_num(abs1),
            csv_field(prot2),
            opt_num(pep_pos2),
            csv_field(&pep2),
            opt_num(link2),
            opt_num(abs2),
            format!("{}", hit.hit.score),
            csv_field(&id),
            hit.scan.map(|s| s.to_string()).unwrap_or_default(),
            csv_field(hit.glycan_composition.as_deref().unwrap_or("")),
            hit.glyco_residue.map(|c| c.to_string()).unwrap_or_default(),
            csv_field(&format_glyco_sites(&hit.glyco_sites)),
            hit.all_sites_plausible.to_string(),
            csv_field(hit.loss_label.as_deref().unwrap_or("")),
            csv_field(&hit.hit.normalized_link_type()),
            opt_float(hit.hit.xlinker_mass),
        ]
        .join(",");
        lines.push(row);
        summary.rows += 1;
    }

    std::fs::write(path, lines.join("\n") + "\n").map_err(|err| err.to_string())?;
    Ok(summary)
}

fn csv_field(value: &str) -> String {
    if value.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}

fn opt_num(value: Option<usize>) -> String {
    value.map(|v| v.to_string()).unwrap_or_default()
}

fn opt_float(value: Option<f64>) -> String {
    value.map(|v| format!("{v:.6}")).unwrap_or_default()
}

// -----------------------------------------------------------------------------
// HTML report
// -----------------------------------------------------------------------------

/// Write a self-contained HTML QC/glycan report to `path`.
pub fn write_html_report(
    path: &Path,
    hits: &[AnnotatedHit],
    stats: &PrefilterStats,
    ctx: &ReportContext,
) -> Result<(), String> {
    let total = hits.len();
    let pass = hits
        .iter()
        .filter(|h| h.postfilter_status == PostfilterStatus::Pass)
        .count();

    let mut body = String::new();
    body.push_str(&summary_section(ctx, total, pass));

    if stats.scans_total > 0 {
        body.push_str(&prefilter_funnel(stats));
    }

    body.push_str("<div class=\"grid\">");
    body.push_str(&hard_status_chart(hits));
    body.push_str(&glycan_chart(hits));
    body.push_str(&glyco_residue_chart(hits));
    body.push_str(&score_histogram(hits));
    body.push_str(&error_histogram(hits));
    body.push_str("</div>");

    body.push_str(&hits_table(hits));

    let html = format!(
        "<!DOCTYPE html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">\
<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\
<title>GlycoQuest report — {project}</title>{style}</head><body>\
<h1>GlycoQuest report</h1>{body}\
<footer>Generated by GlycoQuest. Open <code>viewer/index.html</code> for the interactive crosslink viewer, or import <code>xiview.csv</code> into a CLMS-CSV-compatible network tool.</footer>\
</body></html>\n",
        project = esc(&ctx.project),
        style = STYLE,
        body = body,
    );

    std::fs::write(path, html).map_err(|err| err.to_string())
}

fn summary_section(ctx: &ReportContext, total: usize, pass: usize) -> String {
    let fail = total.saturating_sub(pass);
    format!(
        "<section class=\"summary\"><div class=\"cards\">\
{card_pass}{card_fail}{card_total}</div>\
<table class=\"meta\">\
<tr><th>Project</th><td>{project}</td></tr>\
<tr><th>Input</th><td>{input}</td></tr>\
<tr><th>Crosslinker</th><td>{xl} (sites: {sites})</td></tr>\
<tr><th>Glycan library</th><td>{glycans}</td></tr>\
</table></section>",
        card_pass = card("Passing hits", pass, "pass"),
        card_fail = card("Filtered out", fail, "fail"),
        card_total = card("Total hits", total, "total"),
        project = esc(&ctx.project),
        input = esc(&ctx.input_label),
        xl = esc(&ctx.crosslinker_name),
        sites = esc(&ctx.xlink_sites),
        glycans = esc(&ctx.glycan_library),
    )
}

fn card(label: &str, value: usize, class: &str) -> String {
    format!(
        "<div class=\"card {class}\"><div class=\"num\">{value}</div>\
<div class=\"lbl\">{label}</div></div>",
        class = class,
        value = value,
        label = esc(label),
    )
}

fn prefilter_funnel(stats: &PrefilterStats) -> String {
    let data = vec![
        ("MS/MS scans".to_string(), stats.scans_total as f64),
        (
            "Diagnostic-ion positive".to_string(),
            stats.diagnostic_positive as f64,
        ),
        ("Isotope pairs".to_string(), stats.isotope_pairs as f64),
        ("Passed to xQuest".to_string(), stats.filtered_scans as f64),
    ];
    svg_hbar("Prefilter funnel", &data)
}

fn hard_status_chart(hits: &[AnnotatedHit]) -> String {
    let mut counts: std::collections::BTreeMap<&'static str, usize> =
        std::collections::BTreeMap::new();
    for hit in hits {
        *counts.entry(hit.hard_status.as_str()).or_insert(0) += 1;
    }
    let data: Vec<(String, f64)> = counts
        .into_iter()
        .map(|(label, count)| (label.to_string(), count as f64))
        .collect();
    svg_hbar("Post-filter outcome", &data)
}

fn glycan_chart(hits: &[AnnotatedHit]) -> String {
    let mut counts: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    for hit in hits
        .iter()
        .filter(|h| h.postfilter_status == PostfilterStatus::Pass)
    {
        if let Some(comp) = &hit.glycan_composition {
            *counts.entry(comp.clone()).or_insert(0) += 1;
        }
    }
    let mut data: Vec<(String, f64)> = counts.into_iter().map(|(k, v)| (k, v as f64)).collect();
    data.sort_by(|a, b| {
        b.1.partial_cmp(&a.1)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.0.cmp(&b.0))
    });
    data.truncate(15);
    svg_hbar("Glycans on passing hits (top 15)", &data)
}

fn glyco_residue_chart(hits: &[AnnotatedHit]) -> String {
    let mut counts: std::collections::BTreeMap<char, usize> = std::collections::BTreeMap::new();
    for hit in hits
        .iter()
        .filter(|h| h.postfilter_status == PostfilterStatus::Pass)
    {
        for site in &hit.glyco_sites {
            *counts.entry(site.residue).or_insert(0) += 1;
        }
    }
    let data: Vec<(String, f64)> = counts
        .into_iter()
        .map(|(residue, count)| (glyco_residue_label(residue), count as f64))
        .collect();
    svg_hbar("Glycosylation site (passing hits)", &data)
}

fn glyco_residue_label(residue: char) -> String {
    match residue {
        'N' => "N (N-linked)".to_string(),
        'S' => "S (O-linked)".to_string(),
        'T' => "T (O-linked)".to_string(),
        other => other.to_string(),
    }
}

fn score_histogram(hits: &[AnnotatedHit]) -> String {
    let values: Vec<f64> = hits.iter().map(|h| h.hit.score).collect();
    svg_hist("xQuest score distribution", &values, 20)
}

fn error_histogram(hits: &[AnnotatedHit]) -> String {
    let values: Vec<f64> = hits.iter().map(|h| h.hit.precursor_error_ppm).collect();
    svg_hist("Precursor error (ppm)", &values, 20)
}

fn hits_table(hits: &[AnnotatedHit]) -> String {
    const CAP: usize = 500;
    if hits.is_empty() {
        return "<div class=\"chart\"><h3>Hits</h3><p class=\"empty\">no hits</p></div>"
            .to_string();
    }
    let mut s = String::from("<h2>Hits</h2>");
    if hits.len() > CAP {
        let _ = write!(
            s,
            "<p class=\"note\">Showing first {CAP} of {} hits (sorted by soft score). \
             See <code>glycoquest_xquest.csv</code> for the full table.</p>",
            hits.len()
        );
    }
    s.push_str(
        "<div class=\"tablewrap\"><table class=\"hits\"><thead><tr>\
<th>Scan</th><th>Peptide 1</th><th>Peptide 2</th><th>Glycan</th><th>Site</th>\
<th>Sequon</th><th>z</th><th>ppm</th><th>Score</th><th>Soft</th><th>Status</th>\
</tr></thead><tbody>",
    );
    for hit in hits.iter().take(CAP) {
        let status_class = if hit.postfilter_status == PostfilterStatus::Pass {
            "ok"
        } else {
            "no"
        };
        let sequon = match hit.sequon_present {
            Some(true) => "yes",
            Some(false) => "no",
            None => "-",
        };
        let site = if hit.glyco_sites.is_empty() {
            "-".to_string()
        } else {
            format_glyco_sites(&hit.glyco_sites)
        };
        let _ = write!(
            s,
            "<tr><td>{scan}</td><td class=\"seq\">{seq1}</td><td class=\"seq\">{seq2}</td>\
<td>{glycan}</td><td>{site}</td><td>{sequon}</td><td>{charge}</td>\
<td>{ppm:.1}</td><td>{score:.2}</td><td>{soft:.2}</td>\
<td class=\"st {status_class}\">{status}</td></tr>",
            scan = hit.scan.map(|s| s.to_string()).unwrap_or_default(),
            seq1 = esc(&hit.hit.seq1),
            seq2 = esc(&hit.hit.seq2),
            glycan = esc(hit.glycan_composition.as_deref().unwrap_or("-")),
            site = esc(&site),
            sequon = sequon,
            charge = hit.hit.charge,
            ppm = hit.hit.precursor_error_ppm,
            score = hit.hit.score,
            soft = hit.soft_score,
            status_class = status_class,
            status = hit.postfilter_status.as_str(),
        );
    }
    s.push_str("</tbody></table></div>");
    s
}

// -----------------------------------------------------------------------------
// SVG chart primitives (pure Rust, no dependencies)
// -----------------------------------------------------------------------------

fn svg_hbar(title: &str, data: &[(String, f64)]) -> String {
    if data.is_empty() {
        return format!(
            "<div class=\"chart\"><h3>{}</h3><p class=\"empty\">no data</p></div>",
            esc(title)
        );
    }
    let max = data
        .iter()
        .map(|(_, v)| *v)
        .fold(0.0_f64, f64::max)
        .max(1.0);
    let row_h = 26.0;
    let label_w = 190.0;
    let bar_w = 340.0;
    let val_w = 60.0;
    let width = label_w + bar_w + val_w;
    let height = row_h * data.len() as f64 + 8.0;

    let mut s = format!(
        "<div class=\"chart\"><h3>{}</h3>\
<svg viewBox=\"0 0 {width:.0} {height:.0}\" preserveAspectRatio=\"xMinYMin meet\" \
class=\"bar\" role=\"img\">",
        esc(title)
    );
    for (i, (label, value)) in data.iter().enumerate() {
        let y = i as f64 * row_h + 4.0;
        let w = (value / max) * bar_w;
        let _ = write!(
            s,
            "<text x=\"{lx:.0}\" y=\"{ty:.0}\" class=\"lab\" text-anchor=\"end\">{label}</text>\
<rect x=\"{bx:.0}\" y=\"{by:.1}\" width=\"{w:.1}\" height=\"{bh:.0}\" class=\"barfill\"/>\
<text x=\"{vx:.1}\" y=\"{ty:.0}\" class=\"val\">{value}</text>",
            lx = label_w - 6.0,
            ty = y + 17.0,
            label = esc(label),
            bx = label_w,
            by = y + 4.0,
            w = w,
            bh = row_h - 9.0,
            vx = label_w + w + 5.0,
            value = fmt_num(*value),
        );
    }
    s.push_str("</svg></div>");
    s
}

fn svg_hist(title: &str, values: &[f64], bins: usize) -> String {
    let finite: Vec<f64> = values.iter().copied().filter(|v| v.is_finite()).collect();
    if finite.is_empty() {
        return format!(
            "<div class=\"chart\"><h3>{}</h3><p class=\"empty\">no data</p></div>",
            esc(title)
        );
    }
    let min = finite.iter().copied().fold(f64::INFINITY, f64::min);
    let max = finite.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let span = (max - min).max(1e-9);
    let bins = bins.max(1);
    let mut counts = vec![0usize; bins];
    for &v in &finite {
        let mut idx = (((v - min) / span) * bins as f64) as usize;
        if idx >= bins {
            idx = bins - 1;
        }
        counts[idx] += 1;
    }
    let maxc = counts.iter().copied().max().unwrap_or(1).max(1) as f64;

    let width = 560.0;
    let plot_h = 170.0;
    let pad_l = 34.0;
    let pad_b = 22.0;
    let bar_area = width - pad_l - 6.0;
    let step = bar_area / bins as f64;

    let mut s = format!(
        "<div class=\"chart\"><h3>{}</h3>\
<svg viewBox=\"0 0 {width:.0} {total_h:.0}\" preserveAspectRatio=\"xMinYMin meet\" \
class=\"hist\" role=\"img\">",
        esc(title),
        total_h = plot_h + pad_b + 6.0,
    );
    // Axis baseline.
    let _ = write!(
        s,
        "<line x1=\"{pad_l:.0}\" y1=\"{y:.0}\" x2=\"{x2:.0}\" y2=\"{y:.0}\" class=\"axis\"/>",
        pad_l = pad_l,
        y = plot_h,
        x2 = width - 6.0,
    );
    for (i, &count) in counts.iter().enumerate() {
        if count == 0 {
            continue;
        }
        let h = (count as f64 / maxc) * (plot_h - 6.0);
        let x = pad_l + i as f64 * step;
        let _ = write!(
            s,
            "<rect x=\"{x:.1}\" y=\"{y:.1}\" width=\"{w:.1}\" height=\"{h:.1}\" class=\"barfill\">\
<title>{lo:.2}–{hi:.2}: {count}</title></rect>",
            x = x + 0.5,
            y = plot_h - h,
            w = (step - 1.0).max(0.5),
            h = h,
            lo = min + span * (i as f64 / bins as f64),
            hi = min + span * ((i + 1) as f64 / bins as f64),
            count = count,
        );
    }
    // Min / max labels.
    let _ = write!(
        s,
        "<text x=\"{pad_l:.0}\" y=\"{ly:.0}\" class=\"tick\">{min:.2}</text>\
<text x=\"{x2:.0}\" y=\"{ly:.0}\" class=\"tick\" text-anchor=\"end\">{max:.2}</text>\
<text x=\"{pad_l:.0}\" y=\"12\" class=\"tick\">n={n}</text>",
        pad_l = pad_l,
        ly = plot_h + pad_b - 4.0,
        min = min,
        x2 = width - 6.0,
        max = max,
        n = finite.len(),
    );
    s.push_str("</svg></div>");
    s
}

fn fmt_num(value: f64) -> String {
    if (value.fract()).abs() < 1e-9 {
        format!("{}", value as i64)
    } else {
        format!("{value:.1}")
    }
}

fn esc(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

const STYLE: &str = "<style>\
:root{--bg:#0f172a;--card:#1e293b;--fg:#e2e8f0;--muted:#94a3b8;--accent:#38bdf8;\
--ok:#22c55e;--no:#f87171;--bar:#38bdf8}\
*{box-sizing:border-box}body{margin:0;padding:24px;font:14px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif;\
background:var(--bg);color:var(--fg)}h1{margin:0 0 16px;font-size:22px}\
h2{margin:28px 0 10px;font-size:17px}h3{margin:0 0 8px;font-size:14px;color:var(--muted)}\
.summary{display:flex;flex-wrap:wrap;gap:20px;align-items:flex-start}\
.cards{display:flex;gap:12px}.card{background:var(--card);border-radius:10px;padding:14px 20px;min-width:110px;text-align:center}\
.card .num{font-size:28px;font-weight:700}.card .lbl{color:var(--muted);font-size:12px;margin-top:2px}\
.card.pass .num{color:var(--ok)}.card.fail .num{color:var(--no)}.card.total .num{color:var(--accent)}\
table.meta{background:var(--card);border-radius:10px;border-collapse:collapse;overflow:hidden}\
table.meta th,table.meta td{padding:6px 14px;text-align:left;border-bottom:1px solid #334155}\
table.meta th{color:var(--muted);font-weight:500}\
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:16px;margin-top:20px}\
.chart{background:var(--card);border-radius:10px;padding:14px}\
.chart .empty{color:var(--muted);font-style:italic}\
svg.bar,svg.hist{width:100%;height:auto}\
.barfill{fill:var(--bar)}.axis{stroke:#475569;stroke-width:1}\
text.lab{fill:var(--fg);font-size:12px}text.val{fill:var(--muted);font-size:12px}\
text.tick{fill:var(--muted);font-size:11px}\
.note{color:var(--muted)}.tablewrap{overflow-x:auto;background:var(--card);border-radius:10px}\
table.hits{border-collapse:collapse;width:100%;font-size:13px}\
table.hits th,table.hits td{padding:6px 10px;border-bottom:1px solid #334155;text-align:left;white-space:nowrap}\
table.hits th{color:var(--muted);position:sticky;top:0;background:var(--card)}\
td.seq{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}\
td.st.ok{color:var(--ok)}td.st.no{color:var(--no)}\
footer{margin-top:28px;color:var(--muted);font-size:12px}\
footer a{color:var(--accent)}code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}\
</style>";

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::{FastaDatabase, FastaEntry};
    use crate::results::extract::XQuestHit;
    use std::path::PathBuf;

    fn fasta() -> FastaDatabase {
        FastaDatabase {
            path: PathBuf::from("test.fasta"),
            entries: vec![
                FastaEntry {
                    header: "sp|P00761|TRYP_PIG trypsin".into(),
                    sequence: "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHR".into(),
                },
                FastaEntry {
                    header: "HRP".into(),
                    sequence: "QLTPTFYDNSCPNVSNIVRDTIVNELR".into(),
                },
            ],
        }
    }

    fn pass_hit(seq1: &str, seq2: &str, prot1: &str, prot2: &str, xlink: &str) -> AnnotatedHit {
        AnnotatedHit {
            hit: XQuestHit {
                seq1: seq1.into(),
                seq2: seq2.into(),
                prot1: prot1.into(),
                prot2: prot2.into(),
                xlink_position: xlink.into(),
                charge: 4,
                score: 12.5,
                precursor_error_ppm: 2.1,
                ..Default::default()
            },
            job_id: "HexNAc_1_".into(),
            source_file: None,
            scan: Some(101),
            glycan_name: Some("HexNAc(1)".into()),
            glycan_composition: Some("HexNAc(1)".into()),
            glycan_mass: Some(203.079),
            loss_label: Some("none".into()),
            glyco_residue: Some('N'),
            glyco_peptide: Some(2),
            glyco_sites: vec![super::super::filter::GlycoSite {
                peptide: 2,
                peptide_position: 1,
                residue: 'N',
                sequon_present: Some(true),
                plausible: true,
            }],
            all_sites_plausible: true,
            n_glycan_pseudo: 1,
            matched_families: vec!["HexNAc".into()],
            matched_ions: vec![],
            matched_ion_count: 3,
            sequon_present: Some(true),
            charge_plausible: true,
            hard_status: super::super::filter::HardStatus::Pass,
            soft_score: 13.6,
            postfilter_status: PostfilterStatus::Pass,
        }
    }

    #[test]
    fn resolve_peptide_maps_pseudo_and_strips_nonalpha() {
        // No plan: pseudo residues stay, non-alpha dropped.
        assert_eq!(resolve_peptide("DT-HK", None), "DTHK");
    }

    #[test]
    fn locate_peptide_finds_1_based_start() {
        let proteins = protein_lookup(&fasta());
        // "DTHK" starts at residue 25 of TRYP_PIG (1-based).
        let pos = locate_peptide(&proteins, "P00761", "DTHK");
        assert_eq!(pos, Some(25));
    }

    #[test]
    fn parse_link_positions_splits_pair() {
        assert_eq!(parse_link_positions("3-1"), (Some(3), Some(1)));
        assert_eq!(parse_link_positions("5"), (Some(5), None));
    }

    #[test]
    fn xiview_csv_includes_absolute_positions() {
        let dir = std::env::temp_dir().join(format!("gq_xiview_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("xiview.csv");
        let mut hit = pass_hit("DTHK", "DTIVNELR", "P00761", "HRP", "2-3");
        hit.glyco_sites.push(super::super::filter::GlycoSite {
            peptide: 1,
            peptide_position: 3,
            residue: 'N',
            sequon_present: Some(false),
            plausible: false,
        });
        hit.n_glycan_pseudo = 2;
        hit.all_sites_plausible = false;
        let summary = write_xiview_csv(&path, &[hit], None, &fasta()).unwrap();
        assert_eq!(summary.rows, 1);
        assert_eq!(summary.mapped, 1);
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.starts_with("Protein1,PepPos1,PepSeq1,LinkPos1,AbsPos1"));
        // DTHK starts at 25 -> AbsPos1 = 25 + 2 - 1 = 26. DTIVNELR starts at 20 -> AbsPos2 = 20 + 3 - 1 = 22.
        assert!(content.contains("P00761,25,DTHK,2,26,HRP,20,DTIVNELR,3,22,"));
        assert!(content.contains("pep2:1:N:true:true;pep1:3:N:false:false"));
    }

    #[test]
    fn xiview_csv_skips_failed_hits() {
        let dir = std::env::temp_dir().join(format!("gq_xiview_fail_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("xiview.csv");
        let mut hit = pass_hit("DTHK", "DTIVNELR", "P00761", "HRP", "2-3");
        hit.postfilter_status = PostfilterStatus::Fail;
        let summary = write_xiview_csv(&path, &[hit], None, &fasta()).unwrap();
        assert_eq!(summary.rows, 0);
    }

    /// Render a populated report into `out/report_demo/` for manual/visual inspection.
    /// Ignored by default; run with `cargo test --lib demo_report -- --ignored`.
    #[test]
    #[ignore]
    fn demo_report() {
        use super::super::filter::HardStatus;
        let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("out/report_demo");
        std::fs::create_dir_all(&dir).unwrap();

        let glycans = [
            ("HexNAc(2)Hex(9)", 'N', 12.4, 1.2, true),
            ("HexNAc(2)Hex(9)", 'N', 10.1, -2.0, true),
            ("HexNAc(2)Hex(8)", 'N', 9.7, 3.1, true),
            ("HexNAc(4)Hex(5)NeuAc(2)", 'N', 8.9, 0.4, true),
            ("HexNAc(4)Hex(5)NeuAc(2)", 'N', 7.5, -1.1, false),
            ("HexNAc(1)", 'S', 6.8, 4.9, false),
            ("HexNAc(1)Hex(1)", 'T', 6.1, -3.4, false),
            ("HexNAc(2)Hex(5)", 'N', 5.5, 2.2, true),
        ];
        let mut hits: Vec<AnnotatedHit> = Vec::new();
        for (i, (glycan, residue, score, ppm, seq)) in glycans.iter().enumerate() {
            let mut h = pass_hit("DTHKSEIAHR", "DTIVNELR", "P00761", "HRP", "4-3");
            h.scan = Some(1000 + i as u32);
            h.glycan_composition = Some((*glycan).to_string());
            h.glycan_name = Some((*glycan).to_string());
            h.glyco_residue = Some(*residue);
            h.hit.score = *score;
            h.soft_score = *score + 0.8;
            h.hit.precursor_error_ppm = *ppm;
            h.sequon_present = Some(*seq);
            hits.push(h);
        }
        // A few filtered-out hits to populate the outcome breakdown.
        for status in [HardStatus::FailScore, HardStatus::FailPrecursorError] {
            let mut h = pass_hit("SEIAHR", "DTIVNELR", "P00761", "HRP", "2-3");
            h.hard_status = status;
            h.postfilter_status = PostfilterStatus::Fail;
            h.hit.score = 1.2;
            h.soft_score = 0.5;
            hits.push(h);
        }

        let ctx = ReportContext {
            project: "report_demo".into(),
            input_label: "published_gpx_reference.mzXML".into(),
            crosslinker_name: "dss".into(),
            xlink_sites: "K:K".into(),
            glycan_library: "nglyc309".into(),
        };
        let stats = PrefilterStats {
            scans_total: 24661,
            diagnostic_positive: 9952,
            isotope_pairs: 15,
            filtered_scans: 29,
            rejected: 14709,
        };
        write_html_report(&dir.join("report.html"), &hits, &stats, &ctx).unwrap();
        write_xiview_csv(&dir.join("xiview.csv"), &hits, None, &fasta()).unwrap();
        eprintln!("demo report written to {}", dir.display());
    }

    #[test]
    fn html_report_is_self_contained() {
        let dir = std::env::temp_dir().join(format!("gq_report_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("report.html");
        let ctx = ReportContext {
            project: "demo".into(),
            input_label: "run.mzXML".into(),
            crosslinker_name: "DSS".into(),
            xlink_sites: "K".into(),
            glycan_library: "n-glycans".into(),
        };
        let stats = PrefilterStats {
            scans_total: 100,
            diagnostic_positive: 40,
            isotope_pairs: 12,
            filtered_scans: 30,
            rejected: 60,
        };
        let hit = pass_hit("DTHK", "DTIVNELR", "P00761", "HRP", "2-3");
        write_html_report(&path, &[hit], &stats, &ctx).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.starts_with("<!DOCTYPE html>"));
        assert!(content.contains("Prefilter funnel"));
        assert!(content.contains("HexNAc(1)"));
        // Fully offline: no external script/stylesheet references.
        assert!(!content.contains("<script"));
        assert!(!content.contains("http://"));
    }
}
