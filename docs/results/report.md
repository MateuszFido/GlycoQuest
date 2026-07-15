# HTML report

`results/report.html` is a **self-contained** quality-control and summary report. It requires no network connection and no JavaScript — open it directly in any browser from disk.

```bash
open out/<project>/results/report.html   # macOS
xdg-open out/<project>/results/report.html   # Linux
```

## Header metadata

The report title area shows:

- Project name (from output directory)
- Input mzXML label
- Crosslinker name and site rules
- Glycan library id

## Sections

### Prefilter funnel

Bar or funnel chart of:

1. Total MS2 scans parsed
2. Diagnostic-positive scans
3. Isotope pairs (if enabled)
4. Final filtered scans sent to xQuest

Compare with terminal `prefilter:` lines and `plan.json`.

### Post-filter outcomes

Breakdown of `hard_status` values across all hits — shows how many failed each rule vs `pass`.

Use this to decide whether to relax `min_score`, `max_precursor_error_ppm`, or prefilter tolerances.

### Score distribution

Histogram of xQuest `score` for hits. Passing hits are typically the right tail; adjust `[limits] min_score` based on your data.

### Precursor error distribution

Histogram of `precursor_error_ppm`. Tight clusters near zero indicate good mass calibration; wide tails may warrant stricter `max_precursor_error_ppm`.

### Glycan distribution

Counts of passing hits by `glycan_composition` — which structures were identified most often.

### Glycosylation site distribution

Counts by attachment residue (`N`, `S`, `T`) among passing hits.

### Hits table

Sortable summary of passing crosslinks with key columns (proteins, peptides, glycan, score, scan). For the full table, use `glycoquest_xquest.csv`.

## When to use report.html vs viewer

| Need | Use |
|------|-----|
| Quick QC PDF/screenshot | `report.html` |
| Offline share with no server | `report.html` |
| Interactive exploration, spectra, sequence map | [Viewer](../viewer/using-the-viewer.md) |
| Spreadsheet analysis | `glycoquest_xquest.csv` |
| CLMS-CSV network export | `xiview.csv` |

## Related

- [Output files](output-files.md)
- [Interpreting hits](interpreting-hits.md)
