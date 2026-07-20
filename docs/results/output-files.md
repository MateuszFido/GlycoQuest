# Output files

Project output lives under **`out/<project>/`** when using default `--out out`, or directly under your explicit `--out` path.

```text
out/<project>/
├── plan.json
├── filtered_spectra.tsv
├── isotope_pairs.tsv
├── rejected_spectra.tsv
├── glycan_pruning.tsv
├── input/                            # Staged FASTA (xQuest indexes here, not beside source)
│   └── <database>.fasta
├── spectra/                          # Reduced mzXML (retained MS2)
├── tmp/                              # Ephemeral (removed after success*)
│   ├── jobs/<job_id>/...
│   └── logs/<job_id>.log
└── results/
    ├── glycoquest_xquest.csv
    ├── xiview.csv
    ├── report.html
    ├── failed_jobs.tsv               # If any job failed
    └── viewer/
        ├── index.html
        ├── viewer.json
        ├── database.fasta
        ├── viewer.js
        └── viewer.css
```

\* Set `GLYCOQUEST_KEEP_TMP=1` to retain `tmp/`.

## Prefilter artifacts (project root)

| File | Description |
|------|-------------|
| `plan.json` | Normalized run plan: options, job list, comparison estimates, prefilter stats |
| `filtered_spectra.tsv` | Spectra passing all prefilters with diagnostic metadata |
| `isotope_pairs.tsv` | DSS light/heavy pairs (header-only if isotope prefilter disabled) |
| `rejected_spectra.tsv` | Rejected scans with `no_diagnostic` or `no_isotope_pair` |
| `glycan_pruning.tsv` | Glycan candidates retained per spectrum |
| `input/` | Copy of the `--database` FASTA; xQuest reads this instead of the source path |
| `spectra/` | mzXML containing only prefilter-retained MS2 scans |

### `filtered_spectra.tsv` columns

Key fields: source file, scan number, retention time, precursor m/z, charge, matched diagnostic families, matched ions.

### `rejected_spectra.tsv` columns

`source_file`, `scan_number`, `reason`.

## Ephemeral `tmp/`

| Path | Description |
|------|-------------|
| `tmp/jobs/<job_id>/` | Full xQuest job folder (see [xQuest jobs](../workflow/xquest-jobs.md)) |
| `tmp/logs/<job_id>.log` | Captured stdout/stderr from `run.sh` |

Useful for debugging failed jobs before cleanup.

## Results directory

| File | Description |
|------|-------------|
| `glycoquest_xquest.csv` | All consolidated hits (tab-separated), glycan-annotated, deduplicated |
| `xiview.csv` | Passing crosslinks in CLMS-CSV layout |
| `report.html` | Self-contained QC report (inline CSS/SVG, no network) |
| `failed_jobs.tsv` | Jobs whose `run.sh` exited non-zero |
| `viewer/` | Interactive crosslink viewer bundle |

## `plan.json` (selected fields)

| Field | Meaning |
|-------|---------|
| `job_count` | Number of xQuest jobs |
| `total_comparisons` | Exact spectrum–job assignments (the xQuest progress unit) |
| `isotope_prefilter_enabled` | Whether DSS pairing was required |
| `scans_total`, `filtered_scans` | Prefilter funnel |
| `jobs` | Per-job metadata (id, glycan, spectrum keys) |

## Empty results

| File | When empty |
|------|------------|
| `glycoquest_xquest.csv` | Header only — no hits or all post-filtered |
| `xiview.csv` | No passing mapped crosslinks |
| `viewer/viewer.json` | `crosslinks: []` — viewer still loads |
| `failed_jobs.tsv` | Absent or header-only if all jobs succeeded |

Exit code **0** even when individual jobs fail or hits are empty.

## Related

- [Interpreting hits](interpreting-hits.md)
- [HTML report](report.md)
- [Using the viewer](../viewer/using-the-viewer.md)
