# First run

This guide walks through a complete GlycoQuest workflow: validate configuration, execute the search, and open the interactive viewer.

## 1. Prepare inputs

You need:

1. **mzXML** — MS/MS data in xQuest-compatible mzXML format
2. **FASTA** — protein sequences (**IMPORTANT:** no literal `X`, `U`, `B`, or `J` in sequences; those letters are reserved for xQuest variable modifications)
3. **xQuest root** — path to `xQuest/V2.1.7/xquest` (or your installation)

Optional: copy `settings.ini` from the repository root to your working directory and edit tolerances, crosslinker chemistry, ppm tolerances, etc.

## 2. Dry-run (recommended)

Dry-run validates and writes the prefilter outputs and xQuest job folders **without** actually running xQuest:

```bash
./target/release/glycoquest tests/fixtures/mzxml/dss_pair.mzXML \
  --database data/target_proteins_asf.fasta \
  --xquest-root xQuest/V2.1.7/xquest \
  --out glycoquest_fixture_out \
  --dry-run
```

The two-scan fixture is a synthetic plumbing check, not a scientific benchmark. For a representative full test run, download the Xie et al. spectra from [MassIVE MSV000087442](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000087442) into `data/MSV000087442/` (not shipped in the repo). See [Published reference datasets](../theory/glycopeptide-crosslinking.md#published-reference-datasets) for citations and the search preset.

### What to check

**Terminal — readiness report**

```
readiness:
✓ PASS  MS input  ok
✓ PASS  FASTA database  ok
✓ PASS  Glycan library  ok
✓ PASS  xQuest runtime  ok
✓ PASS  Output directory  ok
overall: ✓ PASS  ready
```

**Terminal — prefilter summary**

```
prefilter: scans=2
prefilter: diagnostic_positive=2
prefilter: isotope_pairs=1
prefilter: filtered_scans=2
prefilter: rejected=0
```

**Terminal — job plan**

```
plan: 864 jobs, 3456 spectrum-job assignments, isotope_prefilter=true
dry-run: job folders and plan.json written (xQuest not executed)
```

**Output directory** (default: `out/<project>/` when `--out` is `out`):

| File / folder | Purpose |
|---------------|---------|
| `plan.json` | Full run plan: jobs, spectrum-job assignments, options |
| `filtered_spectra.tsv` | Spectra that passed all prefilters |
| `rejected_spectra.tsv` | Spectra rejected (with reason) |
| `isotope_pairs.tsv` | DSS light/heavy pairs (empty if isotope prefilter off) |
| `glycan_pruning.tsv` | Glycan candidates per retained spectrum |
| `spectra/` | Reduced mzXML (retained MS2 scans only) |
| `tmp/jobs/` | Generated xQuest job folders (see [xQuest jobs](../workflow/xquest-jobs.md)) |

If `filtered_scans=0`, the run exits with code **2**. Inspect `rejected_spectra.tsv` — see [Troubleshooting](troubleshooting.md).

## 3. Full run

Use your converted experimental data and matching FASTA, then omit `--dry-run` to execute xQuest on all generated jobs:

```bash
./target/release/glycoquest /path/to/experiment.mzXML \
  --database /path/to/proteins.fasta \
  --xquest-root xQuest/V2.1.7/xquest \
  --out glycoquest_out \
  --jobs 8
```

`--jobs 8` runs up to eight xQuest jobs in parallel (overrides `[execution] job_parallelism` in `settings.ini`). Use `0` or leave out for one worker per CPU core.

!!! warning Running unconstrained xQuest jobs
    Spawning many asynchronous Perl jobs will likely be heavy on your machine and may lead to crashes and/or data corruption. Setting `--jobs` to at most `max. avail. cores - 3` is advised.

### Expected output

```
run: executing N xQuest jobs across 8 worker thread(s)
run: wrote N hits from M result file(s) to results/glycoquest_xquest.csv
run: wrote results/xiview.csv (...), results/report.html, and results/viewer/ (...)
```

After success, `tmp/` is cleaned up, unless you specifically set `GLYCOQUEST_KEEP_TMP=1`.

## 4. Inspect results

### Quick QC — HTML report

Open in any browser (offline, no JavaScript required):

```text
out/<project>/results/report.html
```

See [HTML report](../results/report.md).

### Tabular hits — CSV

```text
out/<project>/results/glycoquest_xquest.csv
```

Tab-separated. See [Interpreting hits](../results/interpreting-hits.md).

### Interactive viewer

You can visualize the results locally using the built-in viewer:

```bash
cd out/<project>/results/viewer
python3 -m http.server 8080
```

Open [http://localhost:8080](http://localhost:8080). See [Using the viewer](../viewer/using-the-viewer.md).

## Exit codes

| Code | Meaning | Action |
|------|---------|--------|
| 0 | Success | No action required, output in `results/` |
| 1 | Validation / config error | Fix inputs or paths; check readiness report |
| 2 | No spectra passed prefilter | See `rejected_spectra.tsv` |
| 3 | xQuest job infrastructure failure | Check `tmp/logs/`, Perl, disk space |
| 4 | Result extraction failure | Check job XML in `tmp/jobs/` (use `GLYCOQUEST_KEEP_TMP=1`) |

Individual xQuest job failures are logged as warnings and listed in `results/failed_jobs.tsv`; they do **not** change GlycoQuest's exit code.

## Example: DMTMM (no isotope prefilter)

```bash
./target/release/glycoquest tests/fixtures/mzxml/hexnac_positive.mzXML \
  --database data/target_proteins_asf.fasta \
  --crosslinker dmtmm \
  --xquest-root xQuest/V2.1.7/xquest \
  --out glycoquest_dmtmm_out \
  --dry-run
```

See [Crosslinkers](../configuration/crosslinkers.md).

## Next steps

- [Workflow overview](../workflow/overview.md) — how each stage fits together
- [Settings reference](../configuration/settings-reference.md) — tune search parameters
