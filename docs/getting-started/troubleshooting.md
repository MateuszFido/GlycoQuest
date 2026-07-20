# Troubleshooting

Decision guide for common GlycoQuest failures and empty results.

## Readiness check failures

GlycoQuest prints a **readiness** report before running. Any `✗ FAILED` line stops execution (exit code **1**).

### MS input failed

- Confirm the path exists and is `.mzXML` (case-sensitive on Linux).
- Directories must contain at least one mzXML file.
- Vendor raw formats (`.raw`, `.wiff`, etc.) are rejected — convert with msconvert.

### FASTA database failed

- File must exist and contain at least one `>` header entry.
- Sequences must not contain `X`, `U`, `B`, or `J` (reserved for xQuest variable modifications).

### Glycan library failed

- Bundled ids: `nglyc309`, `oglyc78`, and `msv000087442-sianaz` (plus aliases).
- Custom files must be CSV/TSV with columns: `name`, `composition`, `monoisotopic_mass`, `diagnostic_ions`, `residue_targets`.
- See [Glycan libraries](../configuration/glycan-libraries.md).

### xQuest runtime failed

- `--xquest-root` must point to a directory containing the executable from `settings.ini` (`bin/xquest.pl` by default).

### xQuest Perl modules failed

GlycoQuest readiness and `scripts/run.sh` compile-check the search path (`compare_peaks3.pl`, `xquest.pl`) before launching jobs. Typical gaps on HPC (Euler): **`DB_File`** (`libdb`) and **`XML::Parser`** (`libexpat`). Most pure-Perl libs ship under `V2.1.7/xquest/1209/{lib,share}/perl5`.

**Do not** add `1209/lib64` to `PERL5LIB` — those `.so` files target an old Perl and crash 5.38 (`Perl_Gthr_key_ptr`). GD is optional for GlycoQuest (`drawspectra 0`).

```bash
# ETH Euler (once, login node):
unset PERL5LIB   # if you previously exported lib64 paths
scripts/bootstrap-euler-perl.sh
scripts/check-xquest-perl.pl

# Fedora/RHEL:  sudo dnf install perl-DB_File perl-XML-Parser
# Debian/Ubuntu: sudo apt install libdb-file-perl libxml-parser-perl
```

If a module is installed but `.so` load fails (`libdb` / `libexpat`), the Slurm wrapper adds matching Spack lib dirs to `LD_LIBRARY_PATH`.

### Output directory failed

- Parent path must be writable.
- Explicit `--out` is used as the project root; default `out` creates `out/<project_slug>/`.

## Exit code 2 — no spectra retained

All MS2 scans were rejected by the prefilter.

1. Open `rejected_spectra.tsv` and count reasons:

   | `reason` | Meaning |
   |----------|---------|
   | `no_diagnostic` | No glycan diagnostic ion matched within ppm tolerance |
   | `no_isotope_pair` | Diagnostic-positive but no DSS light/heavy partner (DSS only) |

2. **Many `no_diagnostic`**
   - Loosen `--ppm-tolerance` or `[tolerances] diagnostic_tolerance_ppm`.
   - Confirm the sample contains glycopeptides and the glycan library matches both
     glycosylation type and crosslinker chemistry.
   - Inspect a few MS/MS scans in an external viewer for oxonium ions (e.g. HexNAc ~204 m/z).

3. **Many `no_isotope_pair`**
   - Confirm the experiment used DSS light/heavy duplex labeling.
   - Widen `[tolerances] isotope_pair_ms1_tolerance_ppm` or `isotope_pair_rt_tolerance_min`.
   - For unlabeled crosslinkers, use `--crosslinker dmtmm` or `[crosslinker] label = none`.

4. Compare `prefilter: scans=` vs `diagnostic_positive=` in terminal output.

## Empty or missing hits (exit code 0)

GlycoQuest can succeed with zero rows in `glycoquest_xquest.csv`.

| Situation | What to check |
|-----------|---------------|
| All jobs failed | `results/failed_jobs.tsv`, logs under `tmp/logs/` (set `GLYCOQUEST_KEEP_TMP=1`) |
| Jobs succeeded, no XML hits | Normal for stringent data; lower `min_score` or review prefilter retention |
| Hits in XML but empty CSV | Post-filter removed all hits; open viewer with **Show failed hits** or inspect `hard_status` column |

### Debugging failed xQuest jobs

```bash
export GLYCOQUEST_KEEP_TMP=1
# re-run a full search
less out/<project>/tmp/logs/<job_id>.log
```

Common log issues: Perl module missing, FASTA indexing error, corrupt mzXML symlink.

## Exit code 3 — xQuest infrastructure failure

Cannot create or read the jobs directory, or thread pool failure. Usually a filesystem permissions problem.

## Exit code 4 — result extraction failure

xQuest produced XML that GlycoQuest could not parse, or viewer asset install failed. Keep `tmp/` and inspect `tmp/jobs/*/results/xquest.xml`.

## Viewer issues

### "Failed to load viewer" / blank page

**Cause:** Opening `index.html` via `file://`. Browsers block loading `viewer.json`.

**Fix:** Serve over HTTP:

```bash
cd out/<project>/results/viewer
python3 -m http.server 8080
# or: ./serve-viewer.sh 8080
```

### "No mapped crosslinks to display"

FASTA protein IDs do not match xQuest `prot1`/`prot2` strings. Ensure FASTA headers are consistent with what xQuest reports. Hits may still appear in the table with `mapped: false`.

### Spectrum panel empty

Reduced spectrum not found in `spectra/` for the hit's scan. This usually means scan parsing failed or the result directory is incomplete.

## Performance

| Symptom | Mitigation |
|---------|------------|
| Too many jobs / long runtime | Tighten prefilter; set `[limits] max_jobs` or `max_pruned_spectra` |
| Disk full | Each job builds a peptide index; use limits or smaller FASTA |
| Slow mzXML parse | Expected for large files (~30k scans); dry-run estimates workload via `plan.json` |

## Getting help

When reporting issues, include:

- GlycoQuest version (`glycoquest --version`)
- Command line (minus private paths if needed)
- Readiness report and prefilter summary
- First 20 lines of `rejected_spectra.tsv` or `failed_jobs.tsv`
- `plan.json` job count and `isotope_prefilter_enabled`
