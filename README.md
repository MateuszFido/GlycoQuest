# GlycoQuest

GlycoQuest is a Rust CLI wrapper around xQuest for DSS- and DMTMM-crosslinked glycopeptide–peptide searches. It prefilters mzXML spectra using glycan diagnostic ions, optionally requires DSS light/heavy isotope pairs, prunes glycan candidates, generates xQuest job folders, and consolidates results.

## Requirements

- Rust toolchain (edition 2024)
- xQuest V2.1.6 at `xQuest/V2.1.6/xquest` (or `--xquest-root`)
- Perl modules for xQuest: `DB_File` (Fedora: `dnf install perl-DB_File`), plus bundled libs under `xQuest/.../1209/` (set automatically in `run.sh` via `PERL5LIB`)
- mzXML MS/MS input (convert raw vendor files with msconvert first)

Default xQuest binary in [`settings.ini`](settings.ini): `xquest_bin = bin/xquest.pl` (resolved relative to `--xquest-root`).

## Quick start

### DSS (default, light/heavy isotope pairing)

```bash
cargo build --release
./target/release/glycoquest tests/fixtures/mzxml/dss_pair.mzXML \
  --database data/rcsb_pdb_1HRP_no_contams.fasta \
  --xquest-root /path/to/xQuest/V2.1.6/xquest \
  --out glycoquest_out \
  --dry-run
```

### DMTMM (zero-length, no isotope prefilter)

```bash
./target/release/glycoquest tests/fixtures/mzxml/hexnac_positive.mzXML \
  --database data/rcsb_pdb_1HRP_no_contams.fasta \
  --crosslinker dmtmm \
  --xquest-root /path/to/xQuest/V2.1.6/xquest \
  --out glycoquest_out \
  --dry-run
```

Advanced options (tolerances, crosslinker chemistry, job limits, modifications) live in [`settings.ini`](settings.ini). Copy it to your working directory or pass `--config`.

### Real hCG validation (local data)

Dry-run on the msconvert hCG file (~30k MS2 scans, ~80s parse on a typical workstation):

```bash
./target/release/glycoquest data/260521_LU02_disoic_hCG_01.mzXML \
  --database data/rcsb_pdb_1HRP_no_contams.fasta \
  --xquest-root xQuest/V2.1.6/xquest \
  --out glycoquest_hcg_out \
  --dry-run
```

Observed dry-run stats (Jun 2026): 24 250 MS2 scans, 8 436 diagnostic-positive, 9 DSS isotope pairs → 18 filtered spectra, 576 jobs, 17 496 spectrum comparisons. See `plan.json` for the full run plan.

Full run (executes `compare_peaks3.pl` + `xquest.pl` per job under `jobs/`):

```bash
./target/release/glycoquest data/260521_LU02_disoic_hCG_01.mzXML \
  --database data/rcsb_pdb_1HRP_no_contams.fasta \
  --xquest-root xQuest/V2.1.6/xquest \
  --out glycoquest_hcg_out
```

Requires `perl-DB_File` for the xQuest search step. Each job writes `results/results.spec.xml`, `glycoquest_matched_isotopepairs.txt`, and on success `results/xquest.xml`; consolidated hits land in `results/glycoquest_xquest.csv`.

### Parallel job execution

xQuest jobs are independent and run concurrently on a thread pool. Control the number of concurrent jobs with `--jobs`/`-j` (overrides `[execution] job_parallelism` in `settings.ini`); `0` means one worker per available CPU core.

```bash
./target/release/glycoquest data/260521_LU02_disoic_hCG_01.mzXML \
  --database data/rcsb_pdb_1HRP_no_contams.fasta \
  --xquest-root xQuest/V2.1.6/xquest \
  --out glycoquest_hcg_out \
  --jobs 8
```

Each job builds its own database index in its result directory, so concurrent jobs do not share writable state. Concurrency is also honored by `--resume`.

Optional integration tests (skipped unless env vars are set):

```bash
export GLYCOQUEST_HCG_MZXML=data/260521_LU02_disoic_hCG_01.mzXML   # mzXML parser smoke test
export GLYCOQUEST_XQUEST_ROOT=xQuest/V2.1.6/xquest                  # reserved for future xQuest CI
cargo test
```

Glycan libraries via `--glycans`: a bundled id (`nglyc309` default, `oglyc78`) or a
path to a custom CSV/TSV file with columns
`name,composition,monoisotopic_mass,diagnostic_ions,residue_targets`
(`diagnostic_ions` is a `;`-separated list of `family@mz` with optional `[-loss]`;
`residue_targets` is a `;`-separated list of attachment residues, e.g. `N` or `S;T`).
O-glycan libraries drive the xQuest `variable_mod` onto Ser/Thr; N-glycans onto Asn.

## Exit codes

| Code | Meaning |
|------|---------|
| 0 | Success (results written; may be empty if all xQuest jobs failed) |
| 1 | Validation or configuration error |
| 2 | No spectra passed filters |
| 3 | xQuest job infrastructure failure (cannot read jobs directory or create logs) |
| 4 | Result extraction failure |

Individual xQuest job failures are logged as warnings, listed in `results/failed_jobs.tsv`, and do not change the exit code.

## Output layout (dry-run / run)

- `filtered_spectra.tsv`, `isotope_pairs.tsv`, `rejected_spectra.tsv`, `glycan_pruning.tsv` — prefilter results
- `spectra/` — reduced mzXML with only prefilter-retained MS2 scans (open in any mzXML viewer)
- `tmp/` — ephemeral xQuest job folders, logs, and peptide indexes (removed after a successful run)
- `plan.json` — normalized run plan and resource summary
- `results/glycoquest_xquest.csv` — consolidated, glycan-annotated, de-duplicated hits (run mode; header-only when no hits)
- `results/xiview.csv` — passing crosslinks in xiVIEW CSV layout for the interactive network view
- `results/report.html` — self-contained QC/glycan report (prefilter funnel, score/error distributions, glycan breakdown, hits table)
- `results/viewer/` — interactive crosslink viewer (`index.html`, `viewer.json`, `database.fasta`); open `index.html` in a browser after a run
- `results/failed_jobs.tsv` — jobs whose `run.sh` exited non-zero (when any fail)

### Consolidated results columns

`results/glycoquest_xquest.csv` reports, per hit: `source_file`, `scan`,
`glycan_name`, `glycan_composition`, `glycan_mass`, `loss_label`, `glyco_residue`,
`glyco_peptide`, `n_glycan_pseudo`, `sequon_present`, `charge`, `charge_plausible`,
`matched_families`, `matched_ion_count`, `seq1`, `seq2`, `prot1`, `prot2`,
`topology`, `precursor_mz`, `mr`, `precursor_error_ppm`, `xlink_position`, `score`,
`hard_status`, `soft_score`, `postfilter_status`.

Post-filter hard requirements (all must hold for `postfilter_status=pass`): DSS
crosslink evidence in the hit, exactly one glycan pseudo-residue across the peptide
pair (V1 peptide-glycopeptide class), diagnostic-ion evidence in the originating
spectrum, precursor mass error within `max_precursor_error_ppm`, and `score` at or
above `min_score`. `soft_score` additionally rewards sequon presence, plausible
precursor charge, and diagnostic-ion count.

Note: `--resume` re-consolidates existing job results without prefilter/manifest
state, so it produces reduced annotation (glycan label reconstructed from the job id;
diagnostic-linkage and sequon soft features are skipped). Run a full search for
fully annotated results.

### Visualization outputs

Three artifacts support scientific analysis of crosslinks on the actual proteins:

- `results/viewer/` — MIT-licensed interactive viewer (network, sequence map, MS/MS mirror plot, QC). Open `results/viewer/index.html` in a browser (offline; loads `viewer.json` and bundled FASTA). Build or refresh static assets with `cd viewer && npm install && npm run build` before running GlycoQuest if you change the viewer sources.
- `results/xiview.csv` — passing crosslinks in the [xiVIEW CSV layout](https://www.xiview.org/csv-formats.php)
  (`Protein1,PepPos1,PepSeq1,LinkPos1,AbsPos1,Protein2,PepPos2,PepSeq2,LinkPos2,AbsPos2,Score,…`).
  Peptides are resolved back from xQuest pseudo-residues and located in the FASTA to
  yield 1-based peptide start positions and absolute cross-link residue numbers, so the
  file loads directly at [xiview.org](https://xiview.org) for the 2D network view. Glycan
  composition, glycosylation residue, and loss label are appended as trailing columns
  (xiVIEW ignores unrecognized columns).
- `results/report.html` — a dependency-free report (inline CSS + SVG charts, no network
  or JavaScript) with the prefilter funnel, post-filter outcome breakdown, xQuest score and
  precursor-error distributions, the glycan and glycosylation-site distribution over passing
  hits, and a hits table. Open it directly in any browser.

Default output layout: `out/<project>/` where `<project>` is derived from the first input mzXML filename. Override with `--out path/to/project_dir`.

See [`DESIGN.md`](DESIGN.md) for the full V1 specification.
