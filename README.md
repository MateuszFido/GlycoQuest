# GlycoQuest

GlycoQuest is a Rust CLI wrapper around xQuest for glycopeptide–peptide searches. It supports DSS, DMTMM, and published glycan-to-peptide crosslinkers; prefilters mzXML spectra using glycan diagnostic ions; optionally requires light/heavy isotope pairs; prunes glycan candidates; generates xQuest job folders; and consolidates results.

**Documentation:** full user guide in [`docs/`](docs/index.md) (installation, workflow, settings, results, viewer). A searchable site can be built with `pip install -r docs/requirements.txt && mkdocs serve`.

## Requirements

- Rust toolchain (edition 2024)
- xQuest V2.1.7 at `xQuest/V2.1.7/xquest` (or `--xquest-root`)
- Perl modules for xQuest: `DB_File` (Fedora: `dnf install perl-DB_File`), plus bundled libs under `xQuest/.../1209/` (set automatically in `run.sh` via `PERL5LIB`)
- mzXML MS/MS input (convert raw vendor files with msconvert first)

Default xQuest binary in [`settings.ini`](settings.ini): `xquest_bin = bin/xquest.pl` (resolved relative to `--xquest-root`).

## Quick start

### DSS (default, light/heavy isotope pairing)

```bash
cargo build --release
./target/release/glycoquest tests/fixtures/mzxml/dss_pair.mzXML \
  --database data/target_proteins_asf.fasta \
  --xquest-root /path/to/xQuest/V2.1.7/xquest \
  --out glycoquest_out \
  --dry-run
```

### DMTMM (zero-length, no isotope prefilter)

```bash
./target/release/glycoquest tests/fixtures/mzxml/hexnac_positive.mzXML \
  --database data/target_proteins_asf.fasta \
  --crosslinker dmtmm \
  --xquest-root /path/to/xQuest/V2.1.7/xquest \
  --out glycoquest_out \
  --dry-run
```

Advanced options (tolerances, crosslinker chemistry, job limits, modifications) live in [`settings.ini`](settings.ini). Copy it to your working directory or pass `--config`.

### Example

A bundled publicly available dataset (MassIVE DOI: 10.25345/C5VV5S, MSV000087442) is available for testing.
Citation: Xie, *et al.* Glycan–protein cross-linking mass spectrometry reveals sialic acid-mediated protein networks on cell surfaces. Chem. Sci. 2021; 12 (25): 8767–8777.

```bash
./target/release/glycoquest data/MSV000087442/PNT2-crosslink-in-situ.mzXML \
  --database data/MSV000087442/PNT2-GPx-focused-xquest.fasta \
  --glycans msv000087442-sianaz \
  --crosslinker nhs-cyclooctyne \
  --config configs/msv000087442-full.ini \
  --xquest-root V2.1.7/xquest \
  --out out
```

Its search encoding uses `X:K`: xQuest rewrites the glycan-bearing Asn to its first variable-modification pseudo-residue, `X`, so the glycan and crosslink occupy the same residue. The bond runs through SiaNAz on that N-glycan to a peptide Lys. The glycan library contains ordinary NeuAc, so the preset crosslink mass (205.085126607 Da) includes both the NHS–cyclooctyne residue and the SiaNAz-for-NeuAc mass difference.

### Parallel job execution

xQuest jobs are independent and run concurrently on a thread pool. Control the number of concurrent jobs with `--jobs`/`-j` (overrides `[execution] job_parallelism` in `settings.ini`); `0` means one worker per available CPU core.

Interactive runs show phase-aware progress and timing estimates for spectrum filtering, job preparation, parallel xQuest searches, and result consolidation. The display is enabled automatically for terminals; use `--progress never` to disable it or `--progress always` to force.

```bash
./target/release/glycoquest /path/to/experiment.mzXML \
  --database /path/to/proteins.fasta \
  --xquest-root xQuest/V2.1.7/xquest \
  --out glycoquest_out \
  --jobs 8
```

Each job builds its own database index in its result directory. Starting a new
search in an existing project output clears stopped-run job folders first, so
stale jobs cannot leak into the new output.


Glycan libraries via `--glycans`: a bundled id (`nglyc309` default, `oglyc78`,
`msv000087442-sianaz`) or a
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

Individual xQuest job failures are logged as warnings, listed in `results/failed_jobs.tsv`, and do not change the main exit code.

## Output layout (dry-run / run)

- `filtered_spectra.tsv`, `isotope_pairs.tsv`, `rejected_spectra.tsv`, `glycan_pruning.tsv` — prefilter results
- `spectra/` — reduced mzXML with only the retained MS2 scans (open in any mzXML viewer)
- `tmp/` — ephemeral xQuest job folders, logs, and peptide indexes (removed after a successful run)
- `plan.json` — normalized run plan and resource summary
- `results/glycoquest_xquest.csv` — consolidated, glycan-annotated, de-duplicated hits (run mode; header-only when no hits)
- `results/xiview.csv` — passing crosslinks in CLMS-CSV layout for network visualization tools
- `results/report.html` — self-contained QC/glycan report (prefilter funnel, score/error distributions, glycan breakdown, hits table)
- `results/viewer/` — interactive crosslink viewer (`index.html`, `viewer.json`, `database.fasta`); open `index.html` in a browser after a run
- `results/failed_jobs.tsv` — jobs where `run.sh` exited non-zero (any failure)

### Consolidated results columns

`results/glycoquest_xquest.csv` reports, per hit: `source_file`, `scan`,
`glycan_name`, `glycan_composition`, `glycan_mass`, `loss_label`, `glyco_residue`,
`glyco_peptide`, `n_glycan_pseudo`, `sequon_present`, `charge`, `charge_plausible`,
`matched_families`, `matched_ion_count`, `seq1`, `seq2`, `prot1`, `prot2`,
`topology`, `precursor_mz`, `mr`, `precursor_error_ppm`, `xlink_position`, `score`,
`hard_status`, `soft_score`, `postfilter_status`.

Post-filter hard requirements: configured crosslink evidence in the hit,
exactly one glycan pseudo-residue across the peptide
pair (V1 peptide-glycopeptide class), diagnostic-ion evidence in the originating
spectrum, precursor mass error within `max_precursor_error_ppm`, and `score` at or
above `min_score`. `soft_score` additionally rewards sequon presence, plausible
precursor charge, and diagnostic ion count.

### Visualization

The outputs of GlycoQuest can be visualized using the in-built viewer:

- `results/viewer/` — MIT-licensed interactive viewer (network, sequence map, MS/MS plot, Filtering, QC). Open `results/viewer/index.html` in a browser (offline; loads `viewer.json` and bundled FASTA).
- `results/xiview.csv` — passing crosslinks compatible with the CLMS-CSV layout (https://crosslinkviewer.org/index.php).
  (`Protein1,PepPos1,PepSeq1,LinkPos1,AbsPos1,Protein2,PepPos2,PepSeq2,LinkPos2,AbsPos2,Score,…`).
  Glycan composition, glycosylation residue, and loss label are appended as trailing columns.
- `results/report.html` — a dependency-free report (inline CSS + SVG charts, no network
  or JavaScript) with the prefilter funnel, post-filter outcome breakdown, xQuest score and
  precursor-error distributions, the glycan and glycosylation-site distribution over passing
  hits, and a hits table. Can be opened in any browser.

Default output layout: `out/<project>/` where `<project>` is derived from the first input mzXML filename. Override with `--out path/to/project_dir`.

See [`docs/index.md`](docs/index.md) for the complete user documentation and [`DESIGN.md`](DESIGN.md) for the V1 implementation specification.
