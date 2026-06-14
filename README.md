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

Optional integration tests (skipped unless env vars are set):

```bash
export GLYCOQUEST_HCG_MZXML=data/260521_LU02_disoic_hCG_01.mzXML   # mzXML parser smoke test
export GLYCOQUEST_XQUEST_ROOT=xQuest/V2.1.6/xquest                  # reserved for future xQuest CI
cargo test
```

Bundled glycan libraries: `nglyc309` (default), `oglyc78` via `--glycans`.

## Exit codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Validation or configuration error |
| 2 | No spectra passed filters |
| 3 | xQuest execution failure |
| 4 | Result extraction failure |

## Output layout (dry-run / run)

- `filtered_spectra.tsv`, `isotope_pairs.tsv`, `rejected_spectra.tsv`, `glycan_pruning.tsv` — prefilter results
- `spectra/` — reduced mzXML with only prefilter-retained MS2 scans (open in any mzXML viewer)
- `tmp/` — ephemeral xQuest job folders, logs, and peptide indexes (removed after a successful run)
- `plan.json` — normalized run plan and resource summary
- `results/glycoquest_xquest.csv` — consolidated hits (run mode)

Default output layout: `out/<project>/` where `<project>` is derived from the first input mzXML filename. Override with `--out path/to/project_dir`.

See [`DESIGN.md`](DESIGN.md) for the full V1 specification.
