# xQuest jobs

After the prefilter, GlycoQuest materializes **one xQuest job directory per planned (spectrum × glycan variant)**. Jobs are transparent: you can inspect every file xQuest will use.

## Job directory layout

```text
out/<project>/tmp/jobs/<job_id>/
├── xquest.def              # Main search parameters
├── xmm.def                 # Mass mapping (if used)
├── glycoquest_matched.txt  # Spectrum matchlist for compare_peaks
├── database.fasta          # Protein database for this job
├── <pruned_stem>.mzXML     # Symlink to spectra/<file>
├── run.sh                  # Index → compare_peaks → xquest.pl
└── results/
    ├── xquest.xml          # Search results (on success)
    ├── results.spec.xml    # Intermediate compare_peaks output
    └── glycoquest_matched_isotopepairs.txt
```

After a successful full run, `tmp/` is deleted unless `GLYCOQUEST_KEEP_TMP=1`.

## Job ID and manifest

Each `job_id` encodes the glycan variant (composition, loss). GlycoQuest records a **job manifest** during generation mapping:

- `job_id` → glycan name, composition, mass, loss label
- `variable_mod` plan (pseudo-residue layout)
- Source spectrum keys

The manifest drives glycan annotation when consolidating `xquest.xml` hits.

## What `run.sh` does

1. Set `PERL5LIB` to bundled xQuest Perl modules under `--xquest-root`
2. Build a peptide database index in the job `results/` folder
3. Run `compare_peaks3.pl` to match MS/MS peaks against the matchlist
4. Run `xquest.pl` with `xquest.def`

Jobs do not share writable state — safe to run in parallel (`--jobs N`).

## `xquest.def` highlights

GlycoQuest writes DSS-oriented defaults (adjustable via settings and crosslinker preset):

| Parameter | DSS default | DMTMM (`label=none`) |
|-----------|-------------|----------------------|
| `xlinkermw` | 138.0680796 | −18.0109 |
| `isotopeshift` | 12.075321 | 0 |
| `AArequired` | `K:K` | `K:E,K:D` |
| `ms2_tolerance` | from `ms2_tolerance_da` | same |
| `variable_mod` | per-job glycan (+ optional oxidation) | per-job |

Dry-run tip: inspect `tmp/jobs/<first_job>/xquest.def` before committing to a long run.

## Matchlists

**DSS light/heavy:** matchlist rows pair light and heavy scans with m/z columns for isotope-aware compare_peaks.

**Unlabeled / light-only:** single-scan rows; `isotopeshift 0` in defs.

The matchlist connects prefilter-retained scans to the spectra xQuest actually searches.

## Variable modifications and pseudo-residues

xQuest rewrites modified residues to pseudo-letters:

| Pseudo | Typical GlycoQuest assignment |
|--------|------------------------------|
| `X` | First glycan target (e.g. glycosylated Asn) |
| `U` | Second target or Met oxidation |
| `B`, `J` | Additional mods if configured |

Example `variable_mod` for an N-glycan job with oxidation:

```text
variable_mod N,203.079373,M,15.994915
nvariable_mod 2
```

In xQuest output, `seq1`/`seq2` contain `X` where Asn carried the glycan. GlycoQuest decodes pseudos using the job manifest.

![Pseudo-residue decoding](../assets/pseudo-residue.svg)

**FASTA constraint:** Input sequences must not contain literal `X`, `U`, `B`, or `J`.

## One glycan variant per job

Each job searches **one** glycan mass (parent or water-loss form). Different compositions never share a job. This keeps:

- `nvariable_mod` within xQuest's four-pseudo limit
- Unambiguous glycan assignment in results

Trade-off: job count scales with retained spectra × pruned glycans × loss variants. Check `plan.json`:

```json
{
  "job_count": 576,
  "total_comparisons": 17496,
  "isotope_prefilter_enabled": true
}
```

## Resource limits

Before execution, GlycoQuest can abort if limits are exceeded (`[limits]` in settings):

- `max_jobs`
- `max_pruned_spectra`
- `max_total_job_spectrum_comparisons`

Dry-run always computes the plan without running xQuest.

## Completed-run artifacts

The viewer is written only after a completed run reaches consolidation. If a run is interrupted, start a new full run so prefilter, xQuest, postfilter, and Filtering data remain consistent.

## Related

- [Workflow overview](overview.md)
- [Post-filter](postfilter.md)
- [Crosslinkers](../configuration/crosslinkers.md)
- [Troubleshooting](../getting-started/troubleshooting.md)
