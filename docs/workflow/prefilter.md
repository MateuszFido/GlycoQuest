# Prefilter

The prefilter reduces full mzXML files to a small set of **glycan-supported** MS/MS spectra before xQuest runs. Three gates run in sequence.

![Isotope pairing](../assets/isotope-pair.svg)

## Gate 1: Diagnostic ion matching

**Question:** Does this MS/MS scan contain evidence of a glycan?

For each MS2 scan, GlycoQuest matches peaks against the union of diagnostic targets from the glycan library, within `diagnostic_tolerance_ppm` (default 10 ppm, overridable with `--ppm-tolerance`).

| Behavior | Detail |
|----------|--------|
| Match criterion | Any observed peak within ppm of an expected oxonium ion |
| Intensity threshold | None — any intensity counts |
| Pass condition | At least one diagnostic **family** matched (e.g. `HexNAc`, `NeuAc`) |
| Fail reason | `no_diagnostic` in `rejected_spectra.tsv` |

### Example diagnostic ions (N-glycans)

| Family | Typical m/z | Notes |
|--------|-------------|-------|
| HexNAc | ~204.09 | Core N-acetylhexosamine |
| HexNAc | ~186.08 | −H₂O neutral loss |
| NeuAc | ~292.10 | Sialic acid |

See [Diagnostic ions](../science/diagnostic-ions.md) for the full catalog.

**Output:** `matched_families` and `matched_ions` recorded in `filtered_spectra.tsv`.

## Gate 2: DSS isotope pairing

**Question:** For DSS duplex data, is there a light/heavy precursor partner?

**Enabled when:** `[crosslinker] label = light-heavy` (DSS default).

| Parameter | Default | Role |
|-----------|---------|------|
| `shift_da` | 12.075321 | Expected Δm/z between partners |
| `isotope_pair_ms1_tolerance_ppm` | 10 | MS1 m/z tolerance |
| `isotope_pair_rt_tolerance_min` | 2.0 | Retention time window |

A diagnostic-positive scan must have a partner scan whose precursor m/z differs by the shift (within tolerances) and whose retention time is within the RT window.

| Mode | Isotope gate |
|------|--------------|
| `light-heavy` | Required |
| `light-only`, `none` | Skipped — all diagnostic-positive scans proceed |
| Fail reason | `no_isotope_pair` |

**Output:** `isotope_pairs.tsv` lists light/heavy scan pairs. When disabled, the file has headers only.

## Gate 3: Glycan pruning

**Question:** Which glycans from the library are consistent with this spectrum's diagnostic evidence?

Pruning is **per spectrum**, based on diagnostic **families** (not individual ion m/z):

| Spectrum evidence | Glycans retained |
|-------------------|------------------|
| HexNAc only | All library entries containing HexNAc |
| NeuAc only | All containing NeuAc |
| HexNAc + NeuAc | All containing **both** |

**Output:** `glycan_pruning.tsv` — one row per (scan, glycan candidate).

## Reduced mzXML

Retained scans are written to `spectra/` as xQuest-sized mzXML subsets. These files:

- Feed xQuest job symlinks
- Supply MS/MS peaks for the interactive viewer mirror plots
- Can be opened in any mzXML viewer for manual inspection

## Prefilter statistics

Terminal summary:

```
prefilter: scans=24250
prefilter: diagnostic_positive=8436
prefilter: isotope_pairs=9
prefilter: filtered_scans=18
prefilter: rejected=24232
```

The HTML report and viewer QC panel show the same funnel as a chart.

## Tuning the prefilter

| Symptom | Adjustment |
|---------|------------|
| Everything rejected (`no_diagnostic`) | Increase `diagnostic_tolerance_ppm`; verify glycan library type |
| Diagnostic positives lost at isotope step | Widen MS1/RT tolerances; confirm DSS labeling |
| Too many spectra / jobs | Stricter ppm; add `[limits] max_pruned_spectra` |
| DMTMM / unlabeled data | Use `--crosslinker dmtmm` or `label = none` |

## Related

- [Workflow overview](overview.md)
- [Glycan libraries](../configuration/glycan-libraries.md)
- [Crosslinkers](../configuration/crosslinkers.md)
