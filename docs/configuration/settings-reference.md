# Settings reference

Advanced options live in **`settings.ini`**. Copy the bundled file from the repository root to your working directory, or pass `--config /path/to/settings.ini`.

Boolean values accept: `true`/`false`, `yes`/`no`, `1`/`0`, `on`/`off`.

## `[xquest]`

| Key | Default | Description |
|-----|---------|-------------|
| `xquest_bin` | `bin/xquest.pl` | Path to xQuest executable, relative to `--xquest-root` unless absolute |

## `[tolerances]`

| Key | Default | xQuest / GlycoQuest use |
|-----|---------|-------------------------|
| `diagnostic_tolerance_ppm` | `10` | GlycoQuest: oxonium-ion matching in prefilter. Overridden by CLI `--ppm-tolerance` |
| `ms1_tolerance_ppm` | `10` | MS1 precursor tolerance (general) |
| `ms2_tolerance_da` | `0.2` | Fragment ion tolerance written to `xquest.def` |
| `isotope_pair_ms1_tolerance_ppm` | `10` | DSS light/heavy pair matching in MS1 |
| `isotope_pair_rt_tolerance_min` | `2.0` | Retention-time window for isotope pairs (minutes) |
| `neutral_loss_tolerance_da` | `0.05` | Reserved for neutral-loss annotation (catalog expansion) |

**When to change:** Increase diagnostic ppm if valid glycopeptide spectra are rejected (`no_diagnostic`). Increase isotope RT window if light/heavy pairs are missed on long gradients.

## `[crosslinker]`

| Key | Default | Description |
|-----|---------|-------------|
| `name` | `dss` | Bundled preset (`dss`, `dmtmm`, `nhs-cyclooctyne`, `ssbxl`, `pcbxl`) or custom chemistry name. Overridden by CLI `--crosslinker` |
| `label` | `light-heavy` | `light-heavy`, `light-only`, or `none` — controls isotope prefilter and xQuest isotope settings |
| `shift_da` | `12.075321` | Heavy-label mass shift for DSS isotope coding |
| `xlinkermw` | `138.0680796` | Crosslinker monoisotopic mass in xQuest (Da) |
| `xlink_sites` | `K:K` | Allowed site pairs (`AArequired` syntax), e.g. `K:E,K:D` for DMTMM or `X:K` for GPx (`X` = glycan-modified Asn) |
| `nterm_xlinkable` | `false` | Enable N-terminal crosslinking (`ntermxlinkable 1`, `Z` pseudo-site) |

See [Crosslinkers](crosslinkers.md) for presets and labeling modes.

## `[modifications]`

| Key | Default | xQuest mapping |
|-----|---------|----------------|
| `fixed_carbamidomethyl_cys` | `true` | Fixed carbamidomethylation on Cys |
| `variable_oxidation` | `true` | Variable Met oxidation (`M,15.994915`) as second pseudo-residue when enabled |

Glycan masses are added as **additional** variable modifications per job (not in this section).

## `[glycan]`

| Key | Default | Description |
|-----|---------|-------------|
| `targets` | `N,S,T` | Residues considered for glycan attachment in library validation |

The active library is selected by CLI `--glycans`. See [Glycan libraries](glycan-libraries.md).

## `[limits]`

Resource guards and result post-filters. `0` = no limit.

| Key | Default | Description |
|-----|---------|-------------|
| `max_jobs` | `0` | Abort before xQuest if job count exceeds limit |
| `max_pruned_spectra` | `0` | Cap spectra after glycan pruning |
| `max_total_job_spectrum_comparisons` | `0` | Cap estimated xQuest comparisons (dry-run / plan) |
| `min_score` | `0` | Minimum xQuest score for `postfilter_status=pass` |
| `max_precursor_error_ppm` | `20` | Maximum \|precursor error\| (ppm) for passing hits |

## `[execution]`

| Key | Default | Description |
|-----|---------|-------------|
| `job_parallelism` | `0` | Concurrent xQuest jobs; `0` = one thread per CPU core. Overridden by CLI `--jobs` |

## Example: stricter post-filtering

```ini
[limits]
min_score = 10
max_precursor_error_ppm = 10
```

## Example: DMTMM in settings (alternative to `--crosslinker dmtmm`)

```ini
[crosslinker]
name = dmtmm
label = none
shift_da = 0
xlinkermw = -18.0109
xlink_sites = K:E,K:D
```

## Missing settings file

If `settings.ini` is not found, GlycoQuest prints a warning and uses built-in defaults (same values as the bundled file).
