# CLI reference

```bash
glycoquest [OPTIONS] <INPUT>
```

## Positional arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `INPUT` | Yes | mzXML file or directory of mzXML files (xQuest-compatible) |

## Flags and options

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--database` | | *(required)* | Protein sequence database (FASTA) |
| `--glycans` | | `nglyc309` | Bundled glycan id (`nglyc309`, `oglyc78`, `msv000087442-sianaz`) or path to custom CSV/TSV |
| `--xquest-root` | | `.` | xQuest installation root (contains `bin/xquest.pl`) |
| `--crosslinker` | | from settings | Crosslinker preset name (`dss`, `dmtmm`, `nhs-cyclooctyne`, `ssbxl`, `pcbxl`) |
| `--ppm-tolerance` | | from settings | Diagnostic-ion matching tolerance (ppm); overrides `settings.ini` |
| `--jobs` | `-j` | from settings | Concurrent xQuest jobs; `0` = one per CPU core |
| `--out` | | `out` | Output base; default creates `out/<project>/` from input filename |
| `--config` | | `settings.ini` | Path to settings file |
| `--progress` | | `auto` | Live progress: `auto`, `always`, or `never` |
| `--dry-run` | | off | Validate and write job plan without executing xQuest |
| `--help` | `-h` | | Print help |
| `--version` | `-V` | | Print version |

## Execution modes

- Default mode performs a full completed run: prefilter, xQuest execution, consolidation, report, and viewer.
- `--dry-run` validates inputs and writes the job plan without executing xQuest.

## Live progress

The default `--progress auto` renders phase-aware progress bars and ETAs when stderr is an interactive terminal. It tracks Rust spectrum filtering and job generation, then aggregates the progress of concurrently running xQuest jobs. Redirected output and CI logs retain plain milestone messages without terminal control sequences.

Use `--progress never` to disable the display or `--progress always` to force it when running through a terminal wrapper.

## CLI vs settings.ini

CLI flags override `settings.ini` for:

| CLI | Settings key |
|-----|----------------|
| `--crosslinker` | `[crosslinker] name` (preset fields still from ini unless preset replaces them) |
| `--ppm-tolerance` | `[tolerances] diagnostic_tolerance_ppm` |
| `--jobs` | `[execution] job_parallelism` |

All other advanced options live only in `settings.ini`. See [Settings reference](settings-reference.md).

## Output directory resolution

| `--out` value | Project root |
|---------------|--------------|
| `out` (default) | `out/<slug>/` where slug is derived from the first mzXML filename |
| Any other path | Used directly as the project directory |

Example: input `260607_LU02_disoic_ASF_DSS_1.c.mzXML` → `out/260607_lu02_disoic_asf_dss_1_c/`.

## Examples

**Minimal dry-run**

```bash
glycoquest run.mzXML --database proteins.fasta --dry-run
```

**DSS search with explicit xQuest path and parallelism**

```bash
glycoquest data.mzXML \
  --database proteins.fasta \
  --xquest-root /opt/xQuest/V2.1.7/xquest \
  --out my_project \
  --jobs 8
```

**O-glycan library**

```bash
glycoquest run.mzXML --database proteins.fasta --glycans oglyc78
```

**Published NHS–cyclooctyne GPx chemistry**

```bash
glycoquest run.mzXML \
  --database proteins.fasta \
  --glycans msv000087442-sianaz \
  --crosslinker nhs-cyclooctyne
```

**Custom settings file**

```bash
glycoquest run.mzXML --database proteins.fasta --config /path/to/my_settings.ini
```

## Exit codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Validation or configuration error |
| 2 | No spectra passed prefilter |
| 3 | xQuest job infrastructure failure |
| 4 | Result extraction failure |

See [Troubleshooting](../getting-started/troubleshooting.md).
