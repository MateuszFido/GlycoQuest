# GlycoQuest documentation

**GlycoQuest** finds **glycopeptide–peptide** pairs in LC–MS/MS data using DSS, DMTMM, or published glycan-to-peptide linker models. It prefilters spectra on glycan diagnostic ions, detects light/heavy isotope pairs for labelled crosslinkers, searches retained spectra with [xQuest](https://www.nature.com/articles/nprot.2013.168), and returns annotated hits auditable in the browser.

## Start here

| I want to… | Read |
|------------|------|
| Install and run my first search | [Installation](getting-started/installation.md) → [First run](getting-started/first-run.md) |
| Understand the full pipeline | [Workflow overview](workflow/overview.md) |
| Configure tolerances, crosslinker, limits | [Settings reference](configuration/settings-reference.md) |
| Interpret CSV / report / viewer output | [Interpreting hits](results/interpreting-hits.md) |
| Open the interactive viewer | [Using the viewer](viewer/using-the-viewer.md) |
| Fix a failed or empty run | [Troubleshooting](getting-started/troubleshooting.md) |

## Documentation map

### Getting started

- [Installation](getting-started/installation.md) — pre-built binaries (Windows, macOS, Linux), xQuest, build from source
- [First run](getting-started/first-run.md) — dry-run → full run → viewer
- [Troubleshooting](getting-started/troubleshooting.md) — exit codes, common failures

### Configuration

- [CLI reference](configuration/cli-reference.md) — command-line flags
- [Settings reference](configuration/settings-reference.md) — `settings.ini` sections
- [Crosslinkers](configuration/crosslinkers.md) — DSS, DMTMM, NHS–cyclooctyne, SSBXL, PCBXL, and labeling modes
- [Glycan libraries](configuration/glycan-libraries.md) — bundled and custom libraries

### Workflow

- [Overview](workflow/overview.md) — end-to-end pipeline with diagram
- [Prefilter](workflow/prefilter.md) — diagnostic ions, isotope pairs, glycan pruning
- [xQuest jobs](workflow/xquest-jobs.md) — job folders, defs, matchlists, execution
- [Post-filter](workflow/postfilter.md) — hard/soft rules and scoring

### Results

- [Output files](results/output-files.md) — every artifact in `out/<project>/`
- [Interpreting hits](results/interpreting-hits.md) — scientific meaning of CSV columns
- [HTML report](results/report.md) — `report.html` walkthrough
- [Network CSV export](results/network-csv-export.md) — `xiview.csv` (CLMS-CSV layout)

### Viewer

- [Using the viewer](viewer/using-the-viewer.md) — launch, panels, filters
- [Scientific interpretation](viewer/scientific-interpretation.md) — what the viewer shows and does not claim

### Background

- [Glycopeptide crosslinking](science/glycopeptide-crosslinking.md) — biology and MS context
- [Diagnostic ions](science/diagnostic-ions.md) — oxonium ions and neutral losses

## Building this documentation site

The docs are written in Markdown and can be read directly in the repository. To build a searchable HTML site with [MkDocs Material](https://squidfunk.github.io/mkdocs-material/):

```bash
pip install -r docs/requirements.txt
mkdocs serve    # http://127.0.0.1:8000
mkdocs build    # site/ output
```

See `mkdocs.yml` at the repository root for site configuration.

## Related documents

- `README.md` at the repository root — project overview and quick start
- `DESIGN.md` — V1 implementation specification (contributors)
