# GlycoQuest documentation

**GlycoQuest** finds **glycopeptide–peptide** crosslinks in LC–MS/MS data using published glycopeptide-to-peptide and glycan-to-peptide linker models. It prefilters spectra on glycan diagnostic ions (e.g., Sia, Kdo or HexNAc), detects light/heavy isotope pairs for labelled crosslinkers, searches retained spectra with [xQuest](https://www.nature.com/articles/nprot.2013.168), and returns annotated hits viewable in the browser.

## Start here

| I want to… | Read |
|------------|------|
| Install and run | [Installation](getting-started/installation.md) → [First run](getting-started/first-run.md) |
| Understand the pipeline | [Workflow overview](workflow/overview.md) |
| Configure tolerances, crosslinker properties | [Settings reference](configuration/settings-reference.md) |
| Interpret CSV / report / viewer output | [Interpreting hits](results/interpreting-hits.md) |
| Open the interactive viewer | [Using the viewer](viewer/using-the-viewer.md) |
| Fix a failed or empty run | [Troubleshooting](getting-started/troubleshooting.md) |

## Building this documentation site

The docs are written in Markdown and can be read directly in the repository. To build a searchable HTML site with [MkDocs Material](https://squidfunk.github.io/mkdocs-material/):

```bash
pip install -r docs/requirements.txt
mkdocs serve    # http://127.0.0.1:8000
mkdocs build    # site/ output
```

See `mkdocs.yml` at the repository root for site configuration.

