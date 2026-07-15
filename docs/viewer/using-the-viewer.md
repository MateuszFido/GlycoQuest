# Using the viewer

The GlycoQuest viewer is an offline interactive bundle in `results/viewer/`. It coordinates QC charts, a crosslink network, a hits table, sequence maps, and MS/MS mirror plots.

![Viewer layout](../assets/viewer-layout.svg)

## Launch

**Serve the folder over HTTP:**

```bash
cd out/<project>/results/viewer
python3 -m http.server 8080
```

Open [http://localhost:8080](http://localhost:8080).

Or use the bundled script:

```bash
./serve-viewer.sh 8080
```

## Files loaded

| File | Role |
|------|------|
| `index.html` | Entry page |
| `viewer.json` | All hit, protein, QC, and spectrum data |
| `database.fasta` | Protein sequences for sequence map |
| `viewer.js`, `viewer.css` | Bundled UI assets |

No internet connection required after the page loads.

## Header

Shows project name, input file, crosslinker chemistry, glycan library, xQuest version when available, and passing/total hit counts.

## Toolbar

| Control | Effect |
|---------|--------|
| **Protein** dropdown | Filter crosslinks to those involving the selected protein |
| **Show failed hits** | Include hits with `postfilter_status=fail` |
| **Min score** | Hide hits below xQuest score threshold |

## Panels

### QC panel

Mirrors `report.html` content interactively:

- Prefilter funnel (scans → diagnostic → pairs → filtered)
- Post-filter outcome breakdown
- Top glycan compositions and glycosylation sites
- Score and precursor-error histograms

Use for run-quality assessment before drilling into individual hits.

### Sequence Pair Map

The pair map shows the selected protein pair with compact protein labels, a full-protein coordinate overview, and a focused amino-acid letter track around the selected crosslink context. Crosslink stacks at the same endpoint open a chooser so the selected hit matches the user's click.

Glycan chips and glycan-site markers use SNFG-style inline symbols derived from the glycan composition. Click a glycan chip or marker to expand the representation.

If empty: *"No mapped sequence coordinates"* — FASTA IDs may not match xQuest protein names. Hits can still appear in the table.

### Hits table

List of visible crosslinks. Click a row to select it and update other panels. Glycan cells show SNFG chips when a composition is available.

Columns include scan, retention time (minutes), proteins, glycan, score, and post-filter status.

### Spectrum panel

Experimental MS/MS plot (x = **m/z**, y = relative intensity) from reduced `spectra/` mzXML. Diagnostic markers and xQuest markers are drawn only when the completed-run Filtering record contains the exact peak evidence. Scan metadata above the plot is `Scan | precursor | charge | scan_time` — retention/scan time is never plotted on an axis.

!!! note
    If xQuest did not emit exact matched-ion rows, the viewer reports that absence in Filtering instead of inventing annotations.

### Filtering panel

Step-by-step record for the selected crosslink:

- Input scan
- Diagnostic prefilter
- Isotope evidence
- Glycan pruning
- xQuest search
- Postfilter

Each step shows source artifact names, rows when available, scan ids, peak indices, counts, thresholds, and pass/fail states.

## Selection workflow

1. Open viewer → first passing hit auto-selected (if any)
2. Scan QC panel for funnel health
3. Browse network or table
4. Select a hit → inspect sequence map, spectrum, and Filtering together
5. Toggle **Show failed hits** to audit near-misses

## Rebuilding viewer assets

If you modify `viewer/` TypeScript sources:

```bash
cd viewer && npm install && npm run build
```

Then re-run GlycoQuest to copy assets into `results/viewer/`.

## Related

- [Scientific interpretation](scientific-interpretation.md)
- [HTML report](../results/report.md)
- [Interpreting hits](../results/interpreting-hits.md)
