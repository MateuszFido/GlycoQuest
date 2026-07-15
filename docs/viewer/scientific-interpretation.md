# Scientific interpretation in the viewer

What the GlycoQuest viewer is designed to show — and what you should not over-interpret.

## Biological model

Each hit represents a **crosslinked pair of peptides** from your FASTA database:

- One peptide is a **glycopeptide** carrying an intact glycan on Asn, Ser, or Thr
- The other is typically an unmodified (or oxidized) peptide
- The **crosslinker** (e.g. DSS) connects specific side chains (e.g. Lys–Lys) — this is **not** the glycan–peptide bond

```
     Glycan (variable mod on N/S/T)
            │
    ────Peptide A────┐
                     │ DSS crosslink
    ────Peptide B────┘
```

The viewer's **crosslink edge** shows the crosslinker connection. The **glycan highlight** on the sequence map shows where the carbohydrate was attached in the searched model.

## What each panel proves

| Panel | Strong evidence for | Does not prove |
|-------|---------------------|----------------|
| QC funnel | Sample had glycan MS/MS signatures; prefilter worked | Completeness of glycoproteome |
| Network / table | Database peptides explain the crosslink at stated sites | False discovery rate |
| Sequence map | FASTA-localized sites when `mapped=true` | Localization confidence below residue level |
| Spectrum | Consistency of some b/y ions with sequences | Full de novo validation; glycan fragment assignment |

## Mapped vs unmapped hits

`mapped: true` when both peptides locate uniquely in `database.fasta` with consistent protein IDs.

| State | Network | Sequence map | Scientific use |
|-------|---------|--------------|----------------|
| Mapped | Edge drawn at AbsPos | Sites highlighted | Structural mapping on known proteins |
| Unmapped | Hidden from network | Limited | Peptide-level inspection only; fix FASTA headers |

Common mapping failures: UniProt-style headers in FASTA vs shortened names in xQuest output.

## Glycan fields in the viewer

| Field | Interpretation |
|-------|----------------|
| `glycan_composition` | Composition **searched** in the xQuest job that produced this hit |
| `glyco_residue` | Attachment type (N/S/T) after pseudo-residue decoding |
| `loss_label` | Water-loss search variant if applicable |
| `matched_families` (CSV) | Oxonium-ion families in the **experimental** spectrum (prefilter) |

A passing hit means: diagnostic evidence existed, xQuest found a crosslink consistent with that glycan mass, and hard post-filters passed. It is still a **database search result**, not independent glycan structure confirmation (no MS/MS of oxonium ions re-scored in xQuest).

## Spectrum plot and Filtering

- Peaks from **reduced** `spectra/` mzXML (post-prefilter subset)
- Diagnostic markers come from the completed-run Filtering record
- xQuest markers come from exact matched-ion rows emitted by xQuest V2.1.7
- Missing exact rows are shown as unavailable rather than inferred

## Topology

`topology` values (e.g. `intraprotein`, `interprotein`) classify whether both peptides come from the same FASTA entry or different entries — relevant for protein–protein interaction vs intra-molecular crosslinks.

## Comparing outputs

| Question | Best source |
|----------|-------------|
| Full column set + failed hits | `glycoquest_xquest.csv` |
| Integrated glycan QC + spectra | Viewer |
| CLMS-CSV network export | `xiview.csv` |
| Static QC slide | `report.html` |

## Confidence hierarchy

1. **Experimental glycan presence** — diagnostic ions in prefilter (`matched_families`, `matched_ion_count`)
2. **Crosslink identification** — xQuest `score`, exact xQuest rows when emitted
3. **Glycan identity** — composition consistent with job + mass error
4. **Site assignment** — sequon (N-glycans), database peptide context, not probabilistic localization

## Related

- [Glycopeptide crosslinking](../science/glycopeptide-crosslinking.md)
- [Post-filter](../workflow/postfilter.md)
