# Network CSV export

`results/xiview.csv` contains **passing** crosslinks in the standard **CLMS-CSV** column layout used by crosslink network visualization tools.

## Standard columns

| Column | Description |
|--------|-------------|
| `Protein1` | Protein identifier |
| `PepPos1` | 1-based start position of peptide 1 in protein |
| `PepSeq1` | Peptide 1 sequence (pseudo-residues resolved) |
| `LinkPos1` | Crosslink position within peptide 1 |
| `AbsPos1` | Absolute residue number in protein: `PepPos + LinkPos - 1` |
| `Protein2` … `AbsPos2` | Same for peptide 2 |
| `Score` | xQuest score |

GlycoQuest locates each peptide in the FASTA to compute `PepPos*` and `AbsPos*`.

## GlycoQuest extensions

Additional columns are appended at the end. Conforming CLMS-CSV importers ignore unrecognized columns:

- Glycan composition
- Glycosylation residue
- Loss label

## Mapped vs unmapped

Terminal output reports:

```text
wrote results/xiview.csv (N passing crosslink(s), M with resolved absolute positions)
```

- **N** — rows written (passing post-filter)
- **M** — rows where both peptides mapped to FASTA entries

Unmapped rows occur when protein IDs in xQuest output do not match FASTA headers. Fix FASTA header formatting or use the built-in [viewer](../viewer/using-the-viewer.md) for unmapped hits.

## Related

- [Interpreting hits](interpreting-hits.md)
- [Using the viewer](../viewer/using-the-viewer.md)
