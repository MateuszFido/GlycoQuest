# Diagnostic ions

GlycoQuest's first prefilter gate detects **glycan fragment ions** in MS/MS spectra — typically oxonium ions and their neutral-loss variants.

## Oxonium ions

When glycopeptides fragment in collision-induced dissociation, carbohydrate portions can cleave as charged monosaccharide or disaccharide ions. These appear at characteristic m/z values and define **diagnostic families**.

| Family | Example ion | Approx. m/z | Indicates |
|--------|-------------|-------------|-----------|
| HexNAc | N-acetylhexosamine | ~204.09 | HexNAc-containing structures |
| Hex | Hexose | ~163.06 | Hexose-containing structures |
| NeuAc | N-acetylneuraminic acid | ~292.10 | Sialylated structures |
| Fuc | Fucose | ~147.07 | Fucosylated structures |

Exact masses depend on the catalog and neutral-loss expansion in bundled data (`diagnostic_ion_catalog.txt`).

## Neutral losses

The same family can appear at lower m/z after neutral losses (e.g. −H₂O, −2H₂O). GlycoQuest expands search targets from the catalog:

```text
HexNAc@204.0867
HexNAc@186.0760[-H2O]
HexNAc@168.0654[-2H2O]
```

Matching any expanded target counts toward the **family** for pruning.

## Matching algorithm

For each MS2 peak:

1. Compare observed m/z to each library target
2. Accept if within `diagnostic_tolerance_ppm` (default 10 ppm)
3. Record matched family and ion details
4. Spectrum **passes** if ≥1 family matched

**Intensity:** any peak above baseline counts — no minimum relative abundance.

## Glycan library linkage

Each glycan entry lists which diagnostic ions it **can** produce. The prefilter uses the **union** across the library for matching, then **prunes** to glycans whose **composition** contains every family seen in that spectrum.

Example:

- Spectrum shows HexNAc + NeuAc ions
- Pruning keeps only glycans containing both HexNAc and NeuAc in their composition

See [Prefilter — glycan pruning](../workflow/prefilter.md#gate-3-glycan-pruning).

## Custom libraries

When authoring `diagnostic_ions` in a custom CSV:

```text
NeuAc@292.1027;NeuAc@274.0921[-H2O]
```

Use `family@mz` syntax; optional `[-loss_label]` for catalog-backed losses.

## Scientific limits

Diagnostic ions demonstrate **class-level** glycan presence (e.g. sialylation) — not full structural assignment. GlycoQuest uses them to:

- Gate spectra before database search
- Support hard post-filter (`fail_no_diagnostic`)
- Boost soft score via `matched_ion_count`

They do not replace dedicated glycomics fragmentation analysis for novel structures.

## Related

- [Glycan libraries](../configuration/glycan-libraries.md)
- [Prefilter](../workflow/prefilter.md)
- [Interpreting hits](../results/interpreting-hits.md)
