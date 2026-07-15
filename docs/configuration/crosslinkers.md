# Crosslinkers

GlycoQuest searches **crosslinked peptide pairs** where one peptide carries an intact **glycan**. The crosslinker chemistry determines MS1 isotope behavior, xQuest mass parameters, and allowed amino-acid site pairs.

## Bundled presets

Select with `--crosslinker <name>` or set `[crosslinker] name` in `settings.ini`.

### DSS (default)

Disuccinimidyl suberate — homobifunctional linker targeting primary amines (typically Lys).

| Parameter | Value |
|-----------|-------|
| Preset id | `dss` |
| Label mode | `light-heavy` |
| Mass shift | 12.075321 Da |
| xlinkermw | 138.0680796 Da |
| Sites | `K:K` |
| Isotope prefilter | **Enabled** |

DSS duplex experiments produce light- and heavy-labeled precursor pairs separated by the shift mass. GlycoQuest requires a matching isotope partner before sending a spectrum to xQuest.

```bash
glycoquest run.mzXML --database proteins.fasta --crosslinker dss
```

### DMTMM

4-(4,6-dimethoxy-1,3,5-triazin-2-yl)-4-methylmorpholinium — a coupling reagent that produces a zero-length amide between a Lys amine and an Asp/Glu carboxyl group.

| Parameter | Value |
|-----------|-------|
| Preset id | `dmtmm` |
| Label mode | `none` |
| Mass shift | 0 Da |
| xlinkermw | −18.0109 Da (net water loss) |
| Sites | `K:E,K:D` |
| Isotope prefilter | **Disabled** |

```bash
glycoquest run.mzXML --database proteins.fasta --crosslinker dmtmm
```

### Published glycan-to-peptide linkers

These presets encode a physical SiaNAz–linker–Lys bridge as xQuest `X:K`.
GlycoQuest puts the glycan in the first variable-modification slot, which xQuest
rewrites to pseudo-residue `X`; this forces the glycan and crosslink onto the
same Asn. They are unlabeled, so isotope pairing is disabled.

| Preset | Reference/sample state | `xlinkermw` | Formula relative to NeuAc glycan | Sites |
|--------|------------------------|------------:|-----------------------------------|-------|
| `nhs-cyclooctyne` | Xie et al. 2021 intact NHS–cyclooctyne product | 205.085126607 | C10H11N3O2 | `X:K` |
| `ssbxl` | Chen et al. 2025 after TCEP and IAA | 573.179438173 | C28H27N7O5S | `X:K` |
| `pcbxl` | Chen et al. 2025 after photocleavage | 456.190988659 | C25H24N6O3 | `X:K` |

The formulas and masses assume that the glycan library contains ordinary NeuAc;
the preset mass includes the SiaNAz-for-NeuAc delta. Do not reuse these values
unchanged with a library whose glycan masses already contain SiaNAz. For SSBXL,
573.179438173 + 291.095416527 (NeuAc residue) + 1.007276467 (H+) gives m/z
865.2821, matching the linker's required signature ion in Chen et al.

```bash
glycoquest run.mzXML --database proteins.fasta \
  --crosslinker nhs-cyclooctyne \
  --glycans examples/MSV000087442/glycans.csv
```

See the executed repository example at `examples/MSV000087442/README.md` and
[the chemistry derivation](../science/glycopeptide-crosslinking.md#paper-derived-search-masses).

## Labeling modes (`[crosslinker] label`)

| Mode | Isotope prefilter | xQuest isotope settings | Typical use |
|------|-------------------|-------------------------|-------------|
| `light-heavy` | On | `isotopeshift 12.075321`, isotopic scan pairs | DSS duplex (default) |
| `light-only` | Off | Light-only pair printing | Single-channel DSS |
| `none` | Off | `isotopeshift 0` | DMTMM, GPx linkers, EDC, other unlabeled chemistries |

When `label=none`, `shift_da` is ignored for the isotope prefilter (a warning is printed if shift ≠ 0).

## Custom crosslinker chemistry

For chemistries without a preset, set fields manually in `settings.ini`:

```ini
[crosslinker]
name = mylinker
label = none
shift_da = 0
xlinkermw = 123.4567
xlink_sites = K:K
nterm_xlinkable = false
```

Validate generated `tmp/jobs/*/xquest.def` in a dry-run before long runs.

## N-terminal crosslinking

Set `nterm_xlinkable = true` to enable xQuest N-terminus support (`ntermxlinkable 1`) and `Z` pseudo-sites in `xlink_sites` (e.g. `K:K,K:Z,Z:Z`). Integration testing is recommended before relying on this in production.

## Scientific distinction: crosslink vs glycan

For DSS and DMTMM, the **crosslinker** connects two peptide chains while the
**glycan** is a separate modification on Asn/Ser/Thr. For the GPx presets, the
physical crosslink runs from SiaNAz within that glycan to a peptide Lys; xQuest's
`X:K` site pair is the computational encoding of that relationship.

```
Peptide A ——[DSS]—— Peptide B
              |
         (glycan on Asn
          of Peptide A)
```

See [Glycopeptide crosslinking](../science/glycopeptide-crosslinking.md).

## Related

- [Prefilter — isotope pairs](../workflow/prefilter.md#gate-2-dss-isotope-pairing)
- [xQuest jobs — crosslinker defs](../workflow/xquest-jobs.md)
- [Settings reference](settings-reference.md)
