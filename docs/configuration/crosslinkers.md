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

Dimethyl 3,3′-dithiobispropionimidate — zero-length crosslink between Lys and Asp/Glu side chains.

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

## Labeling modes (`[crosslinker] label`)

| Mode | Isotope prefilter | xQuest isotope settings | Typical use |
|------|-------------------|-------------------------|-------------|
| `light-heavy` | On | `isotopeshift 12.075321`, isotopic scan pairs | DSS duplex (default) |
| `light-only` | Off | Light-only pair printing | Single-channel DSS |
| `none` | Off | `isotopeshift 0` | DMTMM, EDC, other unlabeled chemistries |

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

The **crosslinker** connects two peptide chains (e.g. Lys–Lys for DSS). The **glycan** is a separate modification on Asn/Ser/Thr, modeled as an xQuest variable modification — it is not part of the crosslink bond and is not isotope-coded by DSS.

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
