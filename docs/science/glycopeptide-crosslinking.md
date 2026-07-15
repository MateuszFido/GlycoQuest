# Glycopeptide crosslinking


## The biological question

Glycoproteins carry oligosaccharides on specific amino acid side chains. In **crosslinking mass spectrometry (XL-MS)**, a bifunctional reagent connects two peptides that are close in space. **GlycoQuest** targets the case where:

- One peptide in the pair carries an **intact glycan**
- The other peptide is typically unglycosylated
- A **chemical crosslinker** (e.g., DSS) connects the two peptide backbones

This reveals **spatial proximity** between a glycosylated region and another peptide — useful for epitope mapping, glycoprotein architecture, and interaction studies.

For isotopically labelled crosslinkers, heavy/light isotope prefiltering applies to the **crosslinker**, not the glycan. 

## N-linked vs O-linked glycans

| Type | Attachment | Sequon | Default library |
|------|------------|--------|-----------------|
| N-linked | Asparagine (Asn) | N-X-S/T (X ≠ Pro) | `nglyc309` |
| O-linked | Serine / Threonine | No strict sequon | `oglyc78` |

GlycoQuest's score for N-glycan matches is influenced by sequon presence; O-glycans rely on diagnostic ions and database context.

## Why prefilter before xQuest?

Full crosslink searches are computationally expensive. Glycopeptide spectra often show characteristic **oxonium ions** (e.g. HexNAc at ~204 m/z). Requiring diagnostic evidence:

- Reduces false searches on non-glycan spectra
- Prunes the glycan library per spectrum
- Focuses xQuest on plausible glycopeptide–peptide candidates

## Typical experimental setup (DSS)

1. Proteins crosslinked with DSS (light/heavy mixture)
2. Enzymatic digest (trypsin)
3. LC–MS/MS
4. mzXML → GlycoQuest → xQuest

See [Crosslinkers](../configuration/crosslinkers.md) for other options.

## Related

- [Diagnostic ions](diagnostic-ions.md)
- [Workflow overview](../workflow/overview.md)
- [Interpreting hits](../results/interpreting-hits.md)
