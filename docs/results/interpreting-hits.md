# Interpreting hits

`results/glycoquest_xquest.csv` is tab-separated. Each row is one **deduplicated** crosslink hit after GlycoQuest annotation and post-filtering.

Filter to `postfilter_status=pass` for high-confidence candidates unless you are auditing failures.

## Example row (illustrative)

| Column | Example value | Meaning |
|--------|---------------|---------|
| `source_file` | `sample.mzXML` | Originating mzXML |
| `scan` | `1896` | MS2 scan number |
| `glycan_name` | `HexNAc(2)Hex(5)` | Searched glycan |
| `glycan_composition` | `HexNAc(2)Hex(5)` | Composition string |
| `glycan_mass` | `910.328` | Variable-mod mass (Da) |
| `loss_label` | *(empty)* | Neutral-loss variant if not parent mass |
| `glyco_residue` | `N` | Attachment residue (Asn/Ser/Thr) |
| `glyco_peptide` | `1` | Which peptide carries glycan (1 or 2) |
| `n_glycan_pseudo` | `1` | Glycan pseudo count (must be 1 for pass) |
| `sequon_present` | `true` | N-X-S/T sequon at glycosite |
| `charge` | `5` | Precursor charge |
| `charge_plausible` | `true` | Charge in 2–7 |
| `matched_families` | `HexNAc` | Diagnostic families in spectrum |
| `matched_ion_count` | `3` | Diagnostic peaks matched |
| `seq1` | `MKXLCGR` | Peptide 1 (X = glycosylated Asn) |
| `seq2` | `KVEGLR` | Peptide 2 |
| `prot1` | `sp\|P00433\|PERX_HRP` | Protein accession 1 |
| `prot2` | `sp\|P00433\|PERX_HRP` | Protein accession 2 |
| `topology` | `intraprotein` | Intra- or inter-protein |
| `precursor_mz` | `823.45` | Observed precursor m/z |
| `mr` | `4112.2` | Neutral mass |
| `precursor_error_ppm` | `4.2` | Mass error |
| `xlink_position` | `2-4` | Crosslink sites in peptides |
| `score` | `12.5` | xQuest score |
| `hard_status` | `pass` | Hard filter outcome |
| `soft_score` | `13.600` | Ranking score |
| `postfilter_status` | `pass` | Final pass/fail |

## Column groups

### Provenance

| Column | Scientific meaning |
|--------|-------------------|
| `source_file` | Which LC–MS run |
| `scan` | Which MS/MS event — link to `spectra/` and viewer spectrum panel |

### Glycan identity

| Column | Scientific meaning |
|--------|-------------------|
| `glycan_name` | Library entry searched in the xQuest job |
| `glycan_composition` | Structural composition (HexNAc, Hex, Fuc, NeuAc counts) |
| `glycan_mass` | Monoisotopic mass used as variable modification |
| `loss_label` | If set, a water-loss (or other) search variant matched |

The glycan was **hypothesis-tested** per job: xQuest found a peptide pair consistent with that mass on the attachment residue.

### Glycosylation site

| Column | Scientific meaning |
|--------|-------------------|
| `glyco_residue` | `N`, `S`, or `T` — which residue type bears the glycan |
| `glyco_peptide` | `1` or `2` — which of the two crosslinked peptides is glycosylated |
| `n_glycan_pseudo` | Count of glycan pseudo-residues in `seq1`+`seq2` (V1 requires 1) |
| `sequon_present` | For Asn: classic N-glycan sequon N-X-S/T (X ≠ P) supports N-linked assignment |

Pseudo-residues in `seq1`/`seq2` (`X`, `U`, …) decode to the underlying residue + modification via the job manifest.

### Peptide identification

| Column | Scientific meaning |
|--------|-------------------|
| `seq1`, `seq2` | xQuest peptide sequences (with pseudo-residues) |
| `prot1`, `prot2` | FASTA protein identifiers |
| `topology` | `intraprotein` vs `interprotein` crosslink class |
| `xlink_position` | Residue indices of the crosslink within each peptide |

The **crosslink** connects the configured chemistry sites (e.g. Lys–Lys for DSS). The **glycan** is a separate modification on one peptide.

### Mass spectrometry QC

| Column | Scientific meaning |
|--------|-------------------|
| `precursor_mz` | Experimental precursor |
| `mr` | Calculated neutral mass |
| `precursor_error_ppm` | Parts-per-million error — lower is better |
| `charge` | Precursor charge state |
| `charge_plausible` | Whether charge is in the expected 2–7 range |

### Evidence and confidence

| Column | Scientific meaning |
|--------|-------------------|
| `matched_families` | Oxonium-ion families observed in prefilter for this scan |
| `matched_ion_count` | Number of diagnostic peaks matched — more supports glycan presence |
| `score` | xQuest crosslink score |
| `hard_status` | Which hard rule failed, if any |
| `soft_score` | Composite rank including sequon, charge, diagnostics, mass penalty |
| `postfilter_status` | `pass` or `fail` — use for downstream tables |

## Hard status codes

| Value | Interpretation |
|-------|----------------|
| `pass` | Eligible passing hit |
| `fail_no_xlink` | Not a confident crosslink in xQuest output |
| `fail_no_glycan` | xQuest assignment did not use the glycan modification |
| `fail_multiple_glycans` | Assignment contains more than one glycan modification |
| `fail_no_diagnostic` | Spectrum lacks glycan diagnostic support |
| `fail_precursor_error` | Precursor mass inconsistent |
| `fail_score` | Below `min_score` threshold |

## Related

- [Post-filter](../workflow/postfilter.md)
- [Scientific interpretation in the viewer](../viewer/scientific-interpretation.md)
- [Network CSV export](network-csv-export.md)
