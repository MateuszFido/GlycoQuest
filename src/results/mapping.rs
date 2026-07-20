// Copyright (c) ETH Zurich, Mateusz Fido

//! Shared FASTA / peptide position mapping for network CSV export and the interactive viewer.

use std::collections::HashMap;

use crate::fasta::FastaDatabase;
use crate::jobs::VarModPlan;

/// Map protein identifiers to sequences, indexed by every plausible key xQuest
/// might report (full header, first whitespace token, and `|`-delimited parts).
pub fn protein_lookup(fasta: &FastaDatabase) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for entry in &fasta.entries {
        let seq = entry.sequence.clone();
        map.entry(entry.header.clone())
            .or_insert_with(|| seq.clone());
        if let Some(token) = entry.header.split_whitespace().next() {
            map.entry(token.to_string()).or_insert_with(|| seq.clone());
            for part in token.split('|').filter(|p| !p.is_empty()) {
                map.entry(part.to_string()).or_insert_with(|| seq.clone());
            }
        }
    }
    map
}

/// Resolve xQuest pseudo-residues back to source residues for FASTA lookup.
pub fn resolve_peptide(seq: &str, plan: Option<&VarModPlan>) -> String {
    seq.chars()
        .filter_map(|ch| {
            if !ch.is_ascii_alphabetic() {
                return None;
            }
            let resolved = plan
                .and_then(|plan| plan.entry_for_pseudo(ch))
                .map(|entry| entry.source_residue)
                .unwrap_or(ch);
            Some(resolved.to_ascii_uppercase())
        })
        .collect()
}

/// First protein identifier from a possibly comma/semicolon-delimited list.
pub fn first_protein(prot: &str) -> &str {
    prot.split([',', ';']).next().unwrap_or(prot).trim()
}

/// 1-based start position of `peptide` within its protein, if both are known.
pub fn locate_peptide(
    proteins: &HashMap<String, String>,
    protein: &str,
    peptide: &str,
) -> Option<usize> {
    locate_peptide_resolved(proteins, protein, peptide).map(|located| located.start)
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct LocatedPeptide {
    start: usize,
    sequence: String,
}

fn locate_peptide_resolved(
    proteins: &HashMap<String, String>,
    protein: &str,
    peptide: &str,
) -> Option<LocatedPeptide> {
    if peptide.is_empty() {
        return None;
    }
    let seq = proteins.get(protein)?;
    if let Some(idx) = seq.find(peptide) {
        return Some(LocatedPeptide {
            start: idx + 1,
            sequence: peptide.to_string(),
        });
    }

    locate_with_pseudo_wildcards(seq, peptide)
}

fn locate_with_pseudo_wildcards(protein_seq: &str, peptide: &str) -> Option<LocatedPeptide> {
    let peptide_chars: Vec<char> = peptide.chars().collect();
    if peptide_chars.iter().all(|ch| is_canonical_residue(*ch)) {
        return None;
    }

    let protein_chars: Vec<char> = protein_seq.chars().collect();
    if peptide_chars.len() > protein_chars.len() {
        return None;
    }

    for start in 0..=protein_chars.len() - peptide_chars.len() {
        let candidate = &protein_chars[start..start + peptide_chars.len()];
        let matches = peptide_chars
            .iter()
            .zip(candidate.iter())
            .all(|(pep, prot)| !is_canonical_residue(*pep) || pep.eq_ignore_ascii_case(prot));
        if matches {
            return Some(LocatedPeptide {
                start: start + 1,
                sequence: candidate.iter().collect(),
            });
        }
    }

    None
}

fn is_canonical_residue(ch: char) -> bool {
    matches!(
        ch.to_ascii_uppercase(),
        'A' | 'C'
            | 'D'
            | 'E'
            | 'F'
            | 'G'
            | 'H'
            | 'I'
            | 'K'
            | 'L'
            | 'M'
            | 'N'
            | 'P'
            | 'Q'
            | 'R'
            | 'S'
            | 'T'
            | 'V'
            | 'W'
            | 'Y'
    )
}

/// Parse `"a-b"` crosslink positions within each peptide (1-based).
pub fn parse_link_positions(xlink_position: &str) -> (Option<usize>, Option<usize>) {
    let mut parts = xlink_position.split(['-', ',']).map(str::trim);
    let first = parts.next().and_then(|s| s.parse::<usize>().ok());
    let second = parts.next().and_then(|s| s.parse::<usize>().ok());
    (first, second)
}

/// Absolute 1-based residue position in the protein sequence.
pub fn abs_position(pep_pos: Option<usize>, link_pos: Option<usize>) -> Option<usize> {
    match (pep_pos, link_pos) {
        (Some(p), Some(l)) if l >= 1 => Some(p + l - 1),
        _ => None,
    }
}

/// Resolved crosslink residue mapping for one annotated hit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CrosslinkMapping {
    pub prot1: String,
    pub prot2: String,
    pub pep1: String,
    pub pep2: String,
    pub pep_pos1: Option<usize>,
    pub pep_pos2: Option<usize>,
    pub link1: Option<usize>,
    pub link2: Option<usize>,
    pub abs1: Option<usize>,
    pub abs2: Option<usize>,
}

pub fn map_crosslink(
    seq1: &str,
    seq2: &str,
    prot1: &str,
    prot2: &str,
    xlink_position: &str,
    proteins: &HashMap<String, String>,
    plan: Option<&VarModPlan>,
) -> CrosslinkMapping {
    let pep1 = resolve_peptide(seq1, plan);
    let pep2 = resolve_peptide(seq2, plan);
    let (link1, link2) = parse_link_positions(xlink_position);
    let p1 = first_protein(prot1);
    let p2 = first_protein(prot2);
    let located1 = locate_peptide_resolved(proteins, p1, &pep1);
    let located2 = locate_peptide_resolved(proteins, p2, &pep2);
    let pep_pos1 = located1.as_ref().map(|located| located.start);
    let pep_pos2 = located2.as_ref().map(|located| located.start);
    let pep1 = located1.map(|located| located.sequence).unwrap_or(pep1);
    let pep2 = located2.map(|located| located.sequence).unwrap_or(pep2);
    CrosslinkMapping {
        prot1: p1.to_string(),
        prot2: p2.to_string(),
        pep1,
        pep2,
        pep_pos1,
        pep_pos2,
        link1,
        link2,
        abs1: abs_position(pep_pos1, link1),
        abs2: abs_position(pep_pos2, link2),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::{FastaDatabase, FastaEntry};
    use std::path::PathBuf;

    fn fasta() -> FastaDatabase {
        FastaDatabase {
            path: PathBuf::from("test.fasta"),
            entries: vec![
                FastaEntry {
                    header: "sp|P00761|TRYP_PIG trypsin".into(),
                    sequence: "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHR".into(),
                },
                FastaEntry {
                    header: "HRP".into(),
                    sequence: "QLTPTFYDNSCPNVSNIVRDTIVNELR".into(),
                },
            ],
        }
    }

    #[test]
    fn locate_peptide_finds_1_based_start() {
        let proteins = protein_lookup(&fasta());
        assert_eq!(locate_peptide(&proteins, "P00761", "DTHK"), Some(25));
    }

    #[test]
    fn map_crosslink_resolves_absolute_positions() {
        let proteins = protein_lookup(&fasta());
        let m = map_crosslink("DTHK", "DTIVNELR", "P00761", "HRP", "2-3", &proteins, None);
        assert_eq!(m.abs1, Some(26));
        assert_eq!(m.abs2, Some(22));
    }

    #[test]
    fn map_crosslink_recovers_pseudo_residue_peptides_without_varmod_plan() {
        let proteins = HashMap::from([("P1".to_string(), "MMAKVFKDVFLEMNIPYSVVRQQ".to_string())]);

        let m = map_crosslink(
            "AKVFKDVFLEUXIPYSVVR",
            "QQ",
            "P1",
            "P1",
            "5-1",
            &proteins,
            None,
        );

        assert_eq!(m.pep1, "AKVFKDVFLEMNIPYSVVR");
        assert_eq!(m.pep_pos1, Some(3));
        assert_eq!(m.abs1, Some(7));
    }
}
