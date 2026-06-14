//! FASTA database validation for xQuest compatibility.

use std::path::{Path, PathBuf};

const RESERVED_PSEUDO_RESIDUES: &[char] = &['X', 'U', 'B', 'J'];

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaEntry {
    pub header: String,
    pub sequence: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaDatabase {
    pub path: PathBuf,
    pub entries: Vec<FastaEntry>,
}

/// Parse and validate a protein FASTA file.
pub fn validate_fasta(path: &Path) -> Result<FastaDatabase, String> {
    if !path.is_file() {
        return Err(format!("FASTA database not accessible: {}", path.display()));
    }

    let content = std::fs::read_to_string(path).map_err(|err| {
        format!("cannot read FASTA database {}: {err}", path.display())
    })?;

    let mut entries = Vec::new();
    let mut header: Option<String> = None;
    let mut sequence = String::new();

    for (line_no, raw_line) in content.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }

        if line.starts_with('>') {
            if let Some(prev_header) = header.take() {
                let seq = std::mem::take(&mut sequence);
                validate_entry(&prev_header, &seq, line_no)?;
                entries.push(FastaEntry {
                    header: prev_header,
                    sequence: seq,
                });
            }
            header = Some(line[1..].trim().to_string());
            if header.as_ref().is_some_and(String::is_empty) {
                return Err(format!(
                    "empty FASTA header on line {} in {}",
                    line_no + 1,
                    path.display()
                ));
            }
            continue;
        }

        if header.is_none() {
            return Err(format!(
                "FASTA sequence before header on line {} in {}",
                line_no + 1,
                path.display()
            ));
        }

        if !line.chars().all(|c| c.is_ascii_alphabetic()) {
            return Err(format!(
                "invalid FASTA residue on line {} in {}: {line}",
                line_no + 1,
                path.display()
            ));
        }
        sequence.push_str(line);
    }

    if let Some(prev_header) = header {
        validate_entry(&prev_header, &sequence, content.lines().count())?;
        entries.push(FastaEntry {
            header: prev_header,
            sequence,
        });
    }

    if entries.is_empty() {
        return Err(format!("FASTA database contains no entries: {}", path.display()));
    }

    let path = path
        .canonicalize()
        .unwrap_or_else(|_| path.to_path_buf());

    Ok(FastaDatabase {
        path: path.to_path_buf(),
        entries,
    })
}

fn validate_entry(header: &str, sequence: &str, line_no: usize) -> Result<(), String> {
    if sequence.is_empty() {
        return Err(format!(
            "empty FASTA sequence for entry {header} near line {line_no}"
        ));
    }

    for ch in sequence.chars() {
        if RESERVED_PSEUDO_RESIDUES.contains(&ch) {
            return Err(format!(
                "FASTA entry {header} contains reserved pseudo-residue '{ch}' \
                 on line {line_no} \
                 (xQuest reserves X, U, B, J for variable modifications)"
            ));
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::sync::atomic::{AtomicU64, Ordering};

    static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

    fn temp_fasta(name: &str, content: &str) -> PathBuf {
        let id = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "glycoquest_fasta_test_{}_{}_{}.fasta",
            std::process::id(),
            name,
            id
        ));
        fs::write(&path, content).unwrap();
        path
    }

    #[test]
    fn accepts_valid_fasta() {
        let path = temp_fasta(
            "valid",
            ">protein1\nACDEFG\n>protein2\nHIKLMN\n",
        );
        let db = validate_fasta(&path).unwrap();
        assert_eq!(db.entries.len(), 2);
        assert_eq!(db.entries[0].sequence, "ACDEFG");
    }

    #[test]
    fn accepts_wrapped_sequences() {
        let path = temp_fasta(
            "wrapped",
            ">protein\nACDEF\nGHIKL\n",
        );
        let db = validate_fasta(&path).unwrap();
        assert_eq!(db.entries[0].sequence, "ACDEFGHIKL");
    }

    #[test]
    fn rejects_empty_sequence() {
        let path = temp_fasta("empty", ">protein\n");
        let err = validate_fasta(&path).unwrap_err();
        assert!(err.contains("empty FASTA sequence"));
    }

    #[test]
    fn rejects_reserved_pseudo_residue() {
        let path = temp_fasta("reserved", ">protein\nACDXFG\n");
        let err = validate_fasta(&path).unwrap_err();
        assert!(err.contains("reserved pseudo-residue"));
        assert!(err.contains("'X'"));
    }

    #[test]
    fn rejects_missing_file() {
        let path = std::env::temp_dir().join(format!(
            "glycoquest_fasta_missing_{}_{}.fasta",
            std::process::id(),
            TEMP_COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        let err = validate_fasta(&path).unwrap_err();
        assert!(err.contains("not accessible") || err.contains("not found"));
    }

    #[test]
    fn bundled_hrp_fasta_rejects_contaminant_with_reserved_x() {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/rcsb_pdb_1HRP.fasta");
        if path.is_file() {
            let err = validate_fasta(&path).unwrap_err();
            assert!(err.contains("reserved pseudo-residue"));
        }
    }
}
