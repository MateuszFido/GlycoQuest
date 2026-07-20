// Copyright (c) ETH Zurich, Mateusz Fido

//! Load an explicit user-provided glycan library file (CSV/TSV).
//!
//! Schema (DESIGN.md section 4):
//! `name,composition,monoisotopic_mass,diagnostic_ions,residue_targets`
//! where `diagnostic_ions` is a `;`-separated list of `family@mz` entries with an
//! optional `[-loss]` suffix, and `residue_targets` is a `;`-separated list of
//! residue letters.

use std::path::Path;

use super::composition::canonical_residue;
use super::{DiagnosticIon, GlycanEntry, GlycanLibrary};

const REQUIRED_COLUMNS: [&str; 5] = [
    "name",
    "composition",
    "monoisotopic_mass",
    "diagnostic_ions",
    "residue_targets",
];

/// Load a glycan library from an explicit CSV or TSV file.
pub fn load_glycan_library_file(path: &Path) -> Result<GlycanLibrary, String> {
    let content = std::fs::read_to_string(path)
        .map_err(|err| format!("cannot read glycan library {}: {err}", path.display()))?;

    let mut lines = content
        .lines()
        .map(str::trim_end)
        .filter(|line| !line.trim().is_empty() && !line.trim_start().starts_with('#'));

    let header = lines
        .next()
        .ok_or_else(|| format!("glycan library {} is empty", path.display()))?;

    let delimiter = sniff_delimiter(header).ok_or_else(|| {
        format!(
            "unsupported delimiter in glycan library {} (use CSV or TSV)",
            path.display()
        )
    })?;

    let columns = column_indices(header, delimiter, path)?;

    let mut entries = Vec::new();
    let mut seen = std::collections::HashSet::new();

    for (row_index, line) in lines.enumerate() {
        let fields: Vec<&str> = line.split(delimiter).map(str::trim).collect();
        let line_no = row_index + 2; // header is line 1 (after comment/blank filtering)

        let name = field(&fields, columns.name, "name", line_no, path)?.to_string();
        if !seen.insert(name.clone()) {
            return Err(format!(
                "duplicate glycan name '{name}' on line {line_no} in {}",
                path.display()
            ));
        }

        let composition =
            field(&fields, columns.composition, "composition", line_no, path)?.to_string();

        let mass_str = field(&fields, columns.mass, "monoisotopic_mass", line_no, path)?;
        let monoisotopic_mass: f64 = mass_str.parse().map_err(|_| {
            format!(
                "non-numeric monoisotopic_mass '{mass_str}' on line {line_no} in {}",
                path.display()
            )
        })?;
        if monoisotopic_mass <= 0.0 {
            return Err(format!(
                "non-positive monoisotopic_mass on line {line_no} in {}",
                path.display()
            ));
        }

        let diag_str = field(
            &fields,
            columns.diagnostic_ions,
            "diagnostic_ions",
            line_no,
            path,
        )?;
        let diagnostic_ions = parse_diagnostic_ions(diag_str, line_no, path)?;
        if diagnostic_ions.is_empty() {
            return Err(format!(
                "empty diagnostic_ions on line {line_no} in {} (V1 must not skip diagnostic filtering)",
                path.display()
            ));
        }

        let targets_str = field(
            &fields,
            columns.residue_targets,
            "residue_targets",
            line_no,
            path,
        )?;
        let residue_targets = parse_residue_targets(targets_str, line_no, path)?;

        entries.push(GlycanEntry {
            name,
            composition,
            monoisotopic_mass,
            diagnostic_ions,
            residue_targets,
        });
    }

    if entries.is_empty() {
        return Err(format!(
            "glycan library {} has no glycan rows",
            path.display()
        ));
    }

    let database_id = path
        .file_stem()
        .and_then(|stem| stem.to_str())
        .unwrap_or("custom")
        .to_string();

    Ok(GlycanLibrary {
        database_id,
        entries,
    })
}

struct ColumnIndices {
    name: usize,
    composition: usize,
    mass: usize,
    diagnostic_ions: usize,
    residue_targets: usize,
}

fn sniff_delimiter(header: &str) -> Option<char> {
    if header.contains('\t') {
        Some('\t')
    } else if header.contains(',') {
        Some(',')
    } else {
        None
    }
}

fn column_indices(header: &str, delimiter: char, path: &Path) -> Result<ColumnIndices, String> {
    let headers: Vec<String> = header
        .split(delimiter)
        .map(|h| h.trim().to_ascii_lowercase())
        .collect();

    let find = |name: &str| {
        headers
            .iter()
            .position(|h| h == name)
            .ok_or_else(|| format!("missing required column '{name}' in {}", path.display()))
    };

    for column in REQUIRED_COLUMNS {
        find(column)?;
    }

    Ok(ColumnIndices {
        name: find("name")?,
        composition: find("composition")?,
        mass: find("monoisotopic_mass")?,
        diagnostic_ions: find("diagnostic_ions")?,
        residue_targets: find("residue_targets")?,
    })
}

fn field<'a>(
    fields: &'a [&'a str],
    index: usize,
    name: &str,
    line_no: usize,
    path: &Path,
) -> Result<&'a str, String> {
    let value = fields.get(index).copied().unwrap_or("").trim();
    if value.is_empty() {
        return Err(format!(
            "missing '{name}' value on line {line_no} in {}",
            path.display()
        ));
    }
    Ok(value)
}

fn parse_diagnostic_ions(
    raw: &str,
    line_no: usize,
    path: &Path,
) -> Result<Vec<DiagnosticIon>, String> {
    let mut ions = Vec::new();
    for token in raw.split(';').map(str::trim).filter(|t| !t.is_empty()) {
        let (family_raw, rest) = token.split_once('@').ok_or_else(|| {
            format!(
                "invalid diagnostic ion '{token}' on line {line_no} in {} (expected family@mz)",
                path.display()
            )
        })?;

        let (mz_str, loss_label) = match rest.split_once('[') {
            Some((mz, label)) => {
                let label = label.trim_end_matches(']').trim();
                (mz.trim(), label.to_string())
            }
            None => (rest.trim(), String::new()),
        };

        let mz: f64 = mz_str.parse().map_err(|_| {
            format!(
                "non-numeric diagnostic ion m/z '{mz_str}' on line {line_no} in {}",
                path.display()
            )
        })?;
        if mz <= 0.0 {
            return Err(format!(
                "non-positive diagnostic ion m/z on line {line_no} in {}",
                path.display()
            ));
        }

        ions.push(DiagnosticIon {
            family: canonical_residue(family_raw.trim()).to_string(),
            mz,
            loss_label,
        });
    }
    Ok(ions)
}

fn parse_residue_targets(raw: &str, line_no: usize, path: &Path) -> Result<Vec<String>, String> {
    let targets: Vec<String> = raw
        .split(';')
        .map(str::trim)
        .filter(|t| !t.is_empty())
        .map(str::to_string)
        .collect();
    if targets.is_empty() {
        return Err(format!(
            "empty residue_targets on line {line_no} in {}",
            path.display()
        ));
    }
    Ok(targets)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn temp_file(name: &str, content: &str) -> std::path::PathBuf {
        let path =
            std::env::temp_dir().join(format!("glycoquest_lib_{}_{}", std::process::id(), name));
        let mut file = std::fs::File::create(&path).unwrap();
        file.write_all(content.as_bytes()).unwrap();
        path
    }

    #[test]
    fn loads_valid_tsv_library() {
        let path = temp_file(
            "valid.tsv",
            "name\tcomposition\tmonoisotopic_mass\tdiagnostic_ions\tresidue_targets\n\
             HexNAc1\tHexNAc(1)\t203.079373\tHexNAc@204.0867;HexNAc@186.0760[-H2O]\tN\n\
             NeuAcHexNAc\tNeuAc(1)HexNAc(1)\t494.174789\tNeuAc@292.1027\tN;S;T\n",
        );
        let lib = load_glycan_library_file(&path).unwrap();
        assert_eq!(lib.entries.len(), 2);
        assert_eq!(lib.entries[0].composition, "HexNAc(1)");
        assert_eq!(lib.entries[0].diagnostic_ions.len(), 2);
        assert_eq!(lib.entries[0].diagnostic_ions[1].loss_label, "-H2O");
        assert_eq!(lib.entries[1].residue_targets, vec!["N", "S", "T"]);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn loads_valid_csv_library() {
        let path = temp_file(
            "valid.csv",
            "name,composition,monoisotopic_mass,diagnostic_ions,residue_targets\n\
             HexNAc1,HexNAc(1),203.079373,HexNAc@204.0867,N\n",
        );
        let lib = load_glycan_library_file(&path).unwrap();
        assert_eq!(lib.entries.len(), 1);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn rejects_missing_column() {
        let path = temp_file(
            "missing_col.tsv",
            "name\tcomposition\tmonoisotopic_mass\tdiagnostic_ions\n\
             HexNAc1\tHexNAc(1)\t203.079373\tHexNAc@204.0867\n",
        );
        let err = load_glycan_library_file(&path).unwrap_err();
        assert!(err.contains("missing required column 'residue_targets'"));
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn rejects_duplicate_name() {
        let path = temp_file(
            "dup.tsv",
            "name\tcomposition\tmonoisotopic_mass\tdiagnostic_ions\tresidue_targets\n\
             G\tHexNAc(1)\t203.079373\tHexNAc@204.0867\tN\n\
             G\tHexNAc(2)\t406.16\tHexNAc@204.0867\tN\n",
        );
        let err = load_glycan_library_file(&path).unwrap_err();
        assert!(err.contains("duplicate glycan name"));
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn rejects_non_numeric_mass() {
        let path = temp_file(
            "mass.tsv",
            "name\tcomposition\tmonoisotopic_mass\tdiagnostic_ions\tresidue_targets\n\
             G\tHexNAc(1)\tabc\tHexNAc@204.0867\tN\n",
        );
        let err = load_glycan_library_file(&path).unwrap_err();
        assert!(err.contains("non-numeric monoisotopic_mass"));
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn rejects_non_positive_mass() {
        let path = temp_file(
            "negmass.tsv",
            "name\tcomposition\tmonoisotopic_mass\tdiagnostic_ions\tresidue_targets\n\
             G\tHexNAc(1)\t0\tHexNAc@204.0867\tN\n",
        );
        let err = load_glycan_library_file(&path).unwrap_err();
        assert!(err.contains("non-positive monoisotopic_mass"));
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn rejects_empty_diagnostic_ions() {
        let path = temp_file(
            "nodiag.tsv",
            "name\tcomposition\tmonoisotopic_mass\tdiagnostic_ions\tresidue_targets\n\
             G\tHexNAc(1)\t203.079373\t\tN\n",
        );
        let err = load_glycan_library_file(&path).unwrap_err();
        assert!(err.contains("missing 'diagnostic_ions'") || err.contains("empty diagnostic_ions"));
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn rejects_unsupported_delimiter() {
        let path = temp_file(
            "pipe.txt",
            "name|composition|monoisotopic_mass|diagnostic_ions|residue_targets\n",
        );
        let err = load_glycan_library_file(&path).unwrap_err();
        assert!(err.contains("unsupported delimiter"));
        let _ = std::fs::remove_file(path);
    }
}
