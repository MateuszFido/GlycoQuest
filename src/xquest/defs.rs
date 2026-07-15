//! xQuest definition file generation.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::cli::settings::Settings;
use crate::crosslinker::{CrosslinkerLabel, CrosslinkerProfile};
use crate::jobs::VarModPlan;

#[derive(Debug, Clone, PartialEq)]
pub struct JobDefs {
    pub xquest_def: String,
}

pub fn write_job_defs(
    job_dir: &Path,
    xquest_root: &Path,
    crosslinker: &CrosslinkerProfile,
    settings: &Settings,
    varmod: &VarModPlan,
    fasta_path: &Path,
) -> Result<JobDefs, String> {
    let defs = build_defs(xquest_root, crosslinker, settings, varmod, fasta_path)?;
    fs::create_dir_all(job_dir)
        .map_err(|err| format!("cannot create job directory {}: {err}", job_dir.display()))?;

    let xquest_path = job_dir.join("xquest.def");
    let mut xquest_file = fs::File::create(&xquest_path).map_err(|err| err.to_string())?;
    xquest_file
        .write_all(defs.xquest_def.as_bytes())
        .map_err(|err| err.to_string())?;

    Ok(defs)
}

fn template_path(xquest_root: &Path) -> PathBuf {
    xquest_root.join("deffiles/xQuest/xquest.def")
}

fn build_defs(
    xquest_root: &Path,
    crosslinker: &CrosslinkerProfile,
    settings: &Settings,
    varmod: &VarModPlan,
    fasta_path: &Path,
) -> Result<JobDefs, String> {
    let fasta_path = fasta_path
        .canonicalize()
        .unwrap_or_else(|_| fasta_path.to_path_buf());
    let template_file = template_path(xquest_root);
    let template = if template_file.is_file() {
        fs::read_to_string(&template_file).map_err(|err| {
            format!(
                "cannot read xQuest template {}: {err}",
                template_file.display()
            )
        })?
    } else {
        fallback_template()
    };

    let (isotope_shift, print_pairs, print_light_only, cp_isotope_diff) = match crosslinker.label {
        CrosslinkerLabel::LightHeavy => (crosslinker.shift_da, 1, 0, crosslinker.shift_da),
        CrosslinkerLabel::LightOnly => (0.0, 0, 1, 0.0),
        CrosslinkerLabel::None => (0.0, 0, 0, 0.0),
    };

    let variable_mod = varmod.variable_mod_value();
    let nvariable_mod = variable_modification_limit(settings);

    let aa_required = if crosslinker.nterm_xlinkable || settings.nterm_xlinkable {
        format!("{},K:Z,Z:Z", crosslinker.xlink_sites)
    } else {
        crosslinker.xlink_sites.clone()
    };

    let mut lines: Vec<String> = template
        .lines()
        .map(|line| substitute_template_line(line, &fasta_path, crosslinker, settings, varmod))
        .collect();

    apply_fixed_mods(&mut lines, settings);

    set_or_append(&mut lines, "database", &fasta_path.display().to_string());
    set_or_append(&mut lines, "AArequired", &aa_required);
    set_or_append(&mut lines, "xkinkerID", &crosslinker.name.to_uppercase());
    set_or_append(&mut lines, "crosslinkername", &crosslinker.name);
    set_or_append(
        &mut lines,
        "xlinkermw",
        &format!("{:.7}", crosslinker.xlinkermw),
    );
    set_or_append(
        &mut lines,
        "ms2tolerance",
        &format!("{:.4}", settings.ms2_tolerance_da),
    );
    set_or_append(&mut lines, "variable_mod", &variable_mod);
    set_or_append(&mut lines, "nvariable_mod", &nvariable_mod.to_string());
    set_or_append(&mut lines, "outputpath", "results");
    // Give every job its own database copy + index under results/db. This keeps
    // concurrently-running jobs from contending on a single shared DB_File index
    // next to the source FASTA, which otherwise serializes parallel execution.
    set_or_append(&mut lines, "copydb2resdir", "1");
    set_or_append(&mut lines, "RuntimeDecoys", "0");
    set_or_append(
        &mut lines,
        "cp_isotopediff",
        &format!("{:.6}", cp_isotope_diff),
    );
    set_or_append(&mut lines, "drawspectra", "0");
    set_or_append(&mut lines, "printionmatches", "1");
    set_or_append(&mut lines, "cp_minpeaknumber", "1");

    if crosslinker.nterm_xlinkable || settings.nterm_xlinkable {
        set_or_append(&mut lines, "ntermxlinkable", "1");
    }

    set_or_append(&mut lines, "isotopeshift", &format!("{:.7}", isotope_shift));
    set_or_append(
        &mut lines,
        "printisotopicscanpairs",
        &print_pairs.to_string(),
    );
    set_or_append(
        &mut lines,
        "printlightonlypairs",
        &print_light_only.to_string(),
    );
    set_or_append(&mut lines, "xlinktypes", "1011");

    let xquest_def = lines.join("\n") + "\n";

    Ok(JobDefs { xquest_def })
}

fn variable_modification_limit(settings: &Settings) -> usize {
    let oxidation_slots = usize::from(settings.variable_oxidation);
    settings.max_glycans_per_peptide as usize + oxidation_slots
}

fn substitute_template_line(
    line: &str,
    fasta_path: &Path,
    crosslinker: &CrosslinkerProfile,
    settings: &Settings,
    varmod: &VarModPlan,
) -> String {
    let trimmed = line.trim_start();
    let key = trimmed.split_whitespace().next().unwrap_or("");
    match key {
        "database" => format!("database {}", fasta_path.display()),
        "AArequired" => format!("AArequired {}", crosslinker.xlink_sites),
        "xkinkerID" => format!("xkinkerID {}", crosslinker.name.to_uppercase()),
        "crosslinkername" => format!("crosslinkername {}", crosslinker.name),
        "xlinkermw" => format!("xlinkermw {:.7}", crosslinker.xlinkermw),
        "ms2tolerance" => format!("ms2tolerance {:.4}", settings.ms2_tolerance_da),
        "variable_mod" => format!("variable_mod {}", varmod.variable_mod_value()),
        "nvariable_mod" => format!("nvariable_mod {}", varmod.nvariable_mod()),
        "outputpath" => "outputpath results".into(),
        "RuntimeDecoys" => "RuntimeDecoys 0".into(),
        _ => line.to_string(),
    }
}

/// Rewrite the fixed-modification block so cysteine carbamidomethylation reflects settings.
///
/// The template's `modifications fixed` block lists one residue per line as
/// `<residue>\t<mass>`. We set the `C` row to `57.02146` (on) or `0` (off).
fn apply_fixed_mods(lines: &mut Vec<String>, settings: &Settings) {
    let carbamidomethyl_mass = if settings.fixed_carbamidomethyl_cys {
        "57.02146"
    } else {
        "0"
    };

    let mut in_fixed_block = false;
    for line in lines.iter_mut() {
        let trimmed = line.trim_start();
        if trimmed.starts_with("modifications fixed") {
            in_fixed_block = true;
            continue;
        }
        if !in_fixed_block {
            continue;
        }
        // The fixed block is a flat residue/mass table; stop at the next section header.
        let mut fields = trimmed.split_whitespace();
        match (fields.next(), fields.next()) {
            (Some("C"), Some(_)) => {
                *line = format!("C\t{carbamidomethyl_mass}");
                break;
            }
            (Some(residue), Some(_)) if residue.len() == 1 => continue,
            _ => break,
        }
    }
}

fn set_or_append(lines: &mut Vec<String>, key: &str, value: &str) {
    let replacement = format!("{key} {value}");
    if let Some(line) = lines
        .iter_mut()
        .find(|line| line.trim_start().starts_with(key))
    {
        if let Some(comment_idx) = line.find('#') {
            let comment = &line[comment_idx..];
            *line = format!("{replacement}\t{comment}");
        } else {
            *line = replacement;
        }
    } else {
        lines.push(replacement);
    }
}

fn fallback_template() -> String {
    r#"digestdef
database proteins.fasta
enzyme_num 1
missed_cleavages 2
requiredmissed_cleavages 0
variable_mod 0
nvariable_mod 1
ionseries 010010
xlinktypes 1011
AArequired K:K
xkinkerID DSS
xlinkermw 138.0680796
Iontagmode 1
RuntimeDecoys 0
ms2tolerance 0.2
cp_isotopediff 12.075321
cp_minpeaknumber 2
outputpath results
printionmatches 1
crosslinkername DSS
"#
    .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crosslinker::CrosslinkerProfile;
    use crate::jobs::{GlycanVariant, build_varmod_plan};

    fn n_glycan_varmod(settings: &Settings) -> VarModPlan {
        let variant = GlycanVariant {
            glycan_name: "HexNAc(1)".into(),
            composition: "HexNAc(1)".into(),
            mass: 203.079373,
            loss_label: String::new(),
            residue_targets: vec!['N'],
        };
        build_varmod_plan(&variant, settings).unwrap()
    }

    #[test]
    fn dmtmm_defs_use_zero_isotope_shift() {
        let settings = Settings::defaults();
        let crosslinker = CrosslinkerProfile::resolve(&settings, Some("dmtmm")).unwrap();
        let varmod = n_glycan_varmod(&settings);
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("V2.1.7/xquest");
        let defs = build_defs(
            &root,
            &crosslinker,
            &settings,
            &varmod,
            Path::new("proteins.fasta"),
        )
        .unwrap();
        assert!(defs.xquest_def.contains("isotopeshift 0.0000000"));
        assert!(defs.xquest_def.contains("printisotopicscanpairs 0"));
        assert!(defs.xquest_def.contains("AArequired K:E,K:D"));
        assert!(defs.xquest_def.contains("Iontagmode"));
        assert!(defs.xquest_def.contains("outputpath results"));
        assert!(defs.xquest_def.contains("variable_mod N,203.079373"));
        assert!(defs.xquest_def.contains("xlinktypes 1011"));
        assert!(defs.xquest_def.contains("printionmatches 1"));
    }

    #[test]
    fn o_glycan_variable_mod_targets_ser_thr() {
        let settings = Settings::defaults();
        let crosslinker = CrosslinkerProfile::resolve(&settings, Some("dss")).unwrap();
        let variant = GlycanVariant {
            glycan_name: "HexNAc(1)".into(),
            composition: "HexNAc(1)".into(),
            mass: 203.079373,
            loss_label: String::new(),
            residue_targets: vec!['S', 'T'],
        };
        let varmod = build_varmod_plan(&variant, &settings).unwrap();
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("V2.1.7/xquest");
        let defs = build_defs(
            &root,
            &crosslinker,
            &settings,
            &varmod,
            Path::new("proteins.fasta"),
        )
        .unwrap();
        assert!(
            defs.xquest_def
                .contains("variable_mod S,203.079373,T,203.079373")
        );
        assert!(defs.xquest_def.contains("nvariable_mod 3"));
    }

    #[test]
    fn fixed_carbamidomethyl_toggle_rewrites_cys_row() {
        let crosslinker = CrosslinkerProfile::resolve(&Settings::defaults(), Some("dss")).unwrap();
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("V2.1.7/xquest");

        let mut on = Settings::defaults();
        on.fixed_carbamidomethyl_cys = true;
        let varmod_on = n_glycan_varmod(&on);
        let defs_on = build_defs(
            &root,
            &crosslinker,
            &on,
            &varmod_on,
            Path::new("proteins.fasta"),
        )
        .unwrap();
        assert!(defs_on.xquest_def.contains("C\t57.02146"));

        let crosslinker = CrosslinkerProfile::resolve(&Settings::defaults(), Some("dss")).unwrap();
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("V2.1.7/xquest");
        let varmod = n_glycan_varmod(&Settings::defaults());
        let defs = build_defs(
            &root,
            &crosslinker,
            &Settings::defaults(),
            &varmod,
            Path::new("proteins.fasta"),
        )
        .unwrap();
        assert!(
            defs.xquest_def.contains("copydb2resdir 1"),
            "per-job DB isolation requires copydb2resdir 1 for parallel xQuest jobs"
        );

        let mut off = Settings::defaults();
        off.fixed_carbamidomethyl_cys = false;
        let varmod_off = n_glycan_varmod(&off);
        let defs_off = build_defs(
            &root,
            &crosslinker,
            &off,
            &varmod_off,
            Path::new("proteins.fasta"),
        )
        .unwrap();
        assert!(defs_off.xquest_def.contains("C\t0"));
    }

    #[test]
    fn default_allows_three_glycans_plus_oxidation_per_peptide() {
        let mut settings = Settings::defaults();
        settings.variable_oxidation = true;
        settings.max_glycans_per_peptide = 3;
        let varmod = n_glycan_varmod(&settings);

        let defs = build_defs(
            Path::new("missing-xquest"),
            &CrosslinkerProfile::resolve(&settings, Some("dss")).unwrap(),
            &settings,
            &varmod,
            Path::new("db.fasta"),
        )
        .unwrap();

        assert!(defs.xquest_def.contains("nvariable_mod 4"));
    }
}
