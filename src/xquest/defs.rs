//! xQuest definition file generation.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::cli::settings::Settings;
use crate::crosslinker::{CrosslinkerLabel, CrosslinkerProfile};
use crate::jobs::GlycanVariant;

#[derive(Debug, Clone, PartialEq)]
pub struct JobDefs {
    pub xquest_def: String,
}

pub fn write_job_defs(
    job_dir: &Path,
    xquest_root: &Path,
    crosslinker: &CrosslinkerProfile,
    settings: &Settings,
    variant: &GlycanVariant,
    fasta_path: &Path,
) -> Result<JobDefs, String> {
    let defs = build_defs(xquest_root, crosslinker, settings, variant, fasta_path)?;
    fs::create_dir_all(job_dir).map_err(|err| {
        format!("cannot create job directory {}: {err}", job_dir.display())
    })?;

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
    variant: &GlycanVariant,
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
        CrosslinkerLabel::LightHeavy => (
            crosslinker.shift_da,
            1,
            0,
            crosslinker.shift_da,
        ),
        CrosslinkerLabel::LightOnly => (0.0, 0, 1, 0.0),
        CrosslinkerLabel::None => (0.0, 0, 0, 0.0),
    };

    let mut variable_mod = format!("N,{:.6}", variant.mass);
    let mut nvariable_mod = 1usize;
    if settings.variable_oxidation {
        variable_mod.push_str(&format!(",M,{:.6}", 15.994915));
        nvariable_mod += 1;
    }

    let aa_required = if crosslinker.nterm_xlinkable || settings.nterm_xlinkable {
        format!("{},K:Z,Z:Z", crosslinker.xlink_sites)
    } else {
        crosslinker.xlink_sites.clone()
    };

    let mut lines: Vec<String> = template
        .lines()
        .map(|line| substitute_template_line(line, &fasta_path, crosslinker, settings, variant))
        .collect();

    set_or_append(&mut lines, "database", &fasta_path.display().to_string());
    set_or_append(&mut lines, "AArequired", &aa_required);
    set_or_append(&mut lines, "xkinkerID", &crosslinker.name.to_uppercase());
    set_or_append(&mut lines, "crosslinkername", &crosslinker.name);
    set_or_append(
        &mut lines,
        "xlinkermw",
        &format!("{:.7}", crosslinker.xlinkermw),
    );
    set_or_append(&mut lines, "ms2tolerance", &format!("{:.4}", settings.ms2_tolerance_da));
    set_or_append(&mut lines, "variable_mod", &variable_mod);
    set_or_append(&mut lines, "nvariable_mod", &nvariable_mod.to_string());
    set_or_append(&mut lines, "outputpath", "results");
    set_or_append(&mut lines, "copydb2resdir", "0");
    set_or_append(&mut lines, "RuntimeDecoys", "0");
    set_or_append(&mut lines, "cp_isotopediff", &format!("{:.6}", cp_isotope_diff));
    set_or_append(&mut lines, "drawspectra", "0");
    set_or_append(&mut lines, "cp_minpeaknumber", "1");

    if crosslinker.nterm_xlinkable || settings.nterm_xlinkable {
        set_or_append(&mut lines, "ntermxlinkable", "1");
    }

    set_or_append(
        &mut lines,
        "isotopeshift",
        &format!("{:.7}", isotope_shift),
    );
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
    set_or_append(&mut lines, "xlinktypes", "0011");

    let xquest_def = lines.join("\n") + "\n";

    Ok(JobDefs { xquest_def })
}

fn substitute_template_line(
    line: &str,
    fasta_path: &Path,
    crosslinker: &CrosslinkerProfile,
    settings: &Settings,
    variant: &GlycanVariant,
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
        "variable_mod" => format!("variable_mod N,{:.6}", variant.mass),
        "nvariable_mod" => {
            let count = 1 + usize::from(settings.variable_oxidation);
            format!("nvariable_mod {count}")
        }
        "outputpath" => "outputpath results".into(),
        "RuntimeDecoys" => "RuntimeDecoys 0".into(),
        _ => line.to_string(),
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
xlinktypes 0011
AArequired K:K
xkinkerID DSS
xlinkermw 138.0680796
Iontagmode 1
RuntimeDecoys 0
ms2tolerance 0.2
cp_isotopediff 12.075321
cp_minpeaknumber 2
outputpath results
crosslinkername DSS
"#
    .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crosslinker::CrosslinkerProfile;

    #[test]
    fn dmtmm_defs_use_zero_isotope_shift() {
        let settings = Settings::defaults();
        let crosslinker = CrosslinkerProfile::resolve(&settings, Some("dmtmm")).unwrap();
        let variant = GlycanVariant {
            glycan_name: "HexNAc(1)".into(),
            composition: "HexNAc(1)".into(),
            mass: 203.079373,
            loss_label: String::new(),
        };
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("xQuest/V2.1.6/xquest");
        let defs = build_defs(
            &root,
            &crosslinker,
            &settings,
            &variant,
            Path::new("proteins.fasta"),
        )
        .unwrap();
        assert!(defs.xquest_def.contains("isotopeshift 0.0000000"));
        assert!(defs.xquest_def.contains("printisotopicscanpairs 0"));
        assert!(defs.xquest_def.contains("AArequired K:E,K:D"));
        assert!(defs.xquest_def.contains("Iontagmode"));
        assert!(defs.xquest_def.contains("outputpath results"));
    }
}
