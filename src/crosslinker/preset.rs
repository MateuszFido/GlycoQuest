//! Bundled crosslinker presets.

use super::CrosslinkerProfile;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrosslinkerLabel {
    LightHeavy,
    LightOnly,
    None,
}

impl CrosslinkerLabel {
    pub fn parse(raw: &str) -> Self {
        match raw.trim().to_ascii_lowercase().as_str() {
            "light-heavy" | "light_heavy" | "heavy" => Self::LightHeavy,
            "light-only" | "light_only" | "light" => Self::LightOnly,
            "none" | "unlabeled" | "off" => Self::None,
            _ => Self::LightHeavy,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::LightHeavy => "light-heavy",
            Self::LightOnly => "light-only",
            Self::None => "none",
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct CrosslinkerPreset {
    pub id: &'static str,
    pub label: CrosslinkerLabel,
    pub shift_da: f64,
    pub xlinkermw: f64,
    pub xlink_sites: &'static str,
    pub nterm_xlinkable: bool,
}

impl CrosslinkerPreset {
    const PRESETS: &'static [CrosslinkerPreset] = &[
        CrosslinkerPreset {
            id: "dss",
            label: CrosslinkerLabel::LightHeavy,
            shift_da: 12.075321,
            xlinkermw: 138.0680796,
            xlink_sites: "K:K",
            nterm_xlinkable: false,
        },
        CrosslinkerPreset {
            id: "dmtmm",
            label: CrosslinkerLabel::None,
            shift_da: 0.0,
            xlinkermw: -18.0109,
            xlink_sites: "K:E,K:D",
            nterm_xlinkable: false,
        },
        // Xie et al., Chem. Sci. 2021, DOI 10.1039/D1SC00814E.
        CrosslinkerPreset {
            id: "nhs-cyclooctyne",
            label: CrosslinkerLabel::None,
            shift_da: 0.0,
            xlinkermw: 205.085126607,
            xlink_sites: "X:K",
            nterm_xlinkable: false,
        },
        // Chen et al., Anal. Chem. 2025, DOI 10.1021/acs.analchem.4c04134.
        CrosslinkerPreset {
            id: "ssbxl",
            label: CrosslinkerLabel::None,
            shift_da: 0.0,
            xlinkermw: 573.179438173,
            xlink_sites: "X:K",
            nterm_xlinkable: false,
        },
        // Post-photorelease PCBXL bridge relative to a NeuAc glycan:
        // C25H24N6O3, 456.190988659 Da.
        CrosslinkerPreset {
            id: "pcbxl",
            label: CrosslinkerLabel::None,
            shift_da: 0.0,
            xlinkermw: 456.190988659,
            xlink_sites: "X:K",
            nterm_xlinkable: false,
        },
    ];

    pub fn by_name(name: &str) -> Option<&'static CrosslinkerPreset> {
        let key = name.trim().to_ascii_lowercase();
        Self::PRESETS
            .iter()
            .find(|p| p.id.eq_ignore_ascii_case(&key))
    }

    pub fn to_profile(self) -> CrosslinkerProfile {
        CrosslinkerProfile {
            name: self.id.to_string(),
            label: self.label,
            shift_da: self.shift_da,
            xlinkermw: self.xlinkermw,
            xlink_sites: self.xlink_sites.to_string(),
            nterm_xlinkable: self.nterm_xlinkable,
        }
    }
}
