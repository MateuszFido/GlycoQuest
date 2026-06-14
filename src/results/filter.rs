//! GlycoQuest post-filters for xQuest hits.

use crate::cli::settings::Settings;
use crate::crosslinker::CrosslinkerProfile;
use crate::prefilter::PrefilterResult;
use crate::results::extract::XQuestHit;

#[derive(Debug, Clone, PartialEq)]
pub struct AnnotatedHit {
    pub hit: XQuestHit,
    pub postfilter_status: PostfilterStatus,
    pub soft_score: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PostfilterStatus {
    Pass,
    FailGlycanPseudoResidue,
    FailScore,
    FailPrecursorError,
}

impl PostfilterStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Pass => "pass",
            Self::FailGlycanPseudoResidue => "fail_glycan_pseudo_residue",
            Self::FailScore => "fail_score",
            Self::FailPrecursorError => "fail_precursor_error",
        }
    }
}

pub fn apply_postfilters(
    hits: Vec<XQuestHit>,
    settings: &Settings,
    _crosslinker: &CrosslinkerProfile,
    _prefilter: &PrefilterResult,
) -> Vec<AnnotatedHit> {
    hits.into_iter()
        .map(|hit| annotate_hit(hit, settings))
        .collect()
}

fn annotate_hit(hit: XQuestHit, settings: &Settings) -> AnnotatedHit {
    let glycan_count = hit.seq1.matches('X').count() + hit.seq2.matches('X').count();
    let mut status = PostfilterStatus::Pass;
    let mut soft_score = hit.score;

    if glycan_count != 1 {
        status = PostfilterStatus::FailGlycanPseudoResidue;
        soft_score *= 0.5;
    }
    if hit.score < settings.min_score {
        status = PostfilterStatus::FailScore;
    }
    if hit.precursor_error_ppm.abs() > settings.max_precursor_error_ppm {
        status = PostfilterStatus::FailPrecursorError;
    }

    AnnotatedHit {
        hit,
        postfilter_status: status,
        soft_score,
    }
}

pub fn write_annotated_csv(path: &std::path::Path, rows: &[AnnotatedHit]) -> Result<(), String> {
    let mut lines = vec![
        "spectrum_id\tscore\tseq1\tseq2\tprecursor_error_ppm\txlink_position\tpostfilter_status\tsoft_score"
            .to_string(),
    ];
    for row in rows {
        let hit = &row.hit;
        lines.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}",
            hit.spectrum_id,
            hit.score,
            hit.seq1,
            hit.seq2,
            hit.precursor_error_ppm,
            hit.xlink_position,
            row.postfilter_status.as_str(),
            row.soft_score,
        ));
    }
    std::fs::write(path, lines.join("\n") + "\n").map_err(|err| err.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crosslinker::CrosslinkerProfile;

    #[test]
    fn requires_exactly_one_glycan_pseudo_residue() {
        let settings = Settings::defaults();
        let crosslinker = CrosslinkerProfile::resolve(&settings, Some("dss")).unwrap();
        let prefilter = PrefilterResult {
            filtered: vec![],
            isotope_pairs: vec![],
            rejected: vec![],
            pruning: vec![],
            stats: Default::default(),
        };
        let hits = vec![XQuestHit {
            seq1: "AXCDE".into(),
            seq2: "PEPTIDE".into(),
            ..Default::default()
        }];
        let annotated = apply_postfilters(hits, &settings, &crosslinker, &prefilter);
        assert_eq!(
            annotated[0].postfilter_status,
            PostfilterStatus::Pass
        );
    }
}
