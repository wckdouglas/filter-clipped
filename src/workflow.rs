use crate::clipping::ClipStat;
use log::info;
use rust_htslib::{bam, bam::Read};

/// Workflow to process an input bam file and write the pass-filter alignments
/// into a new bam file
///
/// # Arguments
/// - `in_bam`: input bam file
/// - `out_bam`: output bam file
/// - `inverse`: boolean flag to indicate whether we want to write out the failed-filter alignments only
/// - `max_both_end`: maximum fraction of total clipped bases relative to the read sequence length to consider as pass
/// - `max_single_end`: maximum fraction of of clipped bases on either side relative to the read sequence length to consider as pass
///
/// # Examples
///
/// ```
/// fn count_bam(bam_file: String, expected_count: i32) {
///     let mut bam_reader = bam::Reader::from_path(bam_file).unwrap();
///     let mut aln_count = 0;
///     for r in bam_reader.records() {
///         let _record = r.unwrap();
///         aln_count += 1;
///     }
///     assert_eq!(expected_count, aln_count);
/// }
///
/// let out_bam = "out.sam";
/// workflow(
///     "test/data/test.sam".to_string(),
///     out_bam.to_string(),
///     false,
///     0.2,
///     0.3,
///     );
/// count_bam(out_bam.to_string(), 9);
/// ```
pub fn workflow(
    in_bam: String,
    out_bam: String,
    inverse: bool,
    max_both_end: f64,
    left: f64,
    right: f64,
) {
    let mut out_count = 0;
    let mut in_count = 0;
    let mut in_bam = match in_bam.eq("-") {
        true => bam::Reader::from_stdin().unwrap(),
        _ => bam::Reader::from_path(&in_bam).unwrap(),
    };
    let header = bam::Header::from_template(in_bam.header());

    let mut out_bam = match out_bam.eq("-") {
        true => bam::Writer::from_stdout(&header, bam::Format::Bam).unwrap(),
        _ => bam::Writer::from_path(&out_bam, &header, bam::Format::Bam).unwrap(),
    };

    for r in in_bam.records() {
        in_count += 1;
        let record = r.unwrap();
        let seq_len = record.seq().len() as f64;
        let cigar = record.cigar();

        let leading_clipped: Vec<i64> = vec![cigar.leading_softclips(), cigar.leading_hardclips()];
        let trailing_cliped: Vec<i64> =
            vec![cigar.trailing_softclips(), cigar.trailing_hardclips()];

        let clip_stat: ClipStat = ClipStat::new(leading_clipped, trailing_cliped);

        let keep: bool = clip_stat.total_fraction(seq_len) < max_both_end
            && clip_stat.left_fraction(seq_len) <= left
            && clip_stat.right_fraction(seq_len) <= right;

        info!("{:?} {}", clip_stat, seq_len);
        if (keep && !inverse) || (inverse && !keep) {
            out_bam.write(&record).unwrap();
            out_count += 1;
        }
    }
    info!("Read {} alignments", in_count);
    info!("Written {} alignments", out_count);
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::string::String;

    fn count_bam(bam_file: String, expected_count: i32) {
        let mut bam_reader = bam::Reader::from_path(bam_file).unwrap();
        let mut aln_count = 0;
        for r in bam_reader.records() {
            let _record = r.unwrap();
            aln_count += 1;
        }
        assert_eq!(expected_count, aln_count);
    }

    #[rstest]
    #[case(1, 0.2, 0.3, true, 0)]
    #[case(2, 0.2, 0.3, false, 9)]
    #[case(3, 0.1, 0.1, false, 6)]
    fn test_workflow(
        #[case] test_case: usize,
        #[case] max_both_end: f64,
        #[case] max_single_end: f64,
        #[case] inverse: bool,
        #[case] expected_count: i32,
    ) {
        let out_bam: &str = &format!("test/data/out_{}.bam", test_case);
        workflow(
            "test/data/test.sam".to_string(),
            out_bam.to_string(),
            inverse,
            max_both_end,
            max_single_end,
            max_single_end,
        );
        count_bam(out_bam.to_string(), expected_count);
    }
}
