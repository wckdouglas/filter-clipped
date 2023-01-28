pub mod cli;
pub mod clipping;

use cli::Parser;
use clipping::ClipStat;

use log::{debug, info};
use rust_htslib::{
    bam,
    bam::{record::CigarStringView, Header, Read, Reader, Record},
};

/// Workflow to process an input bam file and write the pass-filter alignments
/// into a new bam file
///
/// # Arguments
/// - `in_bam`: input bam file
/// - `out_bam`: output bam file
/// - `inverse`: boolean flag to indicate whether we want to write out the failed-filter alignments only
/// - `both_end`: maximum fraction of total clipped bases relative to the read sequence length to consider as pass
/// - `left_side`: maximum fraction of of clipped bases on either side relative to the read sequence length to consider as pass
/// - `right_side`: maximum fraction of of clipped bases on either side relative to the read sequence length to consider as pass
///
/// # Examples
///
/// ```
/// use filter_clipped::run;
/// use rust_htslib::bam;
/// use rust_htslib::bam::Read;
/// fn count_bam(bam_file: String, expected_count: i32) {
///     // helper function to verify result
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
/// run(
///     "test/data/test.sam".to_string(),
///     out_bam.to_string(),
///     false,
///     0.1,
///     0.1,
///     0.1,
///     false,
///     );
/// count_bam(out_bam.to_string(), 6);
/// ```
pub fn run(
    in_bam: String,
    out_bam: String,
    inverse: bool,
    both_end: f64,
    left_side: f64,
    right_side: f64,
    unalign: bool,
) -> Result<u8, String> {
    let mut out_count: u32 = 0;
    let mut in_count: u32 = 0;
    let mut unaligned_count: u32 = 0;
    info!("Reading from alignment file: {}", in_bam);
    info!("Writing to alignment file: {}", out_bam);
    info!(
        "Thresholds: trailing clipped: {}, leading clipped: {}, total clipped: {}",
        right_side, left_side, both_end
    );
    let mut in_bam: Reader = match in_bam.eq("-") {
        true => bam::Reader::from_stdin().map_err(|e| e.to_string())?,
        _ => bam::Reader::from_path(&in_bam).map_err(|e| e.to_string())?,
    };
    let header: Header = bam::Header::from_template(in_bam.header());

    let mut out_bam = match out_bam.eq("-") {
        true => bam::Writer::from_stdout(&header, bam::Format::Bam).map_err(|e| e.to_string())?,
        _ => bam::Writer::from_path(&out_bam, &header, bam::Format::Bam)
            .map_err(|e| e.to_string())?,
    };

    for r in in_bam.records() {
        in_count += 1;
        let mut record: Record = r.map_err(|e| e.to_string())?;
        let seq_len: f64 = record.seq().len() as f64;
        let cigar: CigarStringView = record.cigar();

        let leading_clipped: Vec<i64> = vec![cigar.leading_softclips(), cigar.leading_hardclips()];
        let trailing_cliped: Vec<i64> =
            vec![cigar.trailing_softclips(), cigar.trailing_hardclips()];

        let clip_stat: ClipStat = ClipStat::new(leading_clipped, trailing_cliped);

        let keep: bool = clip_stat.total_fraction(seq_len)? < both_end
            && clip_stat.left_fraction(seq_len)? <= left_side
            && clip_stat.right_fraction(seq_len)? <= right_side;

        debug!("{:?} {}", clip_stat, seq_len);
        if !(unalign) {
            if (keep && !inverse) || (inverse && !keep) {
                out_bam.write(&record).map_err(|e| e.to_string())?;
                out_count += 1;
            }
        } else {
            if keep {
                out_bam.write(&record).map_err(|e| e.to_string())?;
            } else {
                record.set_unmapped();
                record.unset_reverse();
                record.unset_proper_pair();
                record.set_tid(-1);
                record.set_pos(-1);
                out_bam.write(&record).map_err(|e| e.to_string())?;
                unaligned_count += 1
            }
            out_count += 1;
        }
    }
    info!(
        "Read {} alignments; Written {} alignments; Making {} to unaligned",
        in_count, out_count, unaligned_count,
    );
    Ok(0) // exit code 0
}

/// Just a wrapper function to read command line arguments and pass it to `run`
///
pub fn wrapper() {
    let args = cli::Command::parse();
    let result = run(
        args.in_bam,
        args.out_bam,
        args.inverse,
        args.both_end,
        args.left_side,
        args.right_side,
        args.unalign,
    );
    match result {
        Ok(_) => (),
        Err(err) => println!("{}", err),
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::string::String;

    fn count_bam(bam_file: String, expected_count: i32, expected_unaligned: i32) {
        let mut bam_reader = bam::Reader::from_path(bam_file).unwrap();
        let mut aln_count = 0;
        let mut unaligned_count = 0;
        for r in bam_reader.records() {
            let _record = r.unwrap();
            aln_count += 1;
            if _record.is_unmapped() {
                unaligned_count += 1;
            }
        }
        assert_eq!(expected_count, aln_count);

        if expected_unaligned > 0 {
            assert_eq!(expected_unaligned, unaligned_count);
        }
    }

    #[rstest]
    #[case(1, 0.2, 0.3, true, 0, 0, false)]
    #[case(2, 0.2, 0.3, false, 9, 0, false)]
    #[case(3, 0.1, 0.1, false, 6, 0, false)]
    #[case(4, 0.1, 0.1, false, 9, 3, true)]
    #[case(4, 0.1, 0.1, false, 9, 3, true)]
    fn test_run(
        #[case] test_case: usize,
        #[case] max_both_end: f64,
        #[case] max_single_end: f64,
        #[case] inverse: bool,
        #[case] expected_count: i32,
        #[case] expected_unaligned: i32,
        #[case] unalign: bool,
    ) {
        let out_bam: &str = &format!("test/data/out_{}.bam", test_case);
        let result = run(
            "test/data/test.sam".to_string(),
            out_bam.to_string(),
            inverse,
            max_both_end,
            max_single_end,
            max_single_end,
            unalign,
        )
        .unwrap();
        assert_eq!(result, 0);
        count_bam(out_bam.to_string(), expected_count, expected_unaligned);
    }
}
