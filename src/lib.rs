pub mod clipping;
pub mod cli;

use cli::Parser;
use clipping::ClipStat;

use log::{debug, info};
use rust_htslib::{bam, bam::Read};

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
///     0.1
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
) -> Result<u8, String>{
    let mut out_count = 0;
    let mut in_count = 0;
    info!("Reading from alignment file: {}", in_bam);
    info!("Writing to alignment file: {}", out_bam);
    info!(
        "Thresholds: trailing clipped: {}, leading clipped: {}, total clipped: {}",
        right_side, left_side, both_end
    );
    let mut in_bam = match in_bam.eq("-") {
        true => bam::Reader::from_stdin().map_err(|e| e.to_string())?,
        _ => bam::Reader::from_path(&in_bam).map_err(|e| e.to_string())?,
    };
    let header = bam::Header::from_template(in_bam.header());

    let mut out_bam = match out_bam.eq("-") {
        true => bam::Writer::from_stdout(&header, bam::Format::Bam).map_err(|e| e.to_string())?,
        _ => bam::Writer::from_path(&out_bam, &header, bam::Format::Bam).map_err(|e| e.to_string())?,
    };

    for r in in_bam.records() {
        in_count += 1;
        let record = r.map_err(|e| e.to_string())?;
        let seq_len = record.seq().len() as f64;
        let cigar = record.cigar();

        let leading_clipped: Vec<i64> = vec![cigar.leading_softclips(), cigar.leading_hardclips()];
        let trailing_cliped: Vec<i64> =
            vec![cigar.trailing_softclips(), cigar.trailing_hardclips()];

        let clip_stat: ClipStat = ClipStat::new(leading_clipped, trailing_cliped);

        let keep: bool = clip_stat.total_fraction(seq_len)? < both_end
            && clip_stat.left_fraction(seq_len)? <= left_side
            && clip_stat.right_fraction(seq_len)? <= right_side;

        debug!("{:?} {}", clip_stat, seq_len);
        if (keep && !inverse) || (inverse && !keep) {
            out_bam.write(&record).map_err(|e| e.to_string())?;
            out_count += 1;
        }
    }
    info!(
        "Read {} alignments; Written {} alignments",
        in_count, out_count
    );
    Ok(0) // exit code 0
}

    
/// Just a wrapper function to read command line arguments and pass it to `run`
/// 
pub fn wrapper(){
    let args = cli::Command::parse();
    let result = run(
        args.in_bam, args.out_bam, args.inverse, args.both_end, args.left_side, args.right_side,
    );
    match result{
        Ok(_) => (),
        Err(err) => println!("{}", err),
    };
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
    fn test_run(
        #[case] test_case: usize,
        #[case] max_both_end: f64,
        #[case] max_single_end: f64,
        #[case] inverse: bool,
        #[case] expected_count: i32,
    ) {
        let out_bam: &str = &format!("test/data/out_{}.bam", test_case);
        let result = run(
            "test/data/test.sam".to_string(),
            out_bam.to_string(),
            inverse,
            max_both_end,
            max_single_end,
            max_single_end,
        ).unwrap();
        assert_eq!(result, 0);
        count_bam(out_bam.to_string(), expected_count);
    }
}
