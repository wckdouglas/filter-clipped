pub use clap::Parser;
use std::string::String;

/// Remove alignments with high number of clipped base. Sometimes aligner has very loose scoring methods and write alignments with
/// high abundant of soft/hard-clipped base into alignment BAM files.
/// This program is for filtering these reads out by gating the number of clipped bases
/// in relative to the read sequence length
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
pub struct Command {
    /// maximum fraction of bases on the sequence being clipped
    /// from the left side (5' end)
    #[clap(short, long, value_parser=check_fraction, default_value_t = 0.1)]
    pub left_side: f64,

    /// maximum fraction of bases on the sequence being clipped
    /// from the right side (3' end)
    #[clap(short, long, value_parser=check_fraction, default_value_t = 0.1)]
    pub right_side: f64,

    /// maximum fraction of total bases on the sequence being clipped
    #[clap(short, long, value_parser=check_fraction, default_value_t = 0.1)]
    pub both_end: f64,

    /// input bam file path  ("-" for stdin)
    #[clap(short, long, value_parser)]
    pub in_bam: String,

    /// output bam file path ("-" for stdout)
    #[clap(short, long, value_parser, default_value = "-")]
    pub out_bam: String,

    /// keeping the failed ones (high-clipped-fraction alignments)
    #[clap(long, action)]
    pub inverse: bool,

    /// make the record to unmapped instead of removing it, ignore --inverse flag
    #[clap(short, long, action)]
    pub unalign: bool,
}

/// check if a give value is between 0 and 1
///
/// # Arguments
/// - val: a file path to test for it's existence
///
/// # Returns
/// - Err if no
///
/// # Example
/// ```
/// use filter_clipped::cli::check_fraction;
/// assert_eq!(check_fraction("0.5").unwrap(), 0.5);
/// ```
pub fn check_fraction(val: &str) -> Result<f64, String> {
    let f_val: f64 = val.parse::<f64>().map_err(|e| e.to_string())?;
    if (0.0..=1.0).contains(&f_val) {
        Ok(f_val)
    } else {
        Err(format!("{} is not within 0 and 1", val))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case("1.0", 1.0)]
    #[case("0.0", 0.0)]
    #[case("0.5", 0.5)]
    fn test_check_fraction(#[case] val: &str, #[case] out: f64) {
        assert_eq!(check_fraction(val).unwrap(), out);
    }

    #[rstest]
    #[case("1.1")]
    #[case("-0.1")]
    #[should_panic]
    fn test_check_fraction_panic(#[case] val: &str) {
        check_fraction(val).unwrap();
    }
}
