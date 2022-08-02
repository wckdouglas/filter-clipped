extern crate clap;
extern crate log;
extern crate rust_htslib;

use clap::Parser;
use std::string::String;

use filter_clipped::run;

/// Remove alignments with high number of clipped base. Sometimes aligner has very loose scoring methods and write alignments with
/// high abundant of soft/hard-clipped base into alignment BAM files.
/// This program is for filtering these reads out by gating the number of clipped bases
/// in relative to the read sequence length
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// maximum fraction of bases on the sequence being clipped
    /// from the left side (5' end)
    #[clap(short, long, value_parser, default_value_t = 0.1)]
    left_side: f64,

    /// maximum fraction of bases on the sequence being clipped
    /// from the right side (3' end)
    #[clap(short, long, value_parser, default_value_t = 0.1)]
    right_side: f64,

    /// maximum fraction of total bases on the sequence being clipped
    #[clap(short, long, value_parser, default_value_t = 0.1)]
    both_end: f64,

    /// input bam file path  ("-" for stdin)
    #[clap(short, long, value_parser)]
    in_bam: String,

    /// output bam file path ("-" for stdout)
    #[clap(short, long, value_parser, default_value = "-")]
    out_bam: String,

    /// keeping the failed ones (high-clipped-fraction alignments)
    #[clap(long, action)]
    inverse: bool,
}

fn main() {
    env_logger::init();

    let args = Args::parse();

    if args.left_side < 0.0 || args.left_side > 1.0 {
        panic!("--left-side should be a postive fraction (0 < s <= 1)");
    }

    if args.right_side < 0.0 || args.right_side > 1.0 {
        panic!("--right-side should be a postive fraction (0 < s <= 1)");
    }

    if args.both_end < 0.0 || args.both_end > 1.0 {
        panic!("--single-end should be a postive fraction (0 < s <= 1)");
    }

    run(
        args.in_bam,
        args.out_bam,
        args.inverse,
        args.both_end,
        args.left_side,
        args.right_side,
    );
}
