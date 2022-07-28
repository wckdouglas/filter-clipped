extern crate clap;
extern crate log;
extern crate rust_htslib;

use clap::Parser;
use std::string::String;

mod workflow;

use workflow::workflow;

/// Remove alignments with high number of clipped base
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// maximum fraction of bases on the sequence being clipped
    /// from either side
    #[clap(short, long, value_parser, default_value_t = 0.1)]
    single_end: f64,

    /// maximum fraction of total bases on the sequence being clipped
    #[clap(short, long, value_parser, default_value_t = 0.1)]
    both_end: f64,

    /// input bam file path  ("-" for stdin)
    #[clap(short, long, value_parser)]
    in_bam: String,

    /// output bam file path ("-" for stdout)
    #[clap(short, long, value_parser, default_value = "-")]
    out_bam: String,

    /// keeping the failed once (high clipped alignments)
    #[clap(long, action)]
    inverse: bool,
}

fn main() {
    env_logger::init();

    let args = Args::parse();

    if args.single_end < 0.0 || args.single_end >= 1.0 {
        panic!("--single-end should be a postive fraction (0 < s <= 1)");
    }

    if args.both_end < 0.0 || args.both_end >= 1.0 {
        panic!("--single-end should be a postive fraction (0 < s <= 1)");
    }

    workflow(
        args.in_bam,
        args.out_bam,
        args.inverse,
        args.both_end,
        args.single_end,
    );
}
