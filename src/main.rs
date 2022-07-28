extern crate clap;
extern crate log;
extern crate rust_htslib;

use clap::Parser;
use log::info;
use rust_htslib::{bam, bam::Read};
use std::string::String;

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

    let mut out_count = 0;
    let mut in_count = 0;
    let mut in_bam = match args.in_bam.eq("-") {
        true => bam::Reader::from_stdin().unwrap(),
        _ => bam::Reader::from_path(&args.in_bam).unwrap(),
    };
    let header = bam::Header::from_template(in_bam.header());

    let mut out_bam = match args.out_bam.eq("-") {
        true => bam::Writer::from_stdout(&header, bam::Format::Bam).unwrap(),
        _ => bam::Writer::from_path(&args.out_bam, &header, bam::Format::Bam).unwrap(),
    };

    for r in in_bam.records() {
        in_count += 1;
        let record = r.unwrap();
        let seq_len = record.seq().len() as f64;
        let cigar = record.cigar();
        let clipped: Vec<i64> = vec![
            cigar.leading_softclips(),
            cigar.trailing_softclips(),
            cigar.leading_hardclips(),
            cigar.trailing_hardclips(),
        ];
        let max_side_clip_result = clipped.iter().max();
        let max_side_clip: f64 = i64::try_from(*match max_side_clip_result {
            Some(min) => min,
            None => &0,
        })
        .unwrap() as f64
            / seq_len;
        let total_clip: f64 = clipped.iter().sum::<i64>() as f64 / seq_len;

        let keep: bool = total_clip < args.both_end && max_side_clip <= args.single_end;
        if (keep && !args.inverse) || (args.inverse && !keep) {
            out_bam.write(&record).unwrap();
            out_count += 1;
        }
    }
    info!("Read {} alignments", in_count);
    info!("Written {} alignments", out_count);
}
