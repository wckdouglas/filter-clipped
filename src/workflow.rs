use log::info;
use rust_htslib::{bam, bam::Read};

pub fn workflow(
    in_bam: String,
    out_bam: String,
    inverse: bool,
    max_both_end: f64,
    max_single_end: f64,
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

        let keep: bool = total_clip < max_both_end && max_side_clip <= max_single_end;
        if (keep && !inverse) || (inverse && !keep) {
            out_bam.write(&record).unwrap();
            info!("{} {}", total_clip, max_side_clip);
            out_count += 1;
        }
    }
    info!("Read {} alignments", in_count);
    info!("Written {} alignments", out_count);
}

#[cfg(test)]
mod tests {
    use super::*;

    const OUT_BAM: &str = "test/data/out.sam";

    fn count_bam(expected_count: i32) {
        let mut bam_reader = bam::Reader::from_path(&OUT_BAM).unwrap();
        let mut aln_count = 0;
        for r in bam_reader.records() {
            let _record = r.unwrap();
            aln_count += 1;
        }
        assert_eq!(expected_count, aln_count);
    }

    #[test]
    fn test_workflow() {
        workflow(
            "test/data/test.sam".to_string(),
            OUT_BAM.to_string(),
            true,
            0.2,
            0.2,
        );
        count_bam(0);
    }
}
