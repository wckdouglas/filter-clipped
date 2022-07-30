#[derive(Debug)]
pub struct ClipStat {
    left: i64,
    right: i64,
    total_clipped: i64,
}

fn vec_to_max(clip_vec: Vec<i64>) -> i64 {
    let max_clip = clip_vec.iter().max();
    return match max_clip {
        Some(n) => *n,
        _ => 0,
    };
}

fn nbase_to_frac(n_base: i64, seq_len: f64) -> f64 {
    return n_base as f64 / seq_len;
}

impl ClipStat {
    pub fn new(leading_clipped: Vec<i64>, tailing_clipped: Vec<i64>) -> Self {
        let all_clipped = leading_clipped.iter().sum::<i64>() + tailing_clipped.iter().sum::<i64>();

        return Self {
            left: vec_to_max(leading_clipped),
            right: vec_to_max(tailing_clipped),
            total_clipped: all_clipped,
        };
    }

    pub fn right_fraction(&self, seq_len: f64) -> f64 {
        return nbase_to_frac(self.right, seq_len);
    }

    pub fn left_fraction(&self, seq_len: f64) -> f64 {
        return nbase_to_frac(self.left, seq_len);
    }

    pub fn total_fraction(&self, seq_len: f64) -> f64 {
        return nbase_to_frac(self.total_clipped, seq_len);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(vec![2,0], vec![0,2], 0.2, 0.2, 0.4)]
    #[case(vec![1,0], vec![0,2], 0.2, 0.1, 0.3)]
    fn test_clip_stat(
        #[case] leading_clipped: Vec<i64>,
        #[case] trailing_cliped: Vec<i64>,
        #[case] expected_r_frac: f64,
        #[case] expected_l_frac: f64,
        #[case] expected_total_frac: f64,
    ) {
        let seq_len = 10.0;
        let clip_stat = ClipStat::new(leading_clipped, trailing_cliped);
        assert_eq!(expected_r_frac, clip_stat.right_fraction(seq_len));
        assert_eq!(expected_l_frac, clip_stat.left_fraction(seq_len));
        assert_eq!(expected_total_frac, clip_stat.total_fraction(seq_len));
    }

    #[rstest]
    #[case(vec![2,3,0], 3)]
    #[case(vec![1,2,3], 3)]
    #[case(vec![1,0], 1)]
    fn test_vec_to_max(#[case] input_vec: Vec<i64>, #[case] expected_out: i64) {
        assert_eq!(expected_out, vec_to_max(input_vec));
    }

    #[rstest]
    #[case(10, 20.0, 0.5)]
    #[case(10, 40.0, 0.25)]
    #[case(2, 40.0, 0.05)]
    fn test_nbase_to_frac(#[case] n_base: i64, #[case] seq_len: f64, #[case] expected_out: f64) {
        assert_eq!(nbase_to_frac(n_base, seq_len), expected_out);
    }
}
