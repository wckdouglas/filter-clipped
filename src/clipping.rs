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
