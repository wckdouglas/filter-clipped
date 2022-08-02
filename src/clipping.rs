#[derive(Debug)]
/// An object to store statistics for base clipping on
/// an alignment
pub struct ClipStat {
    /// number of bases from the 5' end being clipped
    left: i64,
    /// number of bases from the 3' end being clipped
    right: i64,
    /// total number of clipped bases on the alignment,
    /// this should be the sum of left and right
    total_clipped: i64,
}

/// Helper function to find the maximum value in a list
///
/// # Arguments
/// * `clip_vec`: a list of integers
///
/// # Return
/// * the maximum value from the list, 0 if it's None
///
/// # Examples
/// ```
/// use filter_clipped::clipping::vec_to_max;
/// let list_of_numbers = vec![0,1,2,3];
/// assert_eq!(3, vec_to_max(list_of_numbers));
/// ```
pub fn vec_to_max(clip_vec: Vec<i64>) -> i64 {
    let max_clip = clip_vec.iter().max();
    return match max_clip {
        Some(n) => *n,
        _ => 0,
    };
}

/// Helper function to calculate a fraction given two numbers
///
/// # Arguments
/// * `n_base`: the numerator in the fraction
/// * `seq_len`: the denominator in the fraction
///
/// # Return
/// * fraction: n_base / seq_len
///
/// # Examples
/// ```
/// use filter_clipped::clipping::nbase_to_frac;
/// assert_eq!(nbase_to_frac(10, 10.0) , 1.0)
/// ```
pub fn nbase_to_frac(n_base: i64, seq_len: f64) -> f64 {
    return n_base as f64 / seq_len;
}

impl ClipStat {
    /// Creat a new ClipStat object for an alignment
    ///
    /// # Arguments
    /// * `leading_clipped`: a list of [number of 5' soft clipped bases, number of 5' hard clipped bases]
    /// * `trailing_clipped`: a list of [number of 3' soft clipped bases, number of 3' hard clipped bases]
    ///
    /// # Return:
    /// A ClipStat object
    ///
    /// # Example
    /// ```
    /// use filter_clipped::clipping::ClipStat;
    /// let clip_stat = ClipStat::new(
    ///     vec![0,1],
    ///     vec![0,2],
    /// );
    /// assert_eq!(clip_stat.left(), 1);
    /// assert_eq!(clip_stat.right(), 2);
    /// assert_eq!(clip_stat.total_clipped(), 3);
    /// ```
    pub fn new(leading_clipped: Vec<i64>, tailing_clipped: Vec<i64>) -> Self {
        let all_clipped = leading_clipped.iter().sum::<i64>() + tailing_clipped.iter().sum::<i64>();

        return Self {
            left: vec_to_max(leading_clipped),
            right: vec_to_max(tailing_clipped),
            total_clipped: all_clipped,
        };
    }

    /// Return the fraction of 3' clipped base relative to the sequence length
    ///
    /// # Argument
    /// * `seq_len`: sequence length of the alignment
    ///
    /// # Return:
    /// * `f64` fraction of 3' clipped base
    ///
    /// # Example
    /// ```
    /// use filter_clipped::clipping::ClipStat;
    /// let clip_stat = ClipStat::new(
    ///     vec![0,1],
    ///     vec![0,2],
    /// );
    /// assert_eq!(clip_stat.right_fraction(10.0), 0.2);
    /// ```
    pub fn right_fraction(&self, seq_len: f64) -> f64 {
        return nbase_to_frac(self.right, seq_len);
    }

    /// Return the fraction of 5' clipped base relative to the sequence length
    ///
    /// # Argument
    /// * `seq_len`: sequence length of the alignment
    ///
    /// # Return:
    /// * `f64` fraction of 5' clipped base
    ///
    /// # Example
    /// ```
    /// use filter_clipped::clipping::ClipStat;
    /// let clip_stat = ClipStat::new(
    ///     vec![0,1],
    ///     vec![0,2],
    /// );
    /// assert_eq!(clip_stat.left_fraction(10.0), 0.1);
    /// ```    
    pub fn left_fraction(&self, seq_len: f64) -> f64 {
        return nbase_to_frac(self.left, seq_len);
    }
    /// Return the fraction of total clipped base relative to the sequence length
    ///
    /// # Argument
    /// * `seq_len`: sequence length of the alignment
    ///
    /// # Return:
    /// * `f64` fraction of 5' clipped base
    ///
    /// # Example
    /// ```
    /// use filter_clipped::clipping::ClipStat;
    /// let clip_stat = ClipStat::new(
    ///     vec![0,1],
    ///     vec![0,2],
    /// );
    /// assert_eq!(clip_stat.total_fraction(10.0), 0.3);
    /// ```    
    pub fn total_fraction(&self, seq_len: f64) -> f64 {
        return nbase_to_frac(self.total_clipped, seq_len);
    }

    pub fn left(&self) -> i64 {
        return self.left;
    }

    pub fn right(&self) -> i64 {
        return self.right;
    }

    pub fn total_clipped(&self) -> i64 {
        return self.total_clipped;
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
