use crate::errors::MultiSeqAlignError;

#[inline]
pub fn first_sequence_length<T>(sequences: &[Vec<T>]) -> usize {
    match sequences.get(0) {
        Some(seq) => seq.len(),
        None => 0_usize,
    }
}

// Returns a tuple of vector of indices and vector of found lengths
#[inline]
pub fn check_unequal_lengths<T>(
    seqs: &[Vec<T>],
    expected: usize,
) -> Result<(), MultiSeqAlignError> {
    let mismatches: (Vec<usize>, Vec<usize>) = seqs
        .iter()
        .enumerate()
        .filter_map(|(index, collection)| {
            if collection.len() == expected {
                None
            } else {
                Some((index, collection.len()))
            }
        })
        .unzip();

    if mismatches.0.is_empty() {
        Ok(())
    } else {
        Err(MultiSeqAlignError::MultipleSequencesOfDifferentLengths {
            expected_length: expected,
            found_lengths: mismatches.1,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn unequal_lengths_1() {
        let error = crate::utils::check_unequal_lengths(&[b"ILK".to_vec(), b"ILKS".to_vec()], 3)
            .err()
            .unwrap();
        let expected = MultiSeqAlignError::MultipleSequencesOfDifferentLengths {
            expected_length: 3,
            found_lengths: vec![4],
        };
        assert_eq!(error, expected);
    }

    #[test]
    fn unequal_lengths_2() {
        let error = crate::utils::check_unequal_lengths(
            &[
                b"ELK".to_vec(),
                b"ILKS".to_vec(),
                b"ILK".to_vec(),
                b"ILKSS".to_vec(),
            ],
            3,
        )
        .err()
        .unwrap();
        let expected = MultiSeqAlignError::MultipleSequencesOfDifferentLengths {
            expected_length: 3,
            found_lengths: vec![4, 5],
        };
        assert_eq!(error, expected);
    }
}
