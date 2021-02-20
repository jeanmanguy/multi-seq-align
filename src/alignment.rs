use crate::{
    errors::MultiSeqAlignError, utils, Alignment, AlignmentPositionIterator,
    AlignmentSequenceIterator,
};

impl<T> Default for Alignment<T>
where
    T: Clone,
{
    fn default() -> Self {
        Self {
            sequences: Vec::<T>::default(),
            n_sequences: 0_usize,
            length: 0_usize,
        }
    }
}

impl<T> Alignment<T> {
    /// Returns the fixed `length` of the Alignment `self`
    #[must_use]
    pub const fn length(&self) -> &usize {
        &self.length
    }

    /// Returns the number of sequences contained in `self`
    #[must_use]
    pub const fn n_sequences(&self) -> &usize {
        &self.n_sequences
    }

    /// Returns an Iterator over the positions of the alignment
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use multi_seq_align::Alignment;
    /// let align = Alignment::with_sequences(
    ///     &[
    ///         b"AVEQTPRK".to_vec(),
    ///         b"SVEQTPRK".to_vec(),
    ///         b"SVEQTPKK".to_vec(),
    ///     ],
    /// )
    /// .unwrap();
    ///
    /// for position in align.iter_positions() {
    ///     assert_eq!(position.len(), 3)
    /// }
    /// ```
    pub fn iter_positions(
        &self,
    ) -> impl Iterator<Item = Vec<&T>> + ExactSizeIterator<Item = Vec<&T>>
    where
        T: Clone,
    {
        AlignmentPositionIterator {
            alignment: self,
            index: 0_usize,
            size_hint: self.length,
        }
    }

    /// Returns an Iterator over the sequences of the alignment
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use multi_seq_align::Alignment;
    /// let align = Alignment::with_sequences(
    ///     &[
    ///         b"AVEQTPRK".to_vec(),
    ///         b"SVEQTPRK".to_vec(),
    ///         b"SVEQTPKK".to_vec(),
    ///     ],
    /// )
    /// .unwrap();
    ///
    /// for sequence in align.iter_sequences() {
    ///     assert_eq!(sequence.len(), 8)
    /// }
    /// ```
    pub fn iter_sequences(
        &self,
    ) -> impl Iterator<Item = Vec<&T>> + ExactSizeIterator<Item = Vec<&T>>
    where
        T: Clone,
    {
        AlignmentSequenceIterator {
            alignment: self,
            index: 0_usize,
            size_hint: self.n_sequences,
        }
    }

    /// Returns an empty `Alignment` of fixed `length`
    ///
    ///
    /// # Examples
    ///
    /// ```rust
    ///  # use multi_seq_align::Alignment;
    ///  let alignment = Alignment::<char>::new(42);
    ///
    ///  assert_eq!(*alignment.length(), 42 as usize);
    ///  assert_eq!(*alignment.n_sequences(), 0 as usize);
    /// ```
    #[must_use]
    pub const fn new(length: usize) -> Self {
        Self {
            sequences: Vec::new(),
            n_sequences: 0_usize,
            length,
        }
    }

    /// Returns `true` if `self` doesn't contains any sequence
    ///
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use multi_seq_align::Alignment;
    /// let alignment = Alignment::<char>::new(42);
    ///
    /// assert!(alignment.is_empty())
    ///```
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.n_sequences == 0_usize
    }

    /// Create an `Alignment` from same length vectors sequences
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use multi_seq_align::Alignment;
    /// let align = Alignment::with_sequences(
    ///     &[
    ///         b"AVEQTPRK".to_vec(),
    ///         b"SVEQTPRK".to_vec(),
    ///         b"SVEQTPKK".to_vec(),
    ///     ],
    /// )
    /// .unwrap();
    ///
    /// assert_eq!(*align.length(), 8);
    /// assert_eq!(*align.n_sequences(), 3);
    /// ```
    ///
    /// # Errors
    ///
    /// Will return an error if the sequences have different lengths (based on the first sequence).
    pub fn with_sequences(sequences: &[Vec<T>]) -> Result<Self, MultiSeqAlignError>
    where
        T: Clone,
    {
        let length = utils::first_sequence_length(sequences);
        utils::check_unequal_lengths(sequences, length)?;

        let n_sequences = sequences.len();

        let sequences_vec = sequences.iter().flat_map(|x| x.to_vec()).collect();

        Ok(Self {
            sequences: sequences_vec,
            n_sequences,
            length,
        })
    }

    /// Add a sequence to `self`
    ///
    /// The new sequence must have the same length than `self.length`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use multi_seq_align::Alignment;
    /// let mut align = Alignment::new(8);
    ///
    ///  assert_eq!(*align.n_sequences(), 0);
    ///
    /// align
    ///     .add(b"AVEQTPRK".to_vec())
    ///     .unwrap();
    ///
    /// assert_eq!(*align.n_sequences(), 1);
    ///
    /// align
    ///     .add(b"SVEQTPRK".to_vec())
    ///     .unwrap();
    ///
    /// assert_eq!(*align.n_sequences(), 2);
    /// ```
    ///
    /// # Errors
    ///
    /// Will return an error if the length of `sequence` is different from the one of the alignment.
    pub fn add<'a>(&'a mut self, sequence: Vec<T>) -> Result<&'a mut Self, MultiSeqAlignError> {
        if sequence.len() != self.length {
            return Err(MultiSeqAlignError::NewSequenceOfDifferentLength {
                expected_length: self.length,
                found_length: sequence.len(),
            });
        }

        self.sequences.extend(sequence);

        self.n_sequences += 1;

        Ok(self)
    }

    /// Returns all amino acids / bases at a `position` in the alignment `self`. The returned vector has a length equal of number of sequences in `self`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use multi_seq_align::Alignment;
    /// let align = Alignment::<u8>::with_sequences(
    ///     &[b"ELK".to_vec(), b"ILK".to_vec()],
    /// )
    /// .unwrap();
    ///
    /// assert_eq!(align.nth_position(0).unwrap(), &[&b'E', &b'I']);
    /// ```
    /// # Panics
    ///
    /// Panics if `n` is greater or equal to the `length` of the Alignment.
    #[must_use]
    pub fn nth_position(&self, n: usize) -> Option<Vec<&T>> {
        assert!(n < self.length);
        (0..self.n_sequences)
            .map(|i| self.sequences.get(i * self.length + n))
            .collect::<Vec<Option<&T>>>()
            .into_iter()
            .collect::<Option<Vec<&T>>>()
    }

    /// Returns all amino acids / bases of the sequence at the `index` of the Alignment `self`. The returned vector has a length equal to the length of the Alignment `self`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use multi_seq_align::Alignment;
    /// let align = Alignment::<u8>::with_sequences(
    ///     &[b"ELK".to_vec(), b"ILK".to_vec()],
    /// )
    /// .unwrap();
    ///
    /// assert_eq!(align.nth_sequence(1).unwrap(), &[&b'I', &b'L', &b'K']);
    /// ```
    /// # Panics
    ///
    /// Panics if `index` is greater or equal to the `n_sequences` of the Alignment.
    #[must_use]
    pub fn nth_sequence(&self, index: usize) -> Option<Vec<&T>> {
        debug_assert!(index < self.n_sequences);

        (0..self.length)
            .map(|i| self.sequences.get(self.length * index + i))
            .collect::<Vec<Option<&T>>>()
            .into_iter()
            .collect::<Option<Vec<&T>>>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::errors::MultiSeqAlignError;
    use pretty_assertions::assert_eq;

    #[test]
    fn new_align() {
        let x = Alignment::<char>::new(5_usize);
        assert!(x.sequences.is_empty());
        assert_eq!(x.length, 5_usize);
        assert_eq!(x.n_sequences, 0_usize);
    }

    #[test]
    fn new_alignment_with_desc() {
        let x = Alignment::<u8>::with_sequences(&[b"ELK".to_vec(), b"ILK".to_vec()]).unwrap();

        assert_eq!(x.sequences, vec![b'E', b'L', b'K', b'I', b'L', b'K']);
        assert_eq!(x.length, 3);
        assert_eq!(x.n_sequences, 2);
    }

    #[test]
    fn add_1_sequence() {
        let mut align =
            Alignment::with_sequences(&[b"ALKHITAN".to_vec(), b"VLK-ITAN".to_vec()]).unwrap();

        align.add(b"ALRYITAT".to_vec()).unwrap();

        assert_eq!(align.n_sequences, 3_usize);
        assert_eq!(align.nth_position(3).unwrap(), vec![&b'H', &b'-', &b'Y'])
    }

    #[test]
    fn add_1_sequence_wrong_length() {
        let mut x = Alignment::new(3_usize);
        let error = x.add(b"ILKAV".to_vec()).err().unwrap();
        let expected = MultiSeqAlignError::NewSequenceOfDifferentLength {
            expected_length: 3_usize,
            found_length: 5_usize,
        };
        assert_eq!(error, expected);
    }

    #[test]
    fn add_to_new() {
        let mut x = Alignment::new(3_usize);

        x.add(b"ELK".to_vec()).unwrap();
        assert_eq!(x.n_sequences, 1_usize);
        assert_eq!(x.length, 3_usize);
        assert_eq!(x.sequences.len(), 3_usize);

        x.add(b"ILK".to_vec()).unwrap();

        assert_eq!(x.n_sequences, 2_usize);
        assert_eq!(x.length, 3_usize);
        assert_eq!(x.sequences.len(), 6_usize);
    }

    #[test]
    fn empty_align() {
        let mut x = Alignment::new(3_usize);
        assert!(x.is_empty());
        x.add(b"ILK".to_vec()).unwrap();
        assert!(!x.is_empty());
    }

    #[test]
    fn nth_residues_3() {
        let align =
            Alignment::with_sequences(&[b"ALKHITAN".to_vec(), b"VLK-ITAN".to_vec()]).unwrap();
        assert_eq!(align.nth_position(3).unwrap(), vec![&b'H', &b'-'])
    }

    #[test]
    fn nth_residues_more_seqs() {
        let align = Alignment::with_sequences(&[
            b"ALKHITAN".to_vec(),
            b"VLK-ITAN".to_vec(),
            b"ALKWITAN".to_vec(),
            b"VLKMITAN".to_vec(),
        ])
        .unwrap();
        assert_eq!(
            align.nth_position(3).unwrap(),
            vec![&b'H', &b'-', &b'W', &b'M']
        )
    }

    #[test]
    #[should_panic(expected = "assertion failed: n < self.length")]
    fn nth_residues_out() {
        let align =
            Alignment::with_sequences(&[b"ALKHITAN".to_vec(), b"VLK-ITAN".to_vec()]).unwrap();
        let _out_of_bonds = align.nth_position(10);
    }

    #[test]
    fn different_seq_lengths() {
        let error = Alignment::with_sequences(&[b"ALKHITAN".to_vec(), b"VLK-ITAN---".to_vec()])
            .err()
            .unwrap();

        let expected = MultiSeqAlignError::MultipleSequencesOfDifferentLengths {
            expected_length: 8,
            found_lengths: vec![11],
        };
        assert_eq!(error, expected);
    }
}
