/*! ![stability-experimental](https://img.shields.io/badge/stability-experimental-orange.svg)

A crate to manipulate multiple sequences alignments in Rust.

Instead of storing aligned sequences as multiple strings, `multi_seq_align` stores bases or residues in [`Alignment`] using a list of characters, like a matrix. This allows easy access to specific rows and columns of the alignment.

# Usage

```rust
# use multi_seq_align::Alignment;
# use std::error::Error;
# fn main() -> Result<(), Box<dyn Error>> {
let mut kappa_casein_fragments_alignment = Alignment::create(
    vec![
        "P06796".to_string(), // Mouse
        "P07498".to_string(), // Human
        "P02668".to_string(), // Cattle
    ],
    vec![
        Some("CASK_MOUSE".to_string()),
        Some("CASK_HUMAN".to_string()),
        Some("CASK_BOVIN".to_string()),
    ],
    &[
        "PAPISKWQSMP".to_string(),
        "HAQIPQRQYLP".to_string(),
        "PAQILQWQVLS".to_string(),
    ],
)?;

// Let's extract a column of this alignment
assert_eq!(
    kappa_casein_fragments_alignment.nth_position(6),
    [Some(&'W'), Some(&'R'), Some(&'W')]
);

// But we also have the aligned sequence for the Platypus
// Let's add it to the original alignment
kappa_casein_fragments_alignment.add_aligned_sequence(
    "D0QJA9".to_string(),
    Some("D0QJA9_ORNAN".to_string()),
    "EHQRP--YVLP",
)?;

// the new aligned sequence has a gap at the 6th position
assert_eq!(
    kappa_casein_fragments_alignment.nth_position(6),
    [Some(&'W'), Some(&'R'), Some(&'W'), Some(&'-')]
);

// We can also loop over each position of the alignment
for aas in kappa_casein_fragments_alignment.iter_positions() {
    println!("{:?}", aas);
    assert_eq!(aas.len(), 4); // 4 sequences
}

Ok(())
# }
```

# Features

- Create [`Alignment`] from one or multiple aligned sequences at once (see [`add_aligned_sequence()`] and [`create()`]).
- Extract columns of the alignment (see [`nth_position()`]).
This crate is currently in early stage development. I wouldn't recommend using it in production but I am interested in possible ideas to further the developemt of this project. Quite some work needs toi be done to improve the API and make it easy to use in other project.
# Ideas
- Computation of conservation scores
- Identification of conserved sites
- Computation of consensus sequence
- Collapse / trim alignment
- Serialisation / Deserialisation of alignment files
- Extract sub-alignments
    - positions
    - motifs

# Optimisation

My goal is to reduce the footprint of this crate, there is ome work to do to achieve it. The code will eventually be optimised to be faster and to better use memory.

[`Alignment`]: struct.Alignment.html
[`nth_position()`]: struct.Alignment.html#method.nth
[`add_aligned_sequence()`]: struct.Alignment.html#method.add_aligned_sequence
[`create()`]: struct.Alignment.html#method.create
*/
#![warn(clippy::all, clippy::pedantic, clippy::nursery, clippy::cargo)]

mod errors;
mod utils;

use errors::MultiSeqAlignError;

/// An alignment of DNA or amino acids sequences
///
/// Aligned sequences should all have the same length. Each sequence is stored as one vector of `char`s. This allows an easy access to columns and rows of the alignment.
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)] // Use Rc to implement Copy
pub struct Alignment {
    /// Sequence names
    names: Vec<String>,
    /// Sequence descriptions
    descriptions: Vec<Option<String>>,
    /// Sequences (as one)
    sequences: Vec<char>,
    /// The number of sequences in the alignment
    n_sequences: usize,
    /// The length of the alignment
    length: usize,
}

impl Default for Alignment {
    fn default() -> Self {
        Self {
            names: Vec::<String>::default(),
            descriptions: Vec::<Option<String>>::default(),
            sequences: Vec::<char>::default(),
            n_sequences: 0_usize,
            length: 0_usize,
        }
    }
}

struct AlignmentPositionIterator<'a> {
    alignment: &'a Alignment,
    index: usize,
}

struct AlignmentRowIterator<'a> {
    alignment: &'a Alignment,
    index: usize,
}

impl<'a> Iterator for AlignmentPositionIterator<'a> {
    type Item = Vec<Option<&'a char>>;
    fn next(&mut self) -> Option<Vec<Option<&'a char>>> {
        if self.index >= self.alignment.length {
            return None;
        }
        let col = self.alignment.nth_position(self.index);
        self.index += 1;
        Some(col)
    }
}

impl<'a> Iterator for AlignmentRowIterator<'a> {
    type Item = Vec<Option<&'a char>>;
    fn next(&mut self) -> Option<Vec<Option<&'a char>>> {
        if self.index >= self.alignment.n_sequences {
            return None;
        }
        let col = self.alignment.nth_sequence(self.index);
        self.index += 1;
        Some(col)
    }
}

impl Alignment {
    /// Returns a vector of sequence names
    #[must_use]
    pub const fn names(&self) -> &Vec<String> {
        &self.names
    }

    /// Returns a vector of sequence descriptions
    #[must_use]
    pub const fn descriptions(&self) -> &Vec<Option<String>> {
        &self.descriptions
    }

    /// Returns the fixed `length` of the Alignment `self`
    #[must_use]
    pub const fn length(&self) -> usize {
        self.length
    }

    /// Returns the number of sequences contained in `self`
    #[must_use]
    pub const fn n_sequences(&self) -> usize {
        self.n_sequences
    }

    /// Returns an Iterator over the positions of the alignment
    pub fn iter_positions(&self) -> impl Iterator<Item = Vec<Option<&char>>> {
        AlignmentPositionIterator {
            alignment: self,
            index: 0_usize,
        }
    }

    /// Returns an Iterator over the sequences of the alignment
    pub fn iter_sequences(&self) -> impl Iterator<Item = Vec<Option<&char>>> {
        AlignmentRowIterator {
            alignment: self,
            index: 0_usize,
        }
    }

    /// Returns an empty `Alignment` of fixed `length`
    ///
    ///
    /// # Examples
    ///
    /// ```rust
    ///  # use multi_seq_align::Alignment;
    ///  let alignment = Alignment::new(42);
    ///
    ///  assert_eq!(alignment.length(), 42);
    ///  assert_eq!(alignment.n_sequences(), 0);
    /// ```
    #[must_use]
    pub const fn new(length: usize) -> Self {
        Self {
            names: Vec::new(),
            descriptions: Vec::new(),
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
    /// let alignment = Alignment::new(42);
    ///
    /// assert!(alignment.is_empty())
    ///```
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.n_sequences == 0_usize
    }

    /// Create an `Alignment` from same length vectors of names, descriptions, sequences
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use multi_seq_align::Alignment;
    /// let align = Alignment::create(
    ///     vec![
    ///         "ELMI001940".to_string(),
    ///         "ELMI001939".to_string(),
    ///         "ELMI001938".to_string(),
    ///     ],
    ///     vec![None, None, None],
    ///     &[
    ///         "AVEQTPRK".to_string(),
    ///         "SVEQTPRK".to_string(),
    ///         "SVEQTPKK".to_string(),
    ///     ],
    /// )
    /// .unwrap();
    ///
    /// assert_eq!(align.length(), 8);
    /// assert_eq!(align.n_sequences(), 3);
    /// ```
    ///
    /// # Errors
    ///
    /// Will return an error if `names`, `descriptions` and `sequences` have different lengths, and also if the sequences have different lengths (based on the first sequence).
    pub fn create(
        names: Vec<String>,
        descriptions: Vec<Option<String>>,
        sequences: &[String],
    ) -> Result<Self, MultiSeqAlignError> {
        assert!(names.len() == descriptions.len() && names.len() == sequences.len());

        let length = utils::first_sequence_length(sequences);

        utils::check_unequal_lengths(sequences, &names, length)?;

        let n_sequences = sequences.len();
        let mut sequences_char: Vec<char> = Vec::with_capacity(n_sequences * length);
        sequences.iter().for_each(|seq| {
            seq.chars().for_each(|c| sequences_char.push(c));
        });

        Ok(Self {
            names,
            descriptions,
            sequences: sequences_char,
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
    // assert_eq!(align.n_sequences(), 0);
    ///
    /// align
    ///     .add_aligned_sequence("ELMI001940".to_string(), None, "AVEQTPRK")
    ///     .unwrap();
    ///
    /// assert_eq!(align.n_sequences(), 1);
    ///
    /// align
    ///     .add_aligned_sequence("ELMI001939".to_string(), None, "SVEQTPRK")
    ///     .unwrap();
    ///
    /// assert_eq!(align.n_sequences(), 2);
    /// ```
    ///
    /// # Errors
    ///
    /// Will return an error if the length of `sequence` is different from the one of the alignment.
    pub fn add_aligned_sequence<'a>(
        &'a mut self,
        name: String,
        description: Option<String>,
        sequence: &str,
    ) -> Result<&'a mut Self, MultiSeqAlignError> {
        if sequence.len() != self.length {
            return Err(MultiSeqAlignError::NewSequenceOfDifferentLength {
                expected_length: self.length,
                sequences_name: name,
                found_length: sequence.len(),
            });
        }

        self.names.push(name);

        self.descriptions.push(description);

        self.sequences.extend(sequence.chars());

        self.n_sequences += 1;

        Ok(self)
    }

    /// Returns all amino acids / bases at a `position` in the alignment `self`. The returned vector has a length equal of number of sequences in `self`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use multi_seq_align::Alignment;
    /// // TODO
    ///
    /// ```
    /// # Panics
    ///
    /// Panics if `n` is greater or equal to the `length` of the Alignment.
    #[must_use]
    pub fn nth_position(&self, n: usize) -> Vec<Option<&char>> {
        assert!(n < self.length);
        (0..self.n_sequences)
            .map(|i| self.sequences.get(i * self.length + n))
            .collect::<Vec<Option<&char>>>()
    }

    /// Returns all amino acids / bases of the sequence at the `index` of the Alignment `self`. The returned vector has a length equal to the length of the Alignment `self`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use multi_seq_align::Alignment;
    /// // TODO
    ///
    /// ```
    /// # Panics
    ///
    /// Panics if `index` is greater or equal to the `n_sequences` of the Alignment.
    #[must_use]
    pub fn nth_sequence(&self, index: usize) -> Vec<Option<&char>> {
        assert!(index < self.n_sequences);
        (0..self.length)
            .map(|i| self.sequences.get(i * self.length + index))
            .collect::<Vec<Option<&char>>>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn new_align() {
        let x = Alignment::new(5_usize);
        assert!(x.names.is_empty());
        assert!(x.sequences.is_empty());
        assert_eq!(x.length, 5_usize);
        assert_eq!(x.n_sequences, 0_usize);
    }

    #[test]
    fn new_alignment_with_desc() {
        let x = Alignment::create(
            vec![String::from("testname1"), String::from("testname2")],
            vec![Some(String::from("desc1")), Some(String::from("desc2"))],
            &[String::from("ELK"), String::from("ILK")],
        )
        .unwrap();
        assert_eq!(
            x.names,
            vec![String::from("testname1"), String::from("testname2")]
        );
        assert_eq!(
            x.descriptions,
            vec![Some(String::from("desc1")), Some(String::from("desc2"))]
        );
        assert_eq!(x.sequences, vec!['E', 'L', 'K', 'I', 'L', 'K']);
        assert_eq!(x.length, 3);
        assert_eq!(x.n_sequences, 2);
    }

    #[test]
    fn add_1_sequence() {
        let mut align = Alignment::create(
            vec![String::from("NAME1"), String::from("NAME2")],
            vec![Some(String::from("desc1")), Some(String::from("desc2"))],
            &[String::from("ALKHITAN"), String::from("VLK-ITAN")],
        )
        .unwrap();

        align
            .add_aligned_sequence(String::from("added1"), None, "ALRYITAT")
            .unwrap();

        assert_eq!(align.n_sequences, 3_usize);
        assert_eq!(
            align.nth_position(3),
            vec![Some(&'H'), Some(&'-'), Some(&'Y')]
        )
    }

    #[test]
    fn add_1_sequence_wrong_length() {
        let mut x = Alignment::new(3_usize);
        let error = x
            .add_aligned_sequence(
                String::from("too_long"),
                Some(String::from(
                    "add sequence of length 5 to an alignment of length 3",
                )),
                "ILKAV",
            )
            .err()
            .unwrap();
        let expected = MultiSeqAlignError::NewSequenceOfDifferentLength {
            expected_length: 3_usize,
            sequences_name: String::from("too_long"),
            found_length: 5_usize,
        };
        assert_eq!(error, expected);
    }

    #[test]
    fn add_to_new() {
        let mut x = Alignment::new(3_usize);

        x.add_aligned_sequence(String::from("sequence1"), None, "ILK")
            .unwrap();
        assert_eq!(x.n_sequences, 1_usize);
        assert_eq!(x.length, 3_usize);
        assert_eq!(x.names.len(), 1_usize);
        assert_eq!(x.descriptions.len(), 1_usize);
        assert_eq!(x.sequences.len(), 3_usize);

        x.add_aligned_sequence(String::from("sequence2"), None, "ELK")
            .unwrap();

        assert_eq!(x.n_sequences, 2_usize);
        assert_eq!(x.length, 3_usize);
        assert_eq!(x.names.len(), 2_usize);
        assert_eq!(x.descriptions.len(), 2_usize);
        assert_eq!(x.sequences.len(), 6_usize);
    }

    #[test]
    fn empty_align() {
        let mut x = Alignment::new(3_usize);
        assert!(x.is_empty());
        x.add_aligned_sequence(String::from("sequence1"), None, "ILK")
            .unwrap();
        assert!(!x.is_empty());
    }

    #[test]
    fn nth_residues_3() {
        let align = Alignment::create(
            vec![String::from("NAME1"), String::from("NAME2")],
            vec![Some(String::from("desc1")), Some(String::from("desc2"))],
            &[String::from("ALKHITAN"), String::from("VLK-ITAN")],
        )
        .unwrap();
        assert_eq!(align.nth_position(3), vec![Some(&'H'), Some(&'-')])
    }

    #[test]
    fn nth_residues_more_seqs() {
        let align = Alignment::create(
            vec![
                "seq1".to_string(),
                "seq2".to_string(),
                "seq3".to_string(),
                "seq4".to_string(),
            ],
            vec![None, None, None, None],
            &[
                "ALKHITAN".to_string(),
                "VLK-ITAN".to_string(),
                "ALKWITAN".to_string(),
                "VLKMITAN".to_string(),
            ],
        )
        .unwrap();
        assert_eq!(
            align.nth_position(3),
            vec![Some(&'H'), Some(&'-'), Some(&'W'), Some(&'M')]
        )
    }

    #[test]
    #[should_panic(expected = "assertion failed: n < self.length")]
    fn nth_residues_out() {
        let align = Alignment::create(
            vec![String::from("NAME1"), String::from("NAME2")],
            vec![Some(String::from("desc1")), Some(String::from("desc2"))],
            &[String::from("ALKHITAN"), String::from("VLK-ITAN")],
        )
        .unwrap();
        let _out_of_bonds = align.nth_position(10);
    }

    #[test]
    fn different_seq_lengths() {
        let error = Alignment::create(
            vec![String::from("NAME1"), String::from("NAME2")],
            vec![Some(String::from("desc1")), Some(String::from("desc2"))],
            &[String::from("ALKHITAN"), String::from("VLK-ITANLLL")],
        )
        .err()
        .unwrap();
        let expected = MultiSeqAlignError::MultipleSequencesOfDifferentLengths {
            expected_length: 8,
            sequences_names: vec![String::from("NAME2")],
            found_lengths: vec![11],
        };
        assert_eq!(error, expected);
    }

    #[test]
    fn for_names() {
        let align = Alignment::create(
            vec![String::from("NAME1"), String::from("NAME2")],
            vec![Some(String::from("desc1")), Some(String::from("desc2"))],
            &[String::from("ALKHITAN"), String::from("VLK-ITAN")],
        )
        .unwrap();

        let mut x: Vec<&String> = Vec::new();

        for name in align.names() {
            x.push(name);
        }

        assert_eq!(x, vec!["NAME1", "NAME2"]);
    }

    #[test]
    fn for_columns() {
        let align = Alignment::create(
            vec![String::from("NAME1"), String::from("NAME2")],
            vec![Some(String::from("desc1")), Some(String::from("desc2"))],
            &[String::from("ALKHITAN"), String::from("VLK-ITAN")],
        )
        .unwrap();

        let mut x = Vec::new();

        for col in align.iter_positions() {
            x.push(col);
        }

        assert_eq!(x.get(0).unwrap(), &[Some(&'A'), Some(&'V')]);
        assert_eq!(x.get(3).unwrap(), &[Some(&'H'), Some(&'-')]);
    }

    #[test]
    #[should_panic]
    fn for_columns_out_of_bonds() {
        let align = Alignment::create(
            vec![String::from("NAME1"), String::from("NAME2")],
            vec![Some(String::from("desc1")), Some(String::from("desc2"))],
            &[String::from("ALKHITAN"), String::from("VLK-ITAN")],
        )
        .unwrap();

        let mut x = Vec::new();

        for col in align.iter_positions() {
            x.push(col);
        }

        let _ = x.get(22).unwrap();
    }

    #[test]
    fn for_sequences() {
        let align = Alignment::create(
            vec![String::from("NAME1"), String::from("NAME2")],
            vec![Some(String::from("desc1")), Some(String::from("desc2"))],
            &[String::from("ALKHITAN"), String::from("VLK-ITAN")],
        )
        .unwrap();

        let mut x = Vec::new();

        for row in align.iter_sequences() {
            assert_eq!(row.len(), 8);
            x.push(row);
        }

        assert_eq!(x.len(), 2)
    }
}
