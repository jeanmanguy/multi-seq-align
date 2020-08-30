use displaydoc::Display;
use thiserror::Error;

// TODO: manual error handling (https://stevedonovan.github.io/rust-gentle-intro/6-error-handling.html): remove dependencies
#[derive(Debug, Error, Display, PartialEq)]
#[non_exhaustive]
/// Errors
pub enum MultiSeqAlignError {
    /// Expected aligned sequences of length {expected_length}, sequences {sequences_names:?} to have different lengths: {found_lengths:?}
    MultipleSequencesOfDifferentLengths {
        /// Expected length
        expected_length: usize,
        /// Names of the problematic sequences
        sequences_names: Vec<String>,
        /// Found lengths
        found_lengths: Vec<usize>,
    },
    /// Expected new aligned sequence of length {expected_length}, found sequence `{sequences_name}` of length {found_length}
    NewSequenceOfDifferentLength {
        /// Expected length
        expected_length: usize,
        /// Name of the problematic sequence
        sequences_name: String,
        /// Found length
        found_length: usize,
    },
}
