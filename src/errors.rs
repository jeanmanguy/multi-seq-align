use displaydoc::Display;
use thiserror::Error;

// TODO: manual error handling (https://stevedonovan.github.io/rust-gentle-intro/6-error-handling.html): remove dependencies
#[derive(Debug, Error, Display, PartialEq)]
#[non_exhaustive]
/// Errors
pub enum MultiSeqAlignError {
    /// Expected aligned sequences of length {expected_length}, sequences of lengths: {found_lengths:?}
    MultipleSequencesOfDifferentLengths {
        /// Expected length
        expected_length: usize,
        /// Found lengths
        found_lengths: Vec<usize>,
    },
    /// Expected new aligned sequence of length {expected_length}, found sequence of length {found_length}
    NewSequenceOfDifferentLength {
        /// Expected length
        expected_length: usize,
        /// Found length
        found_length: usize,
    },
}
