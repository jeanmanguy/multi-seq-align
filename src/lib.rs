/*! ![stability-experimental](https://img.shields.io/badge/stability-experimental-orange.svg)

A crate to manipulate multiple sequences alignments in Rust.

Instead of storing aligned sequences as multiple strings, `multi_seq_align` stores bases or residues in [`Alignment`] using a list of characters, like a matrix. This allows easy access to specific rows and columns of the alignment.

# Usage

```rust
# use multi_seq_align::Alignment;
# use std::error::Error;
# fn main() -> Result<(), Box<dyn Error>> {
let mut kappa_casein_fragments_alignment = Alignment::with_sequences(
    &[
        b"PAPISKWQSMP".to_vec(),
        b"HAQIPQRQYLP".to_vec(),
        b"PAQILQWQVLS".to_vec(),
    ],
)?;

// Let's extract a column of this alignment
assert_eq!(
    kappa_casein_fragments_alignment.nth_position(6).unwrap(),
    [&b'W', &b'R', &b'W']
);

// But we also have the aligned sequence for the Platypus
// Let's add it to the original alignment
kappa_casein_fragments_alignment.add(
    b"EHQRP--YVLP".to_vec(),
)?;

// the new aligned sequence has a gap at the 6th position
assert_eq!(
    kappa_casein_fragments_alignment.nth_position(6).unwrap(),
    [&b'W', &b'R', &b'W', &b'-']
);

// We can also loop over each position of the alignment
for aas in kappa_casein_fragments_alignment.iter_positions() {
    println!("{:?}", aas);
    assert_eq!(aas.len(), 4); // 4 sequences
}

# Ok(())
# }
```

Here I instancied an alignment using `u8`, but `Alignment` works on generics like numbers, custom or third-party structs.

# Features

- Create [`Alignment`] from one or multiple aligned sequences at once (see [`add()`] and [`with_sequences()`]).
- Extract columns of the alignment (see [`iter_positions()`] and [`iter_sequences(`]).
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
[`iter_positions()`]: struct.Alignment.html#method.iter_positions
[`iter_sequences(`]: struct.Alignment.html#method.iter_sequences
[`add()`]: struct.Alignment.html#method.add
[`with_sequences()`]: struct.Alignment.html#method.with_sequences
*/
#![warn(clippy::all, clippy::pedantic, clippy::nursery, clippy::cargo)]

mod alignment;
mod conservation;
mod errors;
mod iterators;
mod utils;

use displaydoc::Display;

#[cfg(feature = "serde")]
#[macro_use]
extern crate serde;

/// An alignment of DNA or amino acids sequences
///
/// Aligned sequences should all have the same length. Each sequence is stored as one vector of `char`s. This allows an easy access to columns and rows of the alignment.
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)] // Use Rc to implement Copy
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Alignment<T> {
    /// Sequences (as one)
    sequences: Vec<T>,
    /// The number of sequences in the alignment
    n_sequences: usize,
    /// The length of the alignment
    length: usize,
}

struct AlignmentPositionIterator<'a, T> {
    alignment: &'a Alignment<T>,
    index: usize,
    size_hint: usize,
}

struct AlignmentSequenceIterator<'a, T> {
    alignment: &'a Alignment<T>,
    index: usize,
    size_hint: usize,
}

#[derive(Display, Debug)]
pub enum Mutation {
    /// *
    Conserved,
    /// :
    Conservative,
    /// .
    SemiConservative,
    ///
    NonConservative,
}

pub trait Conservation {
    fn positions_are_conserved(&self) -> Vec<Mutation>;
}
