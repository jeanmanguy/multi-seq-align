# Multi-seq-align

![stability-experimental](https://img.shields.io/badge/stability-experimental-orange.svg)

[![Rust](https://github.com/jeanmanguy/multi-seq-align/workflows/Rust/badge.svg?branch=master)](https://github.com/jeanmanguy/multi-seq-align/actions?query=workflow%3ARust)
[![Rust Documentation](https://img.shields.io/badge/api-rustdoc-blue.svg)](https://docs.rs/multi-seq-align)
[![Crates.io version](https://img.shields.io/crates/v/multi-seq-align)](https://crates.io/crates/multi-seq-align/)
[![Crates.io license](https://img.shields.io/crates/l/multi-seq-align)](https://github.com/jeanmanguy/multi-seq-align/blob/master/LICENSE)


A crate to manipulate multiple sequences alignments in Rust.

Instead of storing aligned sequences as multiple strings, `multi_seq_align` stores bases or residues in `Alignment` using a list of characters, like a matrix. This allows easy access to specific rows and columns of the alignment.


## Usage

```rust
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
```

Here I instancied an alignment using `u8`, but `Alignment` works on generics like numbers, custom or third-party structs.

## Features

- Create `Alignment` from one or multiple aligned sequences at once (see `add()` and `create()`).
- Extract columns of the alignment (see `iter_positions()` and `iter_sequences(`).
This crate is currently in early stage development. I wouldn't recommend using it in production but I am interested in possible ideas to further the developemt of this project. Quite some work needs toi be done to improve the API and make it easy to use in other project.


### Ideas
- Computation of conservation scores
- Identification of conserved sites
- Computation of consensus sequence
- Collapse / trim alignment
- Serialisation / Deserialisation of alignment files
- Extract sub-alignments
    - positions
    - motifs

### Optimisation

My goal is to reduce the footprint of this crate, there is ome work to do to achieve it. The code will eventually be optimised to be faster and to better use memory.

### Issues

Assuring that all the sequences have the same lengths in a generic way is chalenging and result in some not so nice code.


## Ideas & bugs

Please create a new issue on the [project repository](https://github.com/jeanmanguy/multi-seq-align/issues).

## License

Aa-regex is distributed under the terms of the Apache License (Version 2.0). See [LICENSE](./LICENSE) for details.