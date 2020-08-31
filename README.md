# Multi-seq-align


## Usage

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

```

## Ideas & bugs

Please create a new issue on the [project repository](https://github.com/jeanmanguy/multi-seq-align/issues).

## License

Aa-regex is distributed under the terms of the Apache License (Version 2.0). See [LICENSE](./LICENSE) for details.