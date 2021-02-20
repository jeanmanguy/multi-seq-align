use std::iter::FromIterator;

use crate::{Alignment, AlignmentPositionIterator, AlignmentSequenceIterator};

impl<'a, T> Iterator for AlignmentPositionIterator<'a, T>
where
    T: Clone,
{
    type Item = Vec<&'a T>;
    fn next(&mut self) -> Option<Vec<&'a T>> {
        if self.index >= self.alignment.length {
            return None;
        }
        match self.alignment.nth_position(self.index) {
            Some(position) => {
                self.index = self.index.saturating_add(1);
                self.size_hint = self.size_hint.saturating_sub(1);
                Some(position)
            }
            None => None,
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        if self.size_hint < usize::max_value() {
            (self.size_hint, Some(self.size_hint))
        } else {
            (usize::max_value(), None)
        }
    }
}

impl<'a, T> ExactSizeIterator for AlignmentPositionIterator<'a, T>
where
    T: Clone,
{
    fn len(&self) -> usize {
        let (lower, upper) = self.size_hint();
        // Note: This assertion is overly defensive, but it checks the invariant
        // guaranteed by the trait. If this trait were rust-internal,
        // we could use debug_assert!; assert_eq! will check all Rust user
        // implementations too.
        assert_eq!(upper, Some(lower));
        lower
    }
}

impl<'a, T> Iterator for AlignmentSequenceIterator<'a, T>
where
    T: Clone,
{
    type Item = Vec<&'a T>;
    fn next(&mut self) -> Option<Vec<&'a T>> {
        if self.index >= self.alignment.n_sequences {
            return None;
        }

        match self.alignment.nth_sequence(self.index) {
            Some(seq) => {
                self.index = self.index.saturating_add(1);
                self.size_hint = self.size_hint.saturating_sub(1);
                Some(seq)
            }
            None => None,
        }
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        self.alignment.nth_sequence(n)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        if self.size_hint < usize::max_value() {
            // ?
            (self.size_hint, Some(self.size_hint))
        } else {
            (usize::max_value(), None)
        }
    }
}

impl<'a, T> ExactSizeIterator for AlignmentSequenceIterator<'a, T>
where
    T: Clone,
{
    fn len(&self) -> usize {
        let (lower, upper) = self.size_hint();
        // Note: This assertion is overly defensive, but it checks the invariant
        // guaranteed by the trait. If this trait were rust-internal,
        // we could use debug_assert!; assert_eq! will check all Rust user
        // implementations too.
        assert_eq!(upper, Some(lower));
        lower
    }
}

impl<A> FromIterator<Vec<A>> for Alignment<A>
where
    A: Clone,
{
    /// # Panics
    ///
    /// Panics if sequences are of different lengths
    fn from_iter<I: IntoIterator<Item = Vec<A>>>(iter: I) -> Self {
        let mut length: Option<usize> = None;
        let mut n_sequences = 0_usize;

        let sequences = iter
            .into_iter()
            .flat_map(|x| {
                if length.is_none() {
                    length = Some(x.len());
                } else if Some(x.len()) != length {
                    panic!("sequences of different lengths");
                }

                n_sequences += 1;
                x.to_vec()
            })
            .collect::<Vec<_>>();

        Self {
            sequences,
            n_sequences,
            length: length.unwrap(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn for_positions() {
        let align =
            Alignment::with_sequences(&[b"ALKHITAN".to_vec(), b"VLK-ITAN".to_vec()]).unwrap();

        let mut x = Vec::new();

        for col in align.iter_positions() {
            x.push(col);
        }

        assert_eq!(x.len(), 8);
        assert_eq!(x.get(0).unwrap(), &[&b'A', &b'V']);
        assert_eq!(x.get(3).unwrap(), &[&b'H', &b'-']);
    }

    #[test]
    #[should_panic]
    fn for_positions_out_of_bonds() {
        let align =
            Alignment::with_sequences(&[b"ALKHITAN".to_vec(), b"VLK-ITAN".to_vec()]).unwrap();
        let mut x = Vec::new();

        for col in align.iter_positions() {
            x.push(col);
        }

        let _ = x.get(22).unwrap();
    }

    #[test]
    fn for_positions_exact() {
        let align =
            Alignment::with_sequences(&[b"ALKHITAN".to_vec(), b"VLK-ITAN".to_vec()]).unwrap();

        assert_eq!(align.iter_positions().len(), 8);
        assert_eq!(align.iter_positions().next().unwrap().len(), 2);
    }

    #[test]
    fn for_sequences() {
        let align =
            Alignment::with_sequences(&[b"ALKHITAN".to_vec(), b"VLK-ITAN".to_vec()]).unwrap();

        let mut x = Vec::new();

        for row in align.iter_sequences() {
            assert_eq!(row.len(), 8);
            x.push(row);
        }

        assert_eq!(x.len(), 2)
    }

    #[test]
    fn for_sequences_exact() {
        let align =
            Alignment::with_sequences(&[b"ALKHITAN".to_vec(), b"VLK-ITAN".to_vec()]).unwrap();

        assert_eq!(align.iter_sequences().len(), 2);
        assert_eq!(align.iter_sequences().next().unwrap().len(), 8);
    }

    #[test]
    fn for_sequences_collect() {
        let align =
            Alignment::with_sequences(&[b"ALKHITAN".to_vec(), b"VLK-ITAN".to_vec()]).unwrap();

        assert_eq!(align.iter_sequences().len(), 2);
        assert_eq!(align.iter_sequences().next().unwrap().len(), 8);
    }
}
