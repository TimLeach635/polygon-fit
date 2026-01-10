use nalgebra::vector;

use crate::{fit::valid_triples, polygon::Polygon};

mod angle_region;
mod fit;
mod polygon;

fn main() {
    let square = Polygon::from_vertices(&[
        vector![0.0f64, 0.0f64],
        vector![1.0f64, 0.0f64],
        vector![1.0f64, 1.0f64],
        vector![0.0f64, 1.0f64],
    ]);

    let triples = valid_triples(&square);

    // There should be 4 total triples, and for the square, all of them should be considered valid
    assert_eq!(triples.len(), 4);
}
