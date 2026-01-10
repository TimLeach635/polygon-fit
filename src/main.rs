use nalgebra::vector;

use crate::{fit::find_fit_range, polygon::Polygon};

mod arc;
mod fit;
mod line;
mod polygon;

fn main() {
    let outer = Polygon::from_vertices(&[
        vector![0.0f64, -1.0f64],
        vector![1.0f64, 0.0f64],
        vector![0.0f64, 1.0f64],
        vector![-1.0f64, 0.0f64],
    ]);
    let inner = Polygon::from_vertices(&[
        vector![0.0f64, 0.0f64],
        vector![1.0f64, 0.0f64],
        vector![1.0f64, 1.0f64],
        vector![0.0f64, 1.0f64],
    ]);

    let fit_range = find_fit_range(&outer, &inner)
        .expect("The inner square should fit within the outer square");

    // There should be 4 total triples, and for the square, all of them should be considered valid
    assert_eq!(
        fit_range.len(),
        4,
        "There should be four distinct sets of rotations for which the inner square fits within the outer one"
    );
}
