use nalgebra::Vector2;

use crate::polygon::{Edge, Polygon};

/// Compute the sine of the anticlockwise angle from `from` to `to`.
///
/// By "unchecked" I don't mean it in the normal Rust sense - here I refer to the fact that I'm
/// assuming you're passing me unit normals, and I don't check.
///
/// This method won't panic or throw some exception if that's not true, it'll just quietly return
/// you a useless quantity. Use carefully!
///
/// To compute this, we use a little thing I've been calling the "2D cross product", but I'm sure
/// it has a better name than that.
///
/// We basically assume that the two vectors are actually 3D vectors that both lie in the xy
/// plane, and then take the cross product. The resulting vector will lie purely in the direction
/// of the z-axis, so we can just take the z-component directly to find its length.
///
/// The cross product is found by a cross b = |a||b|sin(angle)n, where n is a unit vector in the
/// direction perpendicular to both vectors (which for us is the z-axis).
///
/// By assuming both input vectors have unit length, the resulting vector's length is precisely
/// equal to the sine of the angle between them, so we can just read it out directly.
fn sine_of_angle_unchecked(from: &Vector2<f64>, to: &Vector2<f64>) -> f64 {
    from.x * to.y - from.y * to.x
}

/// By "validating a triple", we mean in reference to the step in the polygon fit method whereby we
/// look through every distinct set of three edges in the outer polygon, and determine if the inner
/// polygon needs to be checked against it.
///
/// A set of three edges meets this criterion if and only if it "could" form a triangle (as opposed
/// to it forming a non-closed infinite area), and as a shorthand for determining that, we simply
/// check if the angles between subsequent edge normals are always less than (or equal to?) 180
/// degrees.
///
/// I make a number of assumptions here that are not very good coding practice and that I would
/// love to improve upon:
/// 1. The passed vectors are all of unit length (currently, I just need to remember to only pass
///    it unit vectors)
/// 2. The passed vectors are in anticlockwise order (ditto)
fn validate_triple(n_1: &Vector2<f64>, n_2: &Vector2<f64>, n_3: &Vector2<f64>) -> bool {
    // Angle between vectors <= 180 <=> sine of that angle >= 0
    sine_of_angle_unchecked(n_1, n_2) >= 0.0
        && sine_of_angle_unchecked(n_2, n_3) >= 0.0
        && sine_of_angle_unchecked(n_3, n_1) >= 0.0
}

type IndexTriple = (usize, usize, usize);
type EdgeTriple = (Edge, Edge, Edge);

pub(crate) fn valid_triples(polygon: &Polygon) -> Vec<(IndexTriple, EdgeTriple)> {
    let mut result = Vec::new();

    // Use the fact that the edges are defined in anticlockwise order to simplify this algorithm
    let edges = polygon.edges();
    #[cfg(debug_assertions)]
    let mut count: usize = 0;
    for idx_1 in 0..(polygon.n_vertices() - 2) {
        for idx_2 in (idx_1 + 1)..(polygon.n_vertices() - 1) {
            for idx_3 in (idx_2 + 1)..(polygon.n_vertices()) {
                let edge_1 = edges
                    .get(idx_1)
                    .expect("idx_1 should never be greater than the last index of the vector");
                let edge_2 = edges
                    .get(idx_2)
                    .expect("idx_2 should never be greater than the last index of the vector");
                let edge_3 = edges
                    .get(idx_3)
                    .expect("idx_3 should never be greater than the last index of the vector");
                if validate_triple(&edge_1.n, &edge_2.n, &edge_3.n) {
                    result.push((
                        (idx_1, idx_2, idx_3),
                        (edge_1.clone(), edge_2.clone(), edge_3.clone()),
                    ));
                }

                #[cfg(debug_assertions)]
                {
                    count += 1;
                }
            }
        }
    }

    #[cfg(debug_assertions)]
    {
        // We should have iterated through the above n choose 3 times
        let n = polygon.n_vertices();
        assert_eq!(count, (n * (n - 1) * (n - 2)) / 6);
    }

    result
}

#[cfg(test)]
mod test {
    use nalgebra::vector;

    use crate::{fit::valid_triples, polygon::Polygon};

    #[test]
    fn all_triples_of_square_are_valid() {
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
}
