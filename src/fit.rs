use nalgebra::{Vector2, Vector3, vector};

use crate::{
    angle_region::AngleIntervalSet,
    polygon::{Edge, Polygon},
};

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

/// The "critical region" is the region of (N_A)-dimensional space into which the... okay, I don't
/// have a nice name for it, but there's a point in this space that is generated from the polygons,
/// and if that point falls into the critical region then we have a polygon fit.
///
/// The point is fixed if you don't allow the polygons to rotate (i.e. you are checking for a
/// translational fit only), but if you do allow rotation then the point actually follows a really
/// complicated piecewise path through the space.
///
/// The method, then, hinges on being able to rigorously determine when that poorly-behaved path
/// actually dips into the region. (Or, rather, whether or not it ever does - but it's no fun to
/// just be told "yes, this polygon fits inside the other", you want to be able to see!)
///
/// Anyway. The critical region is actually relatively well-behaved - it is bounded by a number of
/// hyperplanes, and it is those hyperplanes that I represent with this `CriticalRegionBoundary`
/// struct.
///
/// You can represent a hyperplane in any dimension with just two values - a normal vector to the
/// plane, and a constant - this mirrors the dot product formulation of a line (because a line is a
/// hyperplane in 2D!), x dot n = c.
///
/// We can simplify our formulation of this, though. Because of how the region is constructed, each
/// of its boundaries actually really only extends in three dimensions, not through the whole space
/// - this is represented by its normal vector having the value zero in all but three of its
/// components.
///
/// Therefore, we only store a 3D normal vector, along with the indices that each of the three
/// components actually are.
///
/// Finally, we don't need to store a constant, because all of these hyperplanes pass through the
/// origin, and therefore the constant is always equal to zero.
struct CriticalRegionBoundary {
    indices: IndexTriple,
    normal: Vector3<f64>,
    threshold: f64,
}

fn critical_region(polygon: &Polygon) -> Vec<CriticalRegionBoundary> {
    let mut result = Vec::new();
    let valid_triples = valid_triples(polygon);

    for (indices, (edge_1, edge_2, edge_3)) in valid_triples {
        let n_x: Vector3<f64> = vector![edge_1.n.x, edge_2.n.x, edge_3.n.x];
        let n_y: Vector3<f64> = vector![edge_1.n.y, edge_2.n.y, edge_3.n.y];
        let cs: Vector3<f64> = vector![edge_1.c, edge_2.c, edge_3.c];

        let n_p: Vector3<f64> = n_x.cross(&n_y);
        let threshold = n_p.dot(&cs);

        result.push(CriticalRegionBoundary {
            indices,
            normal: n_p,
            threshold,
        });
    }

    result
}

/// Represents a fit of the inner polygon inside the outer polygon.
///
/// The translation assumes that both the inner and outer polygon have already been translated such
/// that their centroids lie on the origin.
///
/// The rotation is in radians, and should be applied anticlockwise.
struct Fit {
    translation: Vector2<f64>,
    rotation: f64,
}

fn cross_2d(v_1: &Vector2<f64>, v_2: &Vector2<f64>) -> f64 {
    v_1.x * v_2.y - v_1.y * v_2.x
}

fn find_single_boundary_fit_range(
    // TODO: It's really messy to have to pass in the normals separately like this. It seems like
    // I'm separating out the parts too much, and actually they remain more linked than I am able
    // to handle.
    (n_1, n_2, n_3): (&Vector2<f64>, &Vector2<f64>, &Vector2<f64>),
    boundary: &CriticalRegionBoundary,
    inner: &Polygon,
) -> Option<AngleIntervalSet> {
    let mut result = AngleIntervalSet::whole_circle();

    // Iterate over all possible choices of inner vertex for each of the three indices of the
    // boundary
    let n_b = inner.n_vertices();
    // TODO: Replace nested `for` loops with itertools shortcut
    for j_1 in 0..n_b {
        let b_1: Vector2<f64> = inner
            .vertices()
            .get(j_1)
            .expect("j_1 should never exceed the largest index of vertices of the inner polygon")
            .clone();

        let mut cos_coeff = n_1.dot(&b_1);
        let mut sin_coeff = cross_2d(&n_1, &b_1);

        for j_2 in 0..n_b {
            let b_2: Vector2<f64> = inner
                .vertices()
                .get(j_2)
                .expect(
                    "j_2 should never exceed the largest index of vertices of the inner polygon",
                )
                .clone();

            cos_coeff += n_2.dot(&b_2);
            sin_coeff += cross_2d(&n_2, &b_2);

            for j_3 in 0..n_b {
                let b_3: Vector2<f64> = inner
                    .vertices()
                    .get(j_3)
                    .expect("j_3 should never exceed the largest index of vertices of the inner polygon")
                    .clone();

                cos_coeff += n_3.dot(&b_3);
                sin_coeff += cross_2d(&n_3, &b_3);

                let single_region = AngleIntervalSet::where_lower(
                    cos_coeff,
                    sin_coeff,
                    boundary.threshold,
                )
                .expect(
                    "Should be a non-empty interval as the threshold should be greater than zero",
                );

                if let Some(new_region) = result.intersect(&single_region) {
                    result = new_region;
                } else {
                    return None;
                }
            }
        }
    }

    Some(result)
}

fn find_fit(outer: &Polygon, inner: &Polygon) -> Option<Fit> {
    let critical_region = critical_region(outer);
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
