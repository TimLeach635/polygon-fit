use nalgebra::{ClosedAddAssign, ClosedMulAssign, Scalar, SimdComplexField, Vector2};
use num_traits::Zero;

/// A line, stored in dot-normal form. `normal` is always kept at unit length.
pub struct Line<T> {
    normal: Vector2<T>,
    constant: T,
}

impl<T: Scalar + Zero + ClosedAddAssign + ClosedMulAssign + SimdComplexField> Line<T> {
    /// Construct a line from a normal vector. The resulting line will be perpendicular to the
    /// provided vector, and pass through the point that you get when you interpret the vector as a
    /// position.
    pub fn from_normal(normal: &Vector2<T>, negative: bool) -> Self {
        let mut unit_normal = normal.normalize();
        if negative {
            unit_normal.neg_mut();
        }

        Self {
            constant: normal.dot(&unit_normal),
            normal: unit_normal,
        }
    }
}
