use nalgebra::{ClosedAddAssign, ClosedMulAssign, Scalar, SimdComplexField, Vector2, vector};
use num_traits::Zero;

use crate::line::Line;

/// Represents an arc of a circle centred at the origin.
///
/// TODO: Write this up, because it's actually pretty complicated.
pub struct Arc<T> {
    radius: T,
    alpha: T,
    beta: T,
    is_negative: bool,
}

impl<T: Scalar + Zero + ClosedAddAssign + ClosedMulAssign + SimdComplexField> Arc<T> {
    pub fn line(&self) -> Line<T> {
        let r_squared = self.radius * self.radius;
        let a_squared = self.alpha * self.alpha;
        let b_squared = self.beta * self.beta;
        let multiplier = r_squared / (a_squared + b_squared);
        Line::from_normal(
            &(multiplier * vector![self.alpha, self.beta]),
            self.is_negative,
        )
    }
}
