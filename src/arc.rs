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
        let r_squared = self.radius.clone() * self.radius.clone();
        let a_squared = self.alpha.clone() * self.alpha.clone();
        let b_squared = self.beta.clone() * self.beta.clone();
        let multiplier = r_squared / (a_squared + b_squared);
        let normal: Vector2<T> = vector![self.alpha.clone(), self.beta.clone()];
        let unit_normal: Vector2<T> = normal * multiplier;
        Line::from_normal(&unit_normal, self.is_negative)
    }
}
