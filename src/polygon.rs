use geo::{Coord, IsConvex, convex_hull::quick_hull};
use nalgebra::{Vector2, vector};

/// A wrapper for geo::Polygon that enforces convexity, and the vertices being defined
/// anti-clockwise (or counter-clockwise, if you're American).
///
/// I've named it just `Polygon` rather than `CcwConvexPolygon` because it's the only type of
/// polygon I need to care about!
#[derive(Clone)]
pub(crate) struct Polygon(geo::Polygon);

/// An `Edge` is a container for the two things required to define a line in two-dimensional
/// Euclidean space using dot-product form: a vector and a constant. The resultant line can be
/// represented as the set of points x in R^2 that satisfy x dot n = c.
#[derive(Clone)]
pub(crate) struct Edge {
    pub(crate) n: Vector2<f64>,
    pub(crate) _c: f64,
}

impl Polygon {
    pub fn from_vertices(vertices: &[Vector2<f64>]) -> Self {
        let mut points: Vec<Coord<f64>> = vertices
            .iter()
            .map(|v| Coord::<f64> { x: v.x, y: v.y })
            .collect();
        let polygon = geo::Polygon::new(quick_hull(&mut points), vec![]);
        assert!(polygon.exterior().is_ccw_convex());
        Self(polygon)
    }

    pub fn n_vertices(&self) -> usize {
        self.0.exterior().coords().count() - 1
    }

    pub fn _vertices(&self) -> Vec<Vector2<f64>> {
        self.0
            .exterior()
            .coords()
            .map(|&c| vector![c.x, c.y])
            .take(self.n_vertices())
            .collect()
    }

    pub fn edges(&self) -> Vec<Edge> {
        let mut result = Vec::new();
        // We can exploit the fact that the final co-ordinate is stored twice
        for idx in 0..self.n_vertices() {
            let from = self
                .0
                .exterior()
                .coords()
                .nth(idx)
                .expect("Vertex list should have n_vertices + 1 elements");
            let to = self
                .0
                .exterior()
                .coords()
                .nth(idx + 1)
                .expect("Vertex list should have n_vertices + 1 elements");
            let edge: Vector2<f64> = vector![to.x, to.y] - vector![from.x, from.y];
            let normal: Vector2<f64> = vector![-edge.y, edge.x];
            let unit_normal: Vector2<f64> = normal / normal.dot(&normal).sqrt();

            // We can easily calculate the constant c using these two facts:
            // 1. For any point x on the edge, x dot n = c
            // 2. Both of the vertices defining this edge lie on it
            // Therefore, c = n dot (either of the vertices)
            let constant = unit_normal.dot(&vector![from.x, from.y]);
            result.push(Edge {
                n: unit_normal,
                _c: constant,
            });
        }
        result
    }
}

#[cfg(test)]
mod test {
    use nalgebra::vector;

    use crate::polygon::Polygon;

    #[test]
    fn constructor_does_not_change_square() {
        // Arrange
        let in_vertices = vec![
            vector![0.0f64, 0.0f64],
            vector![1.0f64, 0.0f64],
            vector![1.0f64, 1.0f64],
            vector![0.0f64, 1.0f64],
        ];

        // Act
        let poly = Polygon::from_vertices(&in_vertices);
        let out_vertices = poly._vertices();

        // Assert
        assert_eq!(out_vertices.len(), 4);
        assert!(out_vertices.contains(in_vertices.first().expect("Vec should have four values")));
        assert!(out_vertices.contains(in_vertices.get(1).expect("Vec should have four values")));
        assert!(out_vertices.contains(in_vertices.get(2).expect("Vec should have four values")));
        assert!(out_vertices.contains(in_vertices.get(3).expect("Vec should have four values")));
    }

    #[test]
    fn constructor_removes_vertices_in_interior() {
        // Arrange
        let in_vertices = vec![
            vector![0.0f64, 0.0f64],
            vector![1.0f64, 0.0f64],
            vector![1.0f64, 1.0f64],
            vector![0.0f64, 1.0f64],
            vector![0.5f64, 0.5f64],
        ];

        // Act
        let poly = Polygon::from_vertices(&in_vertices);
        let out_vertices = poly._vertices();

        // Assert
        assert_eq!(out_vertices.len(), 4);
        assert!(!out_vertices.contains(in_vertices.get(4).expect("Vec should have five values")));
    }

    #[test]
    fn constructor_orders_vertices_anticlockwise() {
        // Arrange
        let in_vertices = vec![
            vector![1.0f64, 1.0f64],
            vector![0.0f64, 0.0f64],
            vector![0.0f64, 1.0f64],
            vector![1.0f64, 0.0f64],
        ];

        // Act
        let poly = Polygon::from_vertices(&in_vertices);
        let out_vertices = poly._vertices();
        // find index of origin - we don't care which vertex comes first, only that the cyclic
        // order is correct
        let origin_idx = out_vertices
            .iter()
            .position(|&v| v == vector![0.0f64, 0.0f64])
            .expect("Origin should be present in vertex list");

        // Assert
        assert_eq!(
            *out_vertices
                .get(origin_idx)
                .expect("Vec should have four values"),
            vector![0.0f64, 0.0f64]
        );
        assert_eq!(
            *out_vertices
                .get((origin_idx + 1) % 4)
                .expect("Vec should have four values"),
            vector![1.0f64, 0.0f64]
        );
        assert_eq!(
            *out_vertices
                .get((origin_idx + 2) % 4)
                .expect("Vec should have four values"),
            vector![1.0f64, 1.0f64]
        );
        assert_eq!(
            *out_vertices
                .get((origin_idx + 3) % 4)
                .expect("Vec should have four values"),
            vector![0.0f64, 1.0f64]
        );
    }
}
