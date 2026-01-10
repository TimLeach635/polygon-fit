use std::f64::consts::PI;

use rust_intervals::{Interval, IntervalSet, interval};

#[derive(Clone, Copy)]
pub(crate) struct AngleInterval {
    from: f64,
    to: f64,
}

impl AngleInterval {
    pub(crate) fn whole_circle() -> Self {
        Self {
            from: 0.0,
            to: 2.0 * PI,
        }
    }

    /// Creates an `AngleInterval` from two `f64` values.
    ///
    /// The inputs are normalised to lie between 0 and 2pi.
    pub(crate) fn from_normalize(from: f64, to: f64) -> Self {
        if from >= to {
            todo!("Handle this case elegantly, probably by returning an `Option`");
        }

        if 2.0 * PI <= to - from {
            return Self::whole_circle();
        }

        Self {
            from: from.rem_euclid(2.0 * PI),
            to: to.rem_euclid(2.0 * PI),
        }
    }

    pub(crate) fn from_real_interval(interval: Interval<f64>) -> Option<Self> {
        if let Some(from) = interval.lower()
            && let Some(to) = interval.upper()
        {
            Some(Self::from_normalize(*from, *to))
        } else {
            None
        }
    }

    /// Find the interval, if one exists, of values of theta for which `cos_coeff * cos(theta) +
    /// sin_coeff * sin(theta) < threshold`.
    ///
    /// Taking everything mod 2*pi, the lower bound is `arctan(cos_coeff/sin_coeff) +
    /// arccos(threshold/sqrt(cos_coeff^2 + sin_coeff^2))`, and the upper bound is the same but
    /// subtracting the arccos.
    pub(crate) fn where_lower(cos_coeff: f64, sin_coeff: f64, threshold: f64) -> Option<Self> {
        if threshold <= 0.0 {
            return None;
        }

        let rms = (cos_coeff * cos_coeff + sin_coeff * sin_coeff).sqrt();
        if threshold >= rms {
            return Some(Self::whole_circle());
        }

        // TODO: Triple-check the signs for this operation
        let center = cos_coeff.atan2(sin_coeff);
        let offset = (threshold / rms).acos();
        debug_assert!(
            offset < PI / 2.0,
            "The offset should be between 0 and pi/2, as the input should be between 0 and 1"
        );

        Some(Self::from_normalize(center + offset, center - offset))
    }

    /// Returns true if the interval "loops around" the circle, which is equivalent to "to" being
    /// less than "from".
    pub(crate) fn does_loop(&self) -> bool {
        self.to < self.from
    }

    /// Returns the angle found right in the centre of this interval, respecting the fact that it's
    /// cyclic.
    ///
    /// Will always be in the range 0 to 2pi.
    pub(crate) fn center(&self) -> f64 {
        let result = if self.does_loop() {
            let mut result = (self.from + self.to) / 2.0 - PI;
            if result < 0.0 {
                result += 2.0 * PI;
            }
            result
        } else {
            (self.from + self.to) / 2.0
        };

        debug_assert!(
            0.0 <= result && result < 2.0 * PI,
            "The center should always be provided in the interval [0, 2pi)"
        );
        result
    }

    pub(crate) fn length(&self) -> f64 {
        let mut to_transformed = self.to;
        if to_transformed < self.from {
            to_transformed += 2.0 * PI;
        }
        to_transformed - self.from
    }

    pub(crate) fn radius(&self) -> f64 {
        self.length() / 2.0
    }

    /// Returns a representation of this `AngleInterval` as if we "unrolled" the circle into the
    /// real line, starting at 0 and ending at, potentially, 4pi.
    pub(crate) fn real_interval(&self) -> Interval<f64> {
        if self.does_loop() {
            interval!(self.from, self.to + 2.0 * PI, "()")
        } else {
            interval!(self.from, self.to, "()")
        }
    }

    pub(crate) fn union(&self, other: &AngleInterval) -> AngleIntervalSet {
        // We can't simply use the interval library as it doesn't know we're on a circle.
        // However, we can convert to the real line, perform the union, and convert back.
        let mut real_result = IntervalSet::<f64>::empty_joining();
        real_result.add(self.real_interval());
        real_result.add(other.real_interval());

        if real_result.len() == 1 {}
    }
}

pub(crate) struct AngleIntervalSet {
    intervals: Vec<AngleInterval>,
}

impl AngleIntervalSet {
    pub(crate) fn whole_circle() -> Self {
        Self {
            intervals: vec![AngleInterval::whole_circle()],
        }
    }

    pub(crate) fn from_single_interval(from: f64, to: f64) -> Self {
        Self {
            intervals: vec![AngleInterval::from_normalize(from, to)],
        }
    }

    /// Returns true if the interval "loops around" the circle, which is equivalent to "to" being
    /// less than "from".
    ///
    /// For multiple intervals, this is true if and only if it is true for any of the member
    /// intervals.
    ///
    /// You can also think of this as being equivalent to zero being a member of the set.
    pub(crate) fn does_loop(&self) -> bool {
        self.intervals.iter().any(|i| i.does_loop())
    }

    /// Mutably add a new interval to this union.
    pub(crate) fn add_interval_mut(&mut self, interval: &AngleInterval) {}

    pub(crate) fn intersect(&self, other: &AngleIntervalSet) -> Option<AngleIntervalSet> {}
}

#[cfg(test)]
mod tests {}
