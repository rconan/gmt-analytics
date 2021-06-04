use super::{Arithmetic, Conic, Vector};
use std::fmt;

/// # Ray definition
///
/// A ray is defined with:
///  - a point of origin: $\vec p = [x,y,z]$,
///  - a direction vector: $\vec u = [k,l,m]$ such as $\| \vec u \|=1$.
///
/// The ray tracing equation is given by: $$\vec{p^\prime} = \vec p + s \vec u,$$ where $s$ is the optical path length.
pub struct Ray {
    /// Ray point of origin
    pub p: Vector,
    /// Ray direction vector
    pub u: Vector,
}
/// # Ray builder
///
/// Build a new [`Ray`](crate::analytic::Ray)
pub struct NewRay {
    /// Ray point of origin
    pub p: Vector,
    /// Ray direction vector
    pub u: Vector,
}
impl Default for NewRay {
    fn default() -> Self {
        Self {
            p: [0f64; 3],
            u: [0f64, 0f64, -1f64],
        }
    }
}
impl NewRay {
    /// Build the [`Ray`](crate::analytic::Ray)
    pub fn build(self) -> Ray {
        Ray {
            p: self.p,
            u: self.u,
        }
    }
    /// Set the [`Ray`](crate::analytic::Ray) point of origin
    pub fn point_of_origin(self, p: Vector) -> Self {
        let u = self.u;
        let s = (p[2] / u[2]).abs();
        Self {
            p: [p[0] - s * u[0], p[1] - s * u[1], p[2]],
            ..self
        }
    }
    /// Set the [`Ray`](crate::analytic::Ray) direction vector
    pub fn direction_vector(self, u: Vector) -> Self {
        let p = self.p;
        let s = (p[2] / u[2]).abs();
        Self {
            p: [p[0] - s * u[0], p[1] - s * u[1], p[2]],
            u,
            ..self
        }
    }
    /// Set the [`Ray`](crate::analytic::Ray) direction vector from polar coordinates
    pub fn polar_direction_vector(self, z: f64, a: f64) -> Self {
        let ca = a.cos();
        let sa = a.sin();
        let sz = z.sin();
        let cz = z.cos();
        let u = [sz * ca, sz * sa, -cz].normalize();
        self.direction_vector(u)
    }
}
/// Create a [`NewRay`](crate::analytic::NewRay) at the origin propagate downward (z<0)
pub fn new_ray() -> NewRay {
    NewRay::default()
}
impl Ray {
    /// Compute the distance $s$ from the ray current location to [`Conic`](crate::analytic::Conic)
    /// We find the distance $s$ from:
    /// $$s=\frac{mR-\vec\alpha\cdot\vec\beta + \sqrt{(\vec\alpha\cdot\vec\beta-mR)^2-\|\vec\beta\|^2(\|\vec\alpha\|^2-2zR)}}{\|\vec\beta\|^2}$$
    /// $$\vec\alpha = [x,y,z\sqrt{\kappa+1}]$$
    /// $$\vec\beta = [k,l,m\sqrt{\kappa+1}]$$
    pub fn distance_to(&self, conic: &Conic) -> f64 {
        let q = (conic.constant + 1f64).sqrt();
        let p: Vector = [
            self.p[0] - conic.origin[0],
            self.p[1] - conic.origin[1],
            self.p[2] - conic.origin[2],
        ];
        let alpha: Vector = [p[0], p[1], p[2] * q];
        let beta: Vector = [self.u[0], self.u[1], self.u[2] * q];
        let a = beta.norm_square();
        let b = 2f64 * (alpha.dot(&beta) - self.u[2] * conic.radius);
        let c = alpha.norm_square() - 2f64 * p[2] * conic.radius;
        0.5 * (-b + (b * b - 4f64 * a * c).sqrt()) / a
    }
    pub fn distance_to_with<F>(&self, conic: &Conic, predicate: F) -> f64
    where
        F: Fn(f64, f64) -> f64,
    {
        let q = (conic.constant + 1f64).sqrt();
        let p: Vector = [
            self.p[0] - conic.origin[0],
            self.p[1] - conic.origin[1],
            self.p[2] - conic.origin[2],
        ];
        let alpha: Vector = [p[0], p[1], p[2] * q];
        let beta: Vector = [self.u[0], self.u[1], self.u[2] * q];
        let a = beta.norm_square();
        let b = 2f64 * (alpha.dot(&beta) - self.u[2] * conic.radius);
        let c = alpha.norm_square() - 2f64 * p[2] * conic.radius;
        let (x1, x2) = (
            0.5 * (-b + (b * b - 4f64 * a * c).sqrt()) / a,
            0.5 * (-b - (b * b - 4f64 * a * c).sqrt()) / a,
        );
        predicate(x1, x2)
    }
    /// Trace ray from ray current position to [`Conic`](crate::analytic::Conic)
    pub fn trace_to(&mut self, conic: &Conic) {
        let s = self.distance_to(conic);
        self.p[0] += self.u[0] * s;
        self.p[1] += self.u[1] * s;
        self.p[2] += self.u[2] * s;
    }
    pub fn trace_to_with<F>(&mut self, conic: &Conic, predicate: F) -> &mut Self
    where
        F: Fn(f64, f64) -> f64,
    {
        let s = self.distance_to_with(conic, predicate);
        self.p[0] += self.u[0] * s;
        self.p[1] += self.u[1] * s;
        self.p[2] += self.u[2] * s;
        self
    }
    pub fn trace(&mut self, s: f64) {
        self.p[0] += self.u[0] * s;
        self.p[1] += self.u[1] * s;
        self.p[2] += self.u[2] * s;
    }
    /// Solve ray tracing equation for $z$ given $x$ and $y$
    pub fn solve_for_z(&self, x: f64, y: f64) -> f64 {
        let x = x - self.p[0];
        let y = y - self.p[1];
        let num = x * x + y * y;
        let denom = self.u[0] * self.u[0] + self.u[1] * self.u[1];
        if denom < 1e-30 {
            std::f64::INFINITY
        } else {
            self.p[2] + self.u[2] * (num / denom).sqrt()
        }
    }
}
impl fmt::Display for Ray {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "P: [{:+15.9},{:+15.9},{:+15.9}] ; U: [{:+.9},{:+.9},{:+.9}]",
            self.p[0], self.p[1], self.p[2], self.u[0], self.u[1], self.u[2],
        )
    }
}
pub trait Trace {
    fn trace_to(&mut self, conic: &Conic) -> &mut Self;
    fn trace_to_with<F>(&mut self, conic: &Conic, predicate: F) -> &mut Self
    where
        F: Fn(f64, f64) -> f64 + Copy;
    fn trace(&mut self, s: f64) -> &mut Self;
    fn coordinates(&self) -> Vec<Vector>;
    fn reflect(&mut self, conic: &Conic) -> &mut Self;
}
impl Trace for Vec<Ray> {
    fn trace_to(&mut self, conic: &Conic) -> &mut Self {
        self.iter_mut().for_each(|ray| {
            ray.trace_to(conic);
        });
        self
    }
    fn trace_to_with<F>(&mut self, conic: &Conic, predicate: F) -> &mut Self
    where
        F: Fn(f64, f64) -> f64 + Copy,
    {
        self.iter_mut().for_each(|ray| {
            ray.trace_to_with(conic, predicate);
        });
        self
    }
    fn trace(&mut self, s: f64) -> &mut Self {
        self.iter_mut().for_each(|ray| {
            ray.trace(s);
        });
        self
    }
    fn reflect(&mut self, conic: &Conic) -> &mut Self {
        self.iter_mut().for_each(|ray| conic.reflect(ray));
        self
    }
    fn coordinates(&self) -> Vec<Vector> {
        self.iter().map(|ray| ray.p).collect()
    }
}
