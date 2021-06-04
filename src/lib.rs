//!
//! # Analytic Ray Tracing
//!

pub mod conic;
pub mod gmt;
pub mod plane;
pub mod ray;
pub use conic::Conic;
pub use gmt::Gmt;
pub use plane::Plane;
pub use ray::{new_ray, Ray, Trace};

pub type Vector = [f64; 3];

pub trait Arithmetic {
    fn dot(&self, other: &[f64]) -> f64;
    fn norm_square(&self) -> f64;
    fn norm(&self) -> f64;
    fn normalize(&mut self) -> Self;
    fn add(&self, other: Self) -> Self;
    fn sub(&self, other: Self) -> Self;
}
impl Arithmetic for Vector {
    fn dot(&self, other: &[f64]) -> f64 {
        self[0] * other[0] + self[1] * other[1] + self[2] * other[2]
    }
    fn norm_square(&self) -> f64 {
        self.dot(self)
    }
    fn norm(&self) -> f64 {
        self.norm_square().sqrt()
    }
    fn normalize(&mut self) -> Self {
        let n = self.norm();
        self[0] /= n;
        self[1] /= n;
        self[2] /= n;
        *self
    }
    fn add(&self, other: Self) -> Self {
        [self[0] + other[0], self[1] + other[1], self[2] + other[2]]
    }
    fn sub(&self, other: Self) -> Self {
        [self[0] - other[0], self[1] - other[1], self[2] - other[2]]
    }
}
pub trait RayTracing {
    fn distance(&self, ray: &Ray) -> f64;
    fn distances(&self, rays: &[Ray]) -> Vec<f64> {
        rays.iter().map(|ray| self.distance(ray)).collect()
    }
    fn traces(&self, rays: &mut [Ray]) -> &Self {
        rays.iter_mut()
            .for_each(|ray| ray.trace(self.distance(ray)));
        self
    }
}
