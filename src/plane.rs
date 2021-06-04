use super::{Arithmetic, Ray, RayTracing, Vector};

pub struct Plane {
    pub normal: Vector,
    pub origin: Vector,
}
impl Plane {
    pub fn new(origin: Vector, normal: Vector) -> Self {
        Self { normal, origin }
    }
}
impl RayTracing for Plane {
    fn distance(&self, ray: &Ray) -> f64 {
        self.normal.dot(&self.origin.sub(ray.p)) / self.normal.dot(&ray.u)
    }
}
