use crate::RayTracing;

use super::{Arithmetic, Ray, Vector};

/// # Conic definition
///
/// A conic surface is defined by the set of coordinates $(x,y,z)$ that satisfies
/// $$ F(x,y,z) = r^2 - 2zR + z^2(\kappa+1)=0,$$
/// where $r^2=x^2+y^2$, $R$ is the radius of curvature and $\kappa$ is the conic constant.
pub struct Conic {
    /// Conic constant $\kappa$
    pub constant: f64,
    /// Radius of curvature
    pub radius: f64,
    /// Origin vector
    pub origin: Vector,
}
impl Conic {
    /// Creates a new `Conic`
    pub fn new(constant: f64, radius: f64) -> Self {
        Self {
            constant,
            radius,
            origin: [0f64; 3],
        }
    }
    pub fn origin(self, origin: Vector) -> Self {
        Self { origin, ..self }
    }
    /// Creates a new `Conic` with GMT M1 prescription:
    ///  - $\kappa=-0.9982857$
    ///  - $R=36$
    pub fn gmt_m1() -> Self {
        Self::new(-0.9982857, 36.0)
    }
    /// Creates a new `Conic` with GMT M2 prescription:
    ///  - $\kappa=-0.71692784$
    ///  - $R=-4.1639009$
    pub fn gmt_m2() -> Self {
        Self {
            constant: -0.71692784,
            radius: -4.1639009,
            origin: [0f64, 0f64, 20.26247614],
        }
    }
    /// Solve conic surface $F(x,y,z)$ for z :
    pub fn height_at(&self, v: Vector) -> Vector {
        let mut _v = v;
        _v[2] = 0f64;
        let r2 = _v.norm_square();
        _v[2] = self.c() * r2 / (1f64 + self.sqrt_(r2));
        _v
    }
    /// Conic $x$ partial derivative
    /// $$
    /// \frac{\partial{F}}{\partial x} = -\frac{x}{R+\sqrt{R^2-(\kappa+1)r^2}}
    /// $$
    ///
    pub fn x_partial_at(&self, v: Vector) -> f64 {
        let mut _v = v;
        _v[2] = 0f64;
        let r2 = _v.norm_square();
        -self.c() * _v[0] / self.sqrt_(r2)
    }
    /// Conic $y$ partial derivative
    /// $$
    /// \frac{\partial{F}}{\partial y} = -\frac{y}{R+\sqrt{R^2-(\kappa+1)r^2}}
    /// $$
    ///
    pub fn y_partial_at(&self, v: Vector) -> f64 {
        let mut _v = v;
        _v[2] = 0f64;
        let r2 = _v.norm_square();
        -self.c() * _v[1] / self.sqrt_(r2)
    }
    /// Conic $z$ partial derivative
    /// $$
    /// \frac{\partial{F}}{\partial z} = 1
    /// $$
    ///
    pub fn z_partial_at(&self, _v: Vector) -> f64 {
        1f64
    }
    /// Normal vector, $\vec n = [\partial_x F,\partial_y F,\partial_z F]$, to conic surface
    pub fn normal_at(&self, v: Vector) -> Vector {
        [
            self.x_partial_at(v),
            self.y_partial_at(v),
            self.z_partial_at(v),
        ]
        .normalize()
    }
    /// Reflect ray from conic surface
    /// The ray direction vector after reflection is given by: $\vec{u^\prime} = \vec u - 2 (\vec u \cdot \vec n)\vec n$
    pub fn reflect(&self, ray: &mut Ray) {
        let n = self.normal_at(ray.p);
        let q = 2f64 * ray.u.dot(&n);
        ray.u[0] -= q * n[0];
        ray.u[1] -= q * n[1];
        ray.u[2] -= q * n[2];
        ray.u = ray.u.normalize();
    }
    fn kp1(&self) -> f64 {
        self.constant + 1f64
    }
    fn c(&self) -> f64 {
        1f64 / self.radius
    }
    fn sqrt_(&self, r2: f64) -> f64 {
        (1f64 - self.kp1() * self.c() * self.c() * r2).sqrt()
    }
}

impl RayTracing for Conic {
    fn distance(&self, ray: &Ray) -> f64 {
        let q = (self.constant + 1f64).sqrt();
        let p: Vector = [
            ray.p[0] - self.origin[0],
            ray.p[1] - self.origin[1],
            ray.p[2] - self.origin[2],
        ];
        let alpha: Vector = [p[0], p[1], p[2] * q];
        let beta: Vector = [ray.u[0], ray.u[1], ray.u[2] * q];
        let a = beta.norm_square();
        let b = 2f64 * (alpha.dot(&beta) - ray.u[2] * self.radius);
        let c = alpha.norm_square() - 2f64 * p[2] * self.radius;
        0.5 * (-b + (b * b - 4f64 * a * c).sqrt()) / a
    }
}
