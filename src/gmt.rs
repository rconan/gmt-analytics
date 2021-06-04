use super::{new_ray, Arithmetic, Conic, Ray, Vector};
use nalgebra as na;
use skyangle::SkyAngle;
use triangle_rs::{Builder, Delaunay};

pub struct Gmt {
    source_zen_az: (f64, f64),
    z0: f64,
}
impl Default for Gmt {
    fn default() -> Self {
        Self {
            source_zen_az: (0f64, 0f64),
            z0: 5_f64,
        }
    }
}
impl Gmt {
    pub fn new() -> Self {
        Default::default()
    }
    pub fn source_zen_arcsec_az_degree(self, source_zen_az: (f64, f64)) -> Self {
        Self {
            source_zen_az: (
                SkyAngle::Arcsecond(source_zen_az.0).to_radians(),
                SkyAngle::Degree(source_zen_az.1).to_radians(),
            ),
            ..self
        }
    }
    fn sample(&self) -> Delaunay {
        let radius = 8.365 * 0.5;
        let perimeter = 2. * std::f64::consts::PI * radius;
        let delta = 10e-2;
        let n = (perimeter / delta).ceil() as usize;
        let nodes: Vec<_> = (0..n)
            .flat_map(|i| {
                let o = 2. * std::f64::consts::PI * i as f64 / n as f64;
                vec![radius * o.cos(), radius * o.sin()]
            })
            .collect();
        let mut builder = Builder::new();
        builder.add_polygon(&nodes).set_switches("QDqa0.0025");
        let tri = builder.build();
        println!(
            "Surface area: {:.3}/{:.3}",
            std::f64::consts::PI * radius * radius,
            tri.area()
        );
        tri
    }
    fn z_xpupil() -> f64 {
        let m1 = Conic::gmt_m1();
        let m2 = Conic::gmt_m2();
        let mut ray = new_ray()
            .polar_direction_vector(SkyAngle::Arcminute(10f64).to_radians(), 0f64)
            .build();
        m1.reflect(&mut ray);
        ray.trace_to(&m2);
        m2.reflect(&mut ray);
        let z_xpupil = ray.solve_for_z(0f64, 0f64);
        println!("Exit pupil: {:.9}m", z_xpupil);
        z_xpupil
    }
    pub fn build(self) -> GmtInner {
        let tri = self.sample();
        let z_xpupil = Self::z_xpupil();
        let gmt = GmtInner {
            m1: Conic::gmt_m1(),
            m2: Conic::gmt_m2(),
            sphere: Conic::new(0f64, 0f64),
            tri,
            rays: Vec::new(),
            chief_opl: None,
        };
        let (zen, azi) = self.source_zen_az;
        let focal_rays = gmt.focal_point(vec![[8., 8., 0.]], zen, azi);
        println!("Chief focal rays: {}", focal_rays[0]);
        println!("Marg. focal rays: {}", focal_rays[1]);
        let ref_sphere_radius = [0., 0., z_xpupil].sub(focal_rays[0].p).norm();
        let ray = new_ray()
            .polar_direction_vector(self.source_zen_az.0, self.source_zen_az.1)
            .point_of_origin([0., 0., self.z0])
            .build();
        let mut gmt = GmtInner {
            m1: Conic::gmt_m1(),
            m2: Conic::gmt_m2(),
            sphere: Conic::new(0., -ref_sphere_radius).origin(focal_rays[0].p),
            tri: gmt.tri,
            rays: vec![ray],
            chief_opl: None,
        };
        let rays: Vec<_> = gmt
            .tri
            .vertex_iter()
            .map(|p| {
                new_ray()
                    .polar_direction_vector(zen, azi)
                    .point_of_origin([p[0], p[1], self.z0])
                    .build()
            })
            .collect();
        println!("# of rays: {}", rays.len());
        let chief_opl = gmt.trace().opd()[0];
        println!("Chief ray opl: {:.6e}", chief_opl);
        GmtInner {
            m1: Conic::gmt_m1(),
            m2: Conic::gmt_m2(),
            sphere: Conic::new(0., -ref_sphere_radius).origin(focal_rays[0].p),
            tri: gmt.tri,
            rays,
            chief_opl: Some(chief_opl),
        }
    }
}
pub struct GmtInner {
    pub m1: Conic,
    pub m2: Conic,
    pub sphere: Conic,
    pub tri: Delaunay,
    pub rays: Vec<Ray>,
    pub chief_opl: Option<f64>,
}
impl GmtInner {
    pub fn trace(&mut self) -> OPD {
        let m1 = &self.m1;
        let m2 = &self.m2;
        let sphere = &self.sphere;
        let chief_opl = &self.chief_opl;
        let opd = self
            .rays
            .iter_mut()
            .map(|mut ray| {
                let p0 = ray.p;
                ray.trace_to(m1);
                m1.reflect(&mut ray);
                let p1 = ray.p;
                ray.trace_to(m2);
                m2.reflect(&mut ray);
                let p2 = ray.p;
                ray.trace_to_with(sphere, |x1, x2| x1.min(x2));
                let p3 = ray.p;
                let distance = |u: &[f64], v: &[f64]| -> f64 {
                    u.iter()
                        .zip(v.iter())
                        .map(|(u, v)| {
                            let d = u - v;
                            d * d
                        })
                        .sum::<f64>()
                        .sqrt()
                };
                let a = distance(&p0, &p1);
                let b = distance(&p1, &p2);
                let c = distance(&p2, &p3);
                match chief_opl {
                    Some(opl) => a + b + c - opl,
                    None => a + b + c,
                }
            })
            .collect();
        OPD::new(&self.tri, self.rays.iter().map(|ray| ray.p).collect(), opd)
    }
    pub fn focal_point(&self, marginals: Vec<Vector>, z: f64, a: f64) -> Vec<Ray> {
        // Chief ray at M1 vertex from field angle (z,a)
        let chief_ray = new_ray().polar_direction_vector(z, a).build();
        // Marginal ray on M1 surface from field angle (z,a)
        let marginal_rays: Vec<Ray> = marginals
            .into_iter()
            .map(|m| {
                new_ray()
                    .point_of_origin(self.m1.height_at(m))
                    .polar_direction_vector(z, a)
                    .build()
            })
            .collect();
        // Collecting the rays
        let mut rays = vec![chief_ray];
        marginal_rays.into_iter().for_each(|m| rays.push(m));
        // Ray tracing to and reflecting from M2
        rays.iter_mut().for_each(|mut ray| {
            ray.trace_to(&self.m1);
            self.m1.reflect(&mut ray);
            ray.trace_to(&self.m2);
            self.m2.reflect(&mut ray);
        });
        // De-structuring chief ray
        let chief_ray = rays.remove(0);
        let p0 = chief_ray.p;
        let u0 = chief_ray.u;
        // De-structuring marginal rays
        let p: Vec<Vector> = rays.iter().map(|m| m.p).collect();
        let u: Vec<Vector> = rays.iter().map(|m| m.u).collect();
        // Build A matrix
        let n_u = u.len();
        let z: Vector = [0f64; 3];
        let mut cols: Vec<Vec<f64>> = vec![u0.to_vec(); n_u];
        for i_row in 0..n_u {
            let mut el = vec![];
            for i_col in 0..n_u {
                if i_row == i_col {
                    el.push(z.sub(u[i_col]).to_vec());
                } else {
                    el.push(z.clone().to_vec());
                }
            }
            cols.push(el.into_iter().flatten().collect());
        }
        let el: Vec<f64> = cols.into_iter().flatten().collect();
        let a = na::DMatrix::from_column_slice(n_u * 3, n_u + 1, &el);
        // Building b vector
        let b = na::DVector::from_vec(
            p.iter()
                .map(|x| x.sub(p0).to_vec())
                .flatten()
                .collect::<Vec<f64>>(),
        );
        // Solving As=b
        let s = a.svd(true, true).solve(&b, std::f64::EPSILON).unwrap();
        // Ray tracing to focal plane
        rays.insert(0, chief_ray);
        rays.iter_mut()
            .zip(s.into_iter())
            .for_each(|x| x.0.trace(*x.1));
        rays
    }
}

pub struct OPD<'a> {
    tri: &'a Delaunay,
    vertices: Vec<[f64; 3]>,
    opd: Vec<f64>,
    area: Option<f64>,
    mean: Option<f64>,
    std: Option<f64>,
}
impl<'a> OPD<'a> {
    pub fn new(tri: &'a Delaunay, vertices: Vec<[f64; 3]>, opd: Vec<f64>) -> Self {
        Self {
            tri,
            vertices,
            opd,
            area: None,
            mean: None,
            std: None,
        }
    }
    pub fn opd(&self) -> &[f64] {
        &self.opd
    }
    pub fn area(&mut self) -> f64 {
        match self.area {
            Some(area) => area,
            None => {
                let area = self.tri.mesh_area(&self.vertices);
                self.area = Some(area);
                area
            }
        }
    }
    pub fn mean(&mut self) -> f64 {
        match self.mean {
            Some(mean) => mean,
            None => {
                let mean = self.tri.average(&self.vertices, &self.opd) / self.area();
                self.mean = Some(mean);
                mean
            }
        }
    }
    pub fn std(&mut self) -> f64 {
        match self.std {
            Some(std) => std,
            None => {
                let mean = self.mean();
                let std = (self
                    .tri
                    .average_with(&self.vertices, &self.opd, |x| (x - mean).powf(2f64))
                    / self.area())
                .sqrt();
                self.std = Some(std);
                std
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn gmt() {
        let mut gmt = Gmt::new().build();
        let mut opd = gmt.trace();
        println!("STD: {:.6e}", opd.std());
    }
}
