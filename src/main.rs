use complot as plt;
use gmt_analytics::*;
use skyangle::Conversion;
use skyangle::SkyAngle::*;
use std::time::Instant;
use triangle_rs::Builder;
use zernike;

fn main() {
    let m1 = Conic::gmt_m1();
    let m2 = Conic::gmt_m2();

    println!("CHIEF RAY:");
    let mut ray = new_ray()
        .polar_direction_vector(Arcminute(10f64).to_radians(), 0f64)
        .build();
    println!("Init   : {}", ray);
    m1.reflect(&mut ray);
    println!("Reflect: {}", ray);
    ray.trace_to(&m2);
    println!("Trace  : {}", ray);
    m2.reflect(&mut ray);
    println!("Reflect: {}", ray);
    let z_xpupil = ray.solve_for_z(0f64, 0f64);
    println!("Exit pupil: {:.9}m", z_xpupil);

    let oa_arcmin = 10f64;

    println!("MARGINAL RAY:");
    let mut ray = new_ray()
        .point_of_origin(m1.height_at([10f64, 0f64, 0f64]))
        .polar_direction_vector(Arcminute(oa_arcmin).to_radians(), 0f64)
        .build();
    println!("Init   : {}", ray);
    m1.reflect(&mut ray);
    println!("Reflect: {}", ray);
    let z_gregorian = ray.solve_for_z(0f64, 0f64);
    println!("Gregorian focus: {:.9}m", z_gregorian);
    ray.trace_to(&m2);
    println!("Trace  : {}", ray);
    m2.reflect(&mut ray);
    println!("Reflect: {}", ray);
    /*let z_focal_plane = ray.solve_for_z(0f64, 0f64);
    println!("Focal plane: {:.9}m", z_focal_plane);*/

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
    let distances = |u: &[[f64; 3]], v: &[[f64; 3]]| -> Vec<f64> {
        u.iter()
            .zip(v.iter())
            .map(|(u, v)| distance(u, v))
            .collect()
    };

    let gmt = Gmt::new().build();
    let focal_rays = gmt.focal_point(vec![[8., 8., 0.]], Arcminute(oa_arcmin).to_radians(), 0f64);
    println!("Chief focal rays: {}", focal_rays[0]);
    println!("Marg. focal rays: {}", focal_rays[1]);

    let z0 = 5.;
    println!("CHIEF RAY:");
    let mut ray = new_ray()
        .polar_direction_vector(Arcminute(oa_arcmin).to_radians(), 0f64)
        .point_of_origin([0., 0., z0])
        .build();
    let p0 = ray.p;
    println!("Init   : {}", ray);
    ray.trace_to(&m1);
    println!("Trace  : {}", ray);
    m1.reflect(&mut ray);
    let p1 = ray.p;
    println!("Reflect: {}", ray);
    ray.trace_to(&m2);
    println!("Trace  : {}", ray);
    m2.reflect(&mut ray);
    let p2 = ray.p;
    println!("Reflect: {}", ray);
    //    let z_xpupil = ray.solve_for_z(0f64, 0f64);
    //    println!("Exit pupil: {:.9}m", z_xpupil);
    let ref_sphere_radius = [0., 0., z_xpupil].sub(focal_rays[0].p).norm();
    println!("Ref. sphere radius: {}", ref_sphere_radius);
    let ref_sphere = Conic::new(0., -ref_sphere_radius).origin(focal_rays[0].p);
    ray.trace_to_with(&ref_sphere, |x1, x2| x1.min(x2));
    let p3 = ray.p;
    println!("Sphere : {}", ray);

    let a = distance(&p0, &p1);
    let b = distance(&p1, &p2);
    let c = distance(&p2, &p3);
    let chief_opl = a + b + c;
    println!("Chief ray opl: {:.6e}", chief_opl);

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
    let tri = {
        let mut builder = Builder::new();
        builder.add_polygon(&nodes).set_switches("QDqa0.0025");
        builder.build()
    };
    println!(
        "Surface area: {:.3}/{:.3}",
        std::f64::consts::PI * radius * radius,
        tri.area()
    );

    let fig = plt::canvas("pupil.svg");
    let mut ax = plt::chart([-radius, radius, -radius, radius], &fig);
    plt::trimesh(&tri.x(), &tri.y(), [0; 3], &mut ax);

    let mut rays: Vec<_> = tri
        .vertex_iter()
        .map(|p| {
            new_ray()
                .polar_direction_vector(Arcminute(oa_arcmin).to_radians(), 0f64)
                .point_of_origin([p[0], p[1], z0])
                .build()
        })
        .collect();
    println!("# of rays: {}", rays.len());
    let now = Instant::now();
    let p0 = rays.coordinates();
    let p1 = rays.trace_to(&m1).reflect(&m1).coordinates();
    let p2 = rays.trace_to(&m2).reflect(&m2).coordinates();
    let p3 = rays
        .trace_to_with(&ref_sphere, |x1, x2| x1.min(x2))
        .coordinates();
    let a = distances(&p0, &p1);
    let b = distances(&p1, &p2);
    let c = distances(&p2, &p3);
    let et = now.elapsed();
    println!(
        "Tracing {} rays in {}ms ({:.3E}rays/s)",
        rays.len(),
        et.as_millis(),
        (rays.len() as f64 / et.as_secs_f64()).round() as u64
    );
    let opd: Vec<_> = a
        .iter()
        .zip(b.iter().zip(c.iter()))
        .map(|(a, (b, c))| a + b + c - chief_opl)
        .collect();
    let tip = oa_arcmin.from_arcmin();
    let tilt = 0f64;

    let vertices: Vec<_> = tri.vertex_iter().collect();

    let opd_tted: Vec<_> = opd
        .iter()
        .zip(tri.vertex_iter())
        .map(|(o, v)| o + tip * v[0] + tilt * v[1])
        .collect();

    let opd_std = (tri.triangle_iter().fold(0., |s, t| {
        let (a, b, c) = (&vertices[t[0]], &vertices[t[1]], &vertices[t[2]]);
        let ta = 0.5 * ((a[0] - c[0]) * (b[1] - a[1]) - (a[0] - b[0]) * (c[1] - a[1])).abs();
        let sa = t.iter().fold(0., |m, i| m + opd_tted[*i] * opd_tted[*i]) / 3 as f64;
        s + ta * sa
    }) / tri.area())
    .sqrt();
    println!("STD: {:.6e}", opd_std);

    let (j, n, m) = zernike::jnm(5);
    println!("zernike j: {:?}", j);
    let modes = {
        let modes: Vec<_> = (j.iter().zip(n.iter()))
            .zip(m.iter())
            .flat_map(|((j, n), m)| {
                tri.vertex_iter().map(move |v| {
                    let r = v[0].hypot(v[1]);
                    let o = v[1].atan2(v[0]);
                    zernike::zernike(*j, *n, *m, r, o)
                })
            })
            .collect();
        zernike::gram_schmidt_with_dot(&modes, j.len(), |u, v| tri.dot(u, v))
    };

    let shape = modes.chunks(tri.n_vertices()).nth(0).unwrap();
    let coefs: Vec<_> = modes
        .chunks(tri.n_vertices())
        .map(|mode| tri.dot(&opd_tted, &mode))
        .collect();
    let zern_var: f64 = coefs.iter().map(|x| x * x).sum();
    println!("Zernike WFE: {:.6e}", zern_var.sqrt());
    println!(
        "Zernike Tip-Tilt: [{:.3},{:.3}]",
        coefs[1].to_arcmin(),
        coefs[2].to_arcmin()
    );
    println!("Nomalized Zernike coefs:");
    coefs
        .iter()
        .zip(j.iter())
        .for_each(|(x, j)| println!("{:2}: {:+0.6}", j, (x / opd_std).powi(2)));

    let fig = plt::png_canvas("opd.png");
    let mut ax = plt::chart([-radius, radius, -radius, radius], &fig);
    plt::trimap(&tri.x(), &tri.y(), &opd_tted, &mut ax);
}
