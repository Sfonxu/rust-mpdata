use std::fs::File;
use std::path::Path;
use ndarray_npy::write_npy;
const GRID_LEN: usize = 45;
//type Field1D = Array::<f64, Ix1>;

#[allow(unused_imports)]
use rust_mpdata::field::*;

fn initial_condition(x: f64, params: &Vec<f64>) -> f64 {
    let a = params[0];
    let x0 = params[1];
    let sigma = params[2];
    return a*(-((x-x0)*(x-x0))/2./(sigma*sigma)).exp();
}

fn main() {
    let vel = 2.;
    let dt = 66./100.;
    let dx = 9.1;
    let linstart = -100.;
    let linend = 300.;
    let nsteps = 45;
    let params = vec! [ 1., 0., 20. ];
    let mut s = Solver1D::assemble(vel, dt, dx, GRID_LEN, linstart, linend, nsteps);
    let mut s_mpdata = Solver1D::assemble(vel, dt, dx, GRID_LEN, linstart, linend, nsteps);
    s.initial_condition(initial_condition, params.clone());
    s_mpdata.initial_condition(initial_condition, params);
    let _ = File::create("initial.npy").expect("failed to create file");
    let path2 = Path::new("initial.npy");
    let _ = write_npy(path2, &s.field);
    s.advance(100, "upwind");
    s_mpdata.advance(100, "MPDATA");
    let _ = File::create("output.npy").expect("failed to load file");
    let _ = File::create("output_mpdata.npy").expect("failed to load file");
    let path = Path::new("output.npy");
    let path3 = Path::new("output_mpdata.npy");
    let _ = write_npy(path, &s.field);
    let _ = write_npy(path3, &s_mpdata.field);
}
