use std::vec::Vec;
use ndarray::{Array};
#[allow(non_snake_case)]
pub fn F(left: f64, right: f64, c: f64) -> f64 {
	return 0.5*(c+c.abs())*left + 0.5*(c-c.abs())*right;
    }

type Field<D> = Array<f64,D>;

pub struct Solver<D> {
    pub field: Field<D>,
    pub courant_number: Field<D>,
    pub c_err: Field<D>,
    steps_done: usize,
}

use ndarray::{Ix1, Ix2, Ix3};
pub type Field1D = Field<Ix1>;
pub type Field2D = Field<Ix2>;
pub type Field3D = Field<Ix3>;

pub type Solver1D = Solver<Ix1>;
pub type Solver2D = Solver<Ix2>;
pub type Solver3D = Solver<Ix3>;


impl Solver1D {
    fn get_c_err (&mut self) {
	let mut c = self.courant_number.clone();
	let len = self.field.len();
	for (x, val) in c.indexed_iter_mut() {
	    if x > 0 && x < len-1 {
		*val = ((*val).abs()-(*val)*(*val)) * (self.field[x]-self.field[x-1])/(self.field[x]+self.field[x-1]);
	    }
	}
	self.c_err = c;
    }
    
    fn upwind(&mut self, mpdata: bool)  {
	let s = &self.field.clone();
	let len = self.field.len();
	let c: &Field1D;
	if !mpdata {
	    c = &self.courant_number;
	} else {
	    c = &self.c_err;
	}
	for (x, psi) in self.field.indexed_iter_mut() {
	    if x > 0 && x < len-1 {
		*psi = *psi - (
		    F(s[x], s[x+1], c[x]) -
		    F(s[x-1], s[x], c[x-1])
		);
	    }
	}
	self.steps_done += 1;
    }

    pub fn initial_condition(&mut self, f: fn(f64, &Vec<f64>) -> f64, params: Vec<f64>) {
	self.field = self.field.map(|x| f(*x, &params));
    }
    
    pub fn advance(&mut self, time_steps: usize, method: &str) -> ()
    {
	match method {
	    "MPDATA" => {
		for _ in 0..time_steps-1 {
		    self.upwind(false);
		    self.get_c_err();
		    self.upwind(true);
		}
	    }
	    "upwind" => {
		for _ in 0..time_steps-1 {
		    self.upwind(false);
		}
	    }
	    _ => { panic!("No method chosen! Availible methods: \"MPDATA\" and \"upwind\".") }
	}
    }    
    pub fn assemble(vel: f64, dt: f64, dx: f64, len: usize, linstart: f64, linend: f64, nsteps: usize) -> Self {
	Solver1D {
	    field: Field1D::linspace(linstart, linend, nsteps),
	    courant_number: (vel*dt)/dx*Field1D::ones(len),
	    c_err: (vel*dt)/dx*Field1D::ones(len),
	    steps_done: 0
	}
    } 
}

