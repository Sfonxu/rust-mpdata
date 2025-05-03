use std::vec::Vec;
use ndarray::{Array, Zip, s};
#[allow(non_snake_case)]
pub fn F(left: f64, right: f64, c: f64) -> f64 {
	return 0.5*(c+c.abs())*left + 0.5*(c-c.abs())*right;
    }

type Field<D> = Array<f64,D>;

pub struct Solver<D> {
    pub field: Field<D>,
    courant_number: Field<D>,
    c_err: Field<D>,
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
	Zip::from(self.c_err.slice_mut(s![1..]))
	    .and(self.field.windows(2))
	    .for_each(|c, w| {
		*c = ((*c).abs()-(*c)*(*c))*(w[0]-w[1])/(w[0]+w[1]);
	    });	
    }

    
    fn upwind(&mut self, mpdata: bool)  {
	let c: &Field1D;
	if !mpdata {
	    c = &self.courant_number;
	} else {
	    self.get_c_err();
	    c = &self.c_err;
	}
	for x in 1..self.field.len()-1 {
	    self.field[x] = self.field[x]-(F(self.field[x],self.field[x+1],c[x+1])-F(self.field[x-1],self.field[x],c[x]));
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
    pub fn assemble(vel: f64, dt: f64, dx: f64, linstart: f64, linend: f64, nsteps: usize) -> Self {
	Solver1D {
	    field: Field1D::linspace(linstart, linend, nsteps),
	    courant_number: (vel*dt)/dx*Field1D::ones(nsteps),
	    c_err: (vel*dt)/dx*Field1D::ones(nsteps),
	    steps_done: 0
	}
    } 
}

