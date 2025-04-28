use rust_mpdata::field::*;

fn main() {
    //First, let's define the conditions of the simulation: velocity, temporal step, spacial step and number of simulated steps: 
    let vel = 2 as f64;
    let sim_time_max = 66 as f64;
    let n_steps = 100 as usize;
    let dt = sim_time_max/(n_steps as f64);
    let first_x = -100 as f64;
    let last_x = 300 as f64;
    let x_step = 45 as usize;
    let dx = (last_x-first_x)/(x_step as f64);

    //Now, we can assemble the solver:
    let mut solver = Solver::assemble(vel, dt, dx, first_x, last_x, x_step);
    //Then let's apply the initial condition:
    let params: Vec<f64> = vec! [ 1., 0., 20., ]; 
    solver.initial_condition(init_cond as fn(f64, &Vec<f64>)-> f64, params);
    //And take a quick look: (TODO: plot?)
    println!("{:?}", solver.field);
    //Now, let's advect the field with MPDATA
    solver.advance(n_steps, "MPDATA");
    //And see the difference:
    println!("{:?}", solver.field);
}

fn init_cond(x: f64, params: &Vec<f64>) -> f64 {
    //This function will help us create an initial state for the simulation
    let a = params[0];
    let x0 = params[1];
    let sigma = params[2];
    a*(-(x-x0)*(x-x0)/2./(sigma*sigma)).exp()
}
