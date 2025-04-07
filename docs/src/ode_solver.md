# Solving the coupled MLM ODEs

We want to solve the MLM in two scenarios: 1) from a specified initial steady-state condition or 2) from an arbitrary initial guess.

## From initial-condition
In this first case, we are taking the steady-state solution from a prior simulation of some climate state (usually present-day CO``_2`` levels of 400 ppm) and perturbing the CO``_2`` and letting the system evolve to reach a new steady state. We will need to do `using OrdinaryDiffEq` to define the `ODEProblem`.

This is one way we could code this with an explicit timestep of 5 hours, running for 10 days. We are using the 
```julia
function run_mlm_from_init(u0, params; dt=3600.0*5.0, tspan=(0.0,3600.0*24.0*10.0))
    prob = ODEProblem(mlm, u0, tspan, params);
    @time begin
        sol = solve(prob, Euler(), dt=dt, progress=true, progress_steps=50);
    end
    return u0, sol
end
```

## No initial condition
On the other hand, we sometimes may want to solve the MLM without having an initial condition in mind and we will need to make a good guess to start. We could do that as follows.
```julia
function run_mlm(params; dt=3600.0*5.0, tspan=(0.0,3600.0*24.0*10.0))
    zi0 = 1000.0
    qtM0 = params.RHsurf * q_sat(0.0, params.SST0);
    sM0 = MixedLayerModel.Cp * params.SST0;
    CF0 = 1.0;
    u0 = [zi0, sM0, qtM0, params.SST0, CF0];
    prob = ODEProblem(mlm, u0, tspan, params);
    @time begin
        sol = solve(prob, Euler(), dt=dt, progress=true, progress_steps=50);
    end
    return u0, sol
end
```
