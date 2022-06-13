# Solving the coupled MLM ODEs

We want to solve the MLM in two scenarios: 1) from a specified initial steady-state condition or 2) from an arbitrary initial guess.

## From initial-condition
In this first case, we are taking the steady-state solution from a prior simulation of some climate state (usually present-day CO``_2`` levels of 400 ppm) and perturbing the CO``_2`` and letting the system evolve to reach a new steady state. We will need to do `using DifferentialEquations` to define the `ODEProblem`.

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
Alternatively, we could specify this as a `SteadyStateProblem` like this:
```julia
function run_mlm_ss_from_init(u0, params; dt=3600.0*5.0, tspan=3600.0*24.0*10.0)
    prob = SteadyStateProblem(mlm, u0, params);
    tol = 1e-9;
    @time begin
        sol = solve(prob, DynamicSS(Euler(); abstol=0.0, reltol=tol, tspan=tspan), dt=dt, progress=true, progress_steps=50);
    end
    return u0, sol
end
```


## No initial condition
On the other hand, we sometimes may want to solve the MLM without having an initial condition in mind and we will need to make a good guess to start. We could do that as follows.
```julia
function run_mlm(params; dt=3600.0*5.0, tspan=(0.0,3600.0*24.0*10.0))
    zi0 = 1200.0
    qtM0 = params.RHsurf * q_sat(0.0, params.SST0);
    sM0 = MixedLayerModel.Cp * params.SST0 + MixedLayerModel.L0 * qtM0;
    u0 = [zi0, sM0, qtM0, params.SST0];
    prob = ODEProblem(mlm, u0, tspan, params);
    @time begin
        sol = solve(prob, Euler(), dt=dt, progress=true, progress_steps=50);
    end
    return u0, sol
end
```
And like above, we could also solve this directly for the steady-state using `SteadyStateProblem` rather than `ODEProblem`.
```julia
function run_mlm_ss(params; dt=3600.0*5.0, tspan=3600.0*24.0*10.0)    
    zi0 = 1200.0;
    qtM0 = params.RHsurf * q_sat(0.0, params.SST0);
    sM0 = MixedLayerModel.Cp * params.SST0 + MixedLayerModel.L0 * qtM0;
    u0 = [zi0, sM0, qtM0, params.SST0];

    prob = SteadyStateProblem(mlm, u0, params);
    tol = 1e-9;
    @time begin
        sol = solve(prob, DynamicSS(Euler(); abstol=0.0, reltol=tol, tspan=tspan), dt=dt, progress=true, progress_steps=50);
    end
    return u0, sol
end
```