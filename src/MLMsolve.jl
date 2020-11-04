export run_mlm, run_mlm_from_init

"""
    run_mlm(params, filename="default.txt")

    run MLM simulation with given parameters
    and save output to file
"""
function run_mlm(params)
    params.qft0 = calc_qft0(params.RHft, params.Gamma_q, params.sft0, params.Gamma_s);
    
    z0 = 0.0;
    qtM0 = params.RHsurf * q_sat(z0, params.SST0);
    hM0 = Cp * (params.SST0 - 2.0) + L0 * qtM0;
    
    zi0 = 1200.0
    u0 = [zi0, hM0, qtM0, params.SST0];

    prob = SteadyStateProblem(mlm, u0, params);
    tspan = 3600.0 * 24.0 * 1000.0;
    tol = 1e-6;

    sol = solve(prob, DynamicSS(Rosenbrock23(autodiff=false);abstol=tol,reltol=0.0,tspan=tspan));
    #sol = solve(prob, DynamicSS(Rosenbrock23(autodiff=false);abstol=1e-10,reltol=1e-10,tspan=tspan));
    #sol = solve(prob, DynamicSS(CVODE_BDF();abstol=1e-10,reltol=1e-10), dt=1.0);

    return u0, sol
end

"""
    run_mlm_from_init(u0, params, filename="default.txt")

    run MLM simulation with given parameters
    from initial state u0
    and save output to file
"""
function run_mlm_from_init(u0, params, filename="default.txt")
    params.qft0 = calc_qft0(params.RHft, params.Gamma_q, params.sft0, params.Gamma_s);

    prob = SteadyStateProblem(mlm, u0, params);
    tspan = 3600.0 * 24.0 * 100.0;
    tol = 1e-6;
    
    sol = solve(prob, DynamicSS(Rosenbrock23(autodiff=false);abstol=tol,reltol=0.0,tspan=tspan));

    return u0, sol
end