"""
    run_mlm(params; dt=x, tspan=(0.0,x))

    run MLM simulation with given parameters
    and save output to file
"""
function run_mlm(params; dt=3600.0*5.0, tspan=(0.0,3600.0*24.0*10.0))
    #params.qft0 = calc_qft0(params.RHft, params.Gamma_q, params.sft0, params.Gamma_s);
    
    z0 = 0.0;
    qtM0 = params.RHsurf * q_sat(z0, params.SST0);
    hM0 = MixedLayerModel.Cp * params.SST0 + MixedLayerModel.L0 * qtM0;
    
    zi0 = 1200.0
    u0 = [zi0, hM0, qtM0, params.SST0];
        
    prob = ODEProblem(mlm, u0, tspan, params);

    @time begin
        sol = solve(prob, Euler(), dt=dt, progress=true, progress_steps=50);
    end

    return u0, sol
end

"""
    run_mlm_from_init(u0, params; dt=x, tspan=(0.0,x))

    run MLM simulation with given parameters
    from initial state u0
    and save output to file
"""
function run_mlm_from_init(u0, params; dt=3600.0*5.0, tspan=(0.0,3600.0*24.0*10.0))
    #params.qft0 = calc_qft0(params.RHft, params.Gamma_q, params.sft0, params.Gamma_s);

    prob = ODEProblem(mlm, u0, tspan, params);

    @time begin
        sol = solve(prob, Euler(), dt=dt, progress=true, progress_steps=50);
    end

    return u0, sol
end

"""
    run_mlm_ss(params; dt=x, tspan=x)

    run MLM simulation with given parameters
    and save output to file
"""
function run_mlm_ss(params; dt=3600.0*5.0, tspan=3600.0*24.0*10.0)
    #params.qft0 = calc_qft0(params.RHft, params.Gamma_q, params.sft0, params.Gamma_s);
    
    z0 = 0.0;
    qtM0 = params.RHsurf * q_sat(z0, params.SST0);
    hM0 = MixedLayerModel.Cp * params.SST0 + MixedLayerModel.L0 * qtM0;
        
    zi0 = 1200.0
    u0 = [zi0, hM0, qtM0, params.SST0];

    prob = SteadyStateProblem(mlm, u0, params);
    tol = 1e-9;

    @time begin
        # println("rootfind");
        # sol = solve(prob, SSRootfind());

        println("euler, dt=4 hrs");
        sol = solve(prob, DynamicSS(Euler(); abstol=0.0, reltol=tol, tspan=tspan), dt=dt, progress=true, progress_steps=50);

        # println("tsit5");
        # sol = solve(prob,DynamicSS(Tsit5()))
        # sol = solve(prob, DynamicSS(Tsit5(); abstol=tol, reltol=0.0, tspan=tspan));

        # println("rosenbrock23");
        # sol = solve(prob, DynamicSS(Rosenbrock23(autodiff=false); abstol=tol, reltol=0.0, tspan=tspan));
        
        # println("cvode_bdf");
        # dt = 0.1;
        # println("dt=",dt);
        # sol = solve(prob, DynamicSS(CVODE_BDF(); abstol=tol, reltol=0.0), dt=dt);
    end

    return u0, sol
end

"""
    run_mlm_ss_from_init(u0, params; dt=x, tspan=x)

    run MLM simulation with given parameters
    from initial state u0
    and save output to file
"""
function run_mlm_ss_from_init(u0, params; dt=3600.0*5.0, tspan=3600.0*24.0*10.0)
    #params.qft0 = calc_qft0(params.RHft, params.Gamma_q, params.sft0, params.Gamma_s);

    prob = SteadyStateProblem(mlm, u0, params);
    tol = 1e-9;

    @time begin
        # println("rootfind");
        # sol = solve(prob, SSRootfind());

        println("euler");
        sol = solve(prob, DynamicSS(Euler(); abstol=0.0, reltol=tol, tspan=tspan), dt=dt, progress=true, progress_steps=50);

        # println("tsit5");
        # sol = solve(prob, DynamicSS(Tsit5(); abstol=0.0, reltol=tol, tspan=tspan));

        # println("rosenbrock23");
        # sol = solve(prob, DynamicSS(Rosenbrock23(autodiff=false); abstol=0.0, reltol=tol, tspan=tspan));
        
        # println("cvode_bdf");
        # sol = solve(prob, DynamicSS(CVODE_BDF(); abstol=0.0, reltol=tol), dt=3600.0, progress=true, progress_steps=50);
    end

    return u0, sol
end