using OrdinaryDiffEq
using StochasticDiffEq, DiffEqNoiseProcess

steptol = 1e-3;
termtol = 1e-6;

"""
    run_mlm(params; dt=x, tspan=(0.0,x))

    run MLM simulation with given parameters
    and save output to file
"""
function run_mlm(params; init=1, dt=3600.0*5.0, tspan=(0.0,3600.0*24.0*10.0), quiet=false)
    qtM0 = 0.7 * q_sat(0.0, params.SST0);
    sM0 = MixedLayerModel.Cp * params.SST0;
    if init == 1
        zi0 = 1000.0;
        CF0 = params.CFmax;
    else
        zi0 = 1000.0;
        CF0 = params.CFmin;
    end 
    u0 = [zi0, sM0, qtM0, params.SST0, CF0]; 

    prob = ODEProblem(ODEFunction(mlm, tgrad=(du, u, p, t) -> fill!(du, 0.0)), 
            u0, 
            tspan, 
            params);
    
    if quiet
        sol = solve(prob, Rodas5(autodiff=false), abstol=0.0, reltol=steptol, dtmax=dt);
    else
        @time begin
            println("Rodas5");
            sol = solve(prob, Rodas5(autodiff=false), abstol=0.0, reltol=steptol, dtmax=dt);
        end
    end

    # sol = solve(prob, Euler(), dt=dt);

    return u0, sol
end


"""
    run_mlm_from_init(u0, params; dt=x, tspan=(0.0,x))

    run MLM simulation with given parameters
    from initial state u0
    and save output to file
"""
function run_mlm_from_init(u0, params; dt=3600.0*5.0, tspan=(0.0,3600.0*24.0*10.0), quiet=false)    
    prob = ODEProblem(ODEFunction(mlm, tgrad=(du, u, p, t) -> fill!(du, 0.0)), 
            u0, 
            tspan, 
            params);

    if quiet
        sol = solve(prob, Rodas5(autodiff=false), abstol=0.0, reltol=steptol, dtmax=dt);
    else
        @time begin
            println("Rodas5");
            sol = solve(prob, Rodas5(autodiff=false), abstol=0.0, reltol=steptol, dtmax=dt);
        end
    end

    # sol = solve(prob, Euler(), dt=dt);

    return u0, sol
end


"""
    run_mlm(params; dt=x, tspan=(0.0,x))

    run MLM simulation with given parameters
    and save output to file
"""
function stochastic_run_mlm(params; init=1, dt=3600.0*5.0, tspan=(0.0,3600.0*24.0*10.0), quiet=false)
    qtM0 = 0.7 * q_sat(0.0, params.SST0);
    sM0 = MixedLayerModel.Cp * params.SST0;
    if init == 1
        zi0 = 1000.0;
        CF0 = params.CFmax;
    else
        zi0 = 1000.0;
        CF0 = params.CFmin;
    end 
    u0 = [zi0, sM0, qtM0, params.SST0, CF0]; 

    τsyn = 10*24*3600. # 5 days
    W = OrnsteinUhlenbeckProcess(1/τsyn, 0., 1., 0., 0.)
    prob = SDEProblem(
        SDEFunction(mlm, g_bounded, tgrad=(du, u, p, t) -> fill!(du, 0.0)),
        g_bounded,
        u0,
        tspan,
        params,
        noise = W,
    );
    ensembleprob = EnsembleProblem(prob);
    sim = solve(ensembleprob, EM(), EnsembleThreads(), trajectories=10, dt=dt);
    sol = EnsembleSummary(sim);

    # prob = SDEProblem(
    #     SDEFunction(mlm_stochastic, g_bounded, tgrad=(du, u, p, t) -> fill!(du, 0.0)),
    #     g_bounded,
    #     u0,
    #     tspan,
    #     params,
    #     seed = 18,
    # );
    # sol = solve(prob, EM(), dt=dt);
    # sol = solve(prob, SOSRI(), abstol=0.0, reltol=steptol, dtmax=dt);

    return u0, sim, sol
end


"""
    run_mlm_from_init(u0, params; dt=x, tspan=(0.0,x))

    run MLM simulation with given parameters
    from initial state u0
    and save output to file
"""
function stochastic_run_mlm_from_init(u0, params; dt=3600.0*5.0, tspan=(0.0,3600.0*24.0*10.0), quiet=false)    
    τsyn = 10*24*3600. # 5 days
    W = OrnsteinUhlenbeckProcess(1/τsyn, 0., 1., 0., 0.)
    prob = SDEProblem(
        SDEFunction(mlm, g_bounded, tgrad=(du, u, p, t) -> fill!(du, 0.0)),
        g_bounded,
        u0,
        tspan,
        params,
        noise = W,
    );
    ensembleprob = EnsembleProblem(prob);
    sim = solve(ensembleprob, EM(), EnsembleThreads(), trajectories=10, dt=dt);
    sol = EnsembleSummary(sim);
    
    # prob = SDEProblem(
    #     SDEFunction(mlm_stochastic, g_bounded, tgrad=(du, u, p, t) -> fill!(du, 0.0)),
    #     g_bounded,
    #     u0,
    #     tspan,
    #     params,
    #     seed = round(u0[5]*10000),
    # );
    # sol = solve(prob, EM(), dt=dt);
    # sol = solve(prob, SOSRI(), abstol=0.0, reltol=0.1, dtmax=dt);
    
    return u0, sim, sol
end
