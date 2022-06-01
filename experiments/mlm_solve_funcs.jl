using OrdinaryDiffEq

# steptol = 1e-6;
# termtol = 1e-9;

steptol = 1e-3;
termtol = 1e-6;

"""
    run_mlm(params; dt=x, tspan=(0.0,x))

    run MLM simulation with given parameters
    and save output to file
"""
function run_mlm(params; dt=3600.0*5.0, tspan=(0.0,3600.0*24.0*10.0))
    if params.fttype == fixedFT()
        params.qft0 = calc_qft0(params.RHft, params.Gamma_q, params.sft0, params.Gamma_s)
    end
    qtM0 = 0.7 * q_sat(0.0, params.SST0);
    hM0 = MixedLayerModel.Cp * params.SST0 + MixedLayerModel.L0 * qtM0;
    zi0 = 1000.0;
    CF0 = 1.0;
    u0 = [zi0, hM0, qtM0, params.SST0, CF0]; 
    prob = ODEProblem(ODEFunction(mlm, tgrad=(du, u, p, t) -> fill!(du, 0.0)), 
            u0, 
            tspan, 
            params);
    
    @time begin
        println("Rodas5");
        sol = solve(prob, Rodas5(autodiff=false), abstol=0.0, reltol=steptol, dtmax=dt);
    end

    # @time begin
    #     println("Euler");
    #     sol = solve(prob, Euler(), dt=dt);
    # end

    return u0, sol
end


"""
    run_mlm_from_init(u0, params; dt=x, tspan=(0.0,x))

    run MLM simulation with given parameters
    from initial state u0
    and save output to file
"""
function run_mlm_from_init(u0, params; dt=3600.0*5.0, tspan=(0.0,3600.0*24.0*10.0))
    if params.fttype == fixedFT()
        params.qft0 = calc_qft0(params.RHft, params.Gamma_q, params.sft0, params.Gamma_s)
    end
    prob = ODEProblem(ODEFunction(mlm, tgrad=(du, u, p, t) -> fill!(du, 0.0)), 
            u0, 
            tspan, 
            params);

    @time begin
        println("Rodas5");
        sol = solve(prob, Rodas5(autodiff=false), abstol=0.0, reltol=steptol, dtmax=dt);
    end

    # @time begin
    #     println("Euler");
    #     sol = solve(prob, Euler(), dt=dt);
    # end

    return u0, sol
end
