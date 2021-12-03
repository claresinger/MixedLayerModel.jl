using OrdinaryDiffEq
using SteadyStateDiffEq

# steptol = 1e-6;
# termtol = 1e-9;

steptol = 1e-3;
termtol = 1e-6;

RHsurf = 0.65;

"""
    run_mlm(params; dt=x, tspan=(0.0,x))

    run MLM simulation with given parameters
    and save output to file
"""
function run_mlm(params; dt=3600.0*5.0, tspan=(0.0,3600.0*24.0*10.0))
    if params.fttype == fixedFT()
        params.qft0 = calc_qft0(params.RHft, params.Gamma_q, params.sft0, params.Gamma_s)
    end
    qtM0 = RHsurf * q_sat(0.0, params.SST0);
    hM0 = MixedLayerModel.Cp * params.SST0 + MixedLayerModel.L0 * qtM0;
    zi0 = 1100.0;
    CF0 = 1.0;
    u0 = [zi0, hM0, qtM0, params.SST0, CF0]; 
    prob = ODEProblem(ODEFunction(mlm, tgrad=(du, u, p, t) -> fill!(du, 0.0)), 
            u0, 
            tspan, 
            params);
    
    @time begin
        println("Rodas5");
        sol = solve(prob, Rodas5(autodiff=false), abstol=0.0, reltol=steptol);
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
        sol = solve(prob, Rodas5(autodiff=false), abstol=0.0, reltol=steptol);
    end

    # @time begin
    #     println("Euler");
    #     sol = solve(prob, Euler(), dt=dt);
    # end

    return u0, sol
end

# """
#     run_mlm_ss(params; dt=x, tspan=x)

#     run MLM simulation with given parameters
#     and save output to file
# """
# function run_mlm_ss(params; dt=3600.0*5.0, tspan=3600.0*24.0*10.0)
#     if params.fttype == fixedFT()
#         params.qft0 = calc_qft0(params.RHft, params.Gamma_q, params.sft0, params.Gamma_s)
#     end
#     qtM0 = RHsurf * q_sat(0.0, params.SST0);
#     hM0 = MixedLayerModel.Cp * params.SST0 + MixedLayerModel.L0 * qtM0; 
#     zi0 = 1200.0;
#     CF0 = 1.0;
#     u0 = [zi0, hM0, qtM0, params.SST0, CF0];
#     prob = SteadyStateProblem(mlm, u0, params);

#     @time begin
#         println("Rodas5");
#         sol = solve(prob, DynamicSS(Rodas5(); abstol=0.0, reltol=termtol, tspan=tspan), abstol=0.0, reltol=steptol);
#     end

#     return u0, sol
# end

# """
#     run_mlm_ss_from_init(u0, params; dt=x, tspan=x)

#     run MLM simulation with given parameters
#     from initial state u0
#     and save output to file
# """
# function run_mlm_ss_from_init(u0, params; dt=3600.0*5.0, tspan=3600.0*24.0*10.0)
#     if params.fttype == fixedFT()
#         params.qft0 = calc_qft0(params.RHft, params.Gamma_q, params.sft0, params.Gamma_s)
#     end
#     prob = SteadyStateProblem(mlm, u0, params);

#     @time begin
#         # println("euler, dt="*string(dt/3600.0)*" hrs, tmax="*string(tspan/3600.0/24.0)*" days");
#         # sol = solve(prob, DynamicSS(Euler(); abstol=0.0, reltol=tol, tspan=tspan), dt=dt, progress=true, progress_steps=50);

#         println("Rodas5");
#         sol = solve(prob, DynamicSS(Rodas5(); abstol=0.0, reltol=termtol, tspan=tspan), abstol=0.0, reltol=steptol);

#         # println("tsit5");
#         # sol = solve(prob,DynamicSS(Tsit5()))
#         # sol = solve(prob, DynamicSS(Tsit5(); abstol=tol, reltol=0.0, tspan=tspan));

#         # println("rosenbrock23");
#         # sol = solve(prob, DynamicSS(Rosenbrock23(autodiff=false); abstol=tol, reltol=0.0, tspan=tspan));
        
#         # println("cvode_bdf");
#         # dt = 0.1;
#         # println("dt=",dt);
#         # sol = solve(prob, DynamicSS(CVODE_BDF(); abstol=tol, reltol=0.0), dt=dt);
#     end

#     return u0, sol
# end