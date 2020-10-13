export run_mlm, run_mlm_from_init

"""
    run_mlm(params, filename="default.txt")

    run MLM simulation with given parameters
    and save output to file
"""
function run_mlm(params, filename="default.txt")
    params.qft0 = calc_qft0(params.RHft, params.Gamma_q, params.sft0, params.Gamma_s);
    
    z0 = 0.0;
    qtM0 = params.RHsurf * q_sat(z0, params.SST0);
    hM0 = Cp * (params.SST0 - 2.0) + L0 * qtM0;
    
    zi0 = 900.0
    u0 = [zi0, hM0, qtM0, params.SST0];

    prob = SteadyStateProblem(mlm, u0, params);    
    sol = solve(prob, DynamicSS(Rosenbrock23(autodiff=false);abstol=1e-6,reltol=0.0,tspan=3600.0 * 24.0 * 100.0));
    
    #sol = solve(prob, DynamicSS(Rosenbrock23(autodiff=false);abstol=1e-10,reltol=1e-10,tspan=tspan));
    #sol = solve(prob, DynamicSS(CVODE_BDF();abstol=1e-10,reltol=1e-10), dt=1.0);
    
    uf = sol.u;
    println(sol.retcode);
    
    du = zeros(4);
    mlm(du, uf, params, 0.0);
    zi,hM,qM,SST = uf;

    filename = string("experiments/output/",filename)
    open(filename, "a") do io
        write(io, string(sol.retcode,"\n"))
        write(io, string("u0: ",u0,"\n"))
        write(io, string("uf: ",uf,"\n"))
        write(io, string("du/u: ", du ./ uf,"\n"))
        write(io, string("we = ", we(uf,params,params.etype)*1000, " (mm/s)\n"))
        write(io, string("zc = ", zi - calc_LCL(zi,hM,qM)," (m)\n"))
        write(io, string("RH0 = ", RH(0.0,hM,qM)*100.0," (%)\n"))
        write(io, string("LHF = ", calc_LHF(uf,params)," (W/m2)\n"))
        write(io, string("SHF = ", calc_SHF(uf,params)," (W/m2)\n"))
        write(io, string("dR = ", calc_cloudtop_RAD(uf,params,params.rtype)," (W/m^2)\n"))
        write(io, string("OHU = ", calc_OHU(uf,params,params.stype)," (W/m^2)\n"))
    end;

    return u0, uf
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
    tspan = 3600.0 * 24.0 * 10.0;
    tol = 1e-6;
    
    sol = solve(prob, DynamicSS(Rosenbrock23(autodiff=false);abstol=tol,reltol=0.0,tspan=tspan));
    uf = sol.u;
    println(sol.retcode);
    
    du = zeros(4);
    mlm(du, uf, params, 0.0);
    zi,hM,qM,SST = uf;

    filename = string("experiments/output/",filename)
    open(filename, "a") do io
        write(io, string(sol.retcode,"\n"))
        write(io, string("u0: ",u0,"\n"))
        write(io, string("uf: ",uf,"\n"))
        write(io, string("du/u: ", du ./ uf,"\n"))
        write(io, string("we = ", we(uf,params,params.etype)*1000, " (mm/s)\n"))
        write(io, string("zc = ", zi - calc_LCL(zi,hM,qM)," (m)\n"))
        write(io, string("RH0 = ", RH(0.0,hM,qM)*100.0," (%)\n"))
        write(io, string("LHF = ", calc_LHF(uf,params)," (W/m2)\n"))
        write(io, string("SHF = ", calc_SHF(uf,params)," (W/m2)\n"))
        write(io, string("dR = ", calc_cloudtop_RAD(uf,params,params.rtype)," (W/m^2)\n"))
        write(io, string("OHU = ", calc_OHU(uf,params,params.stype)," (W/m^2)\n"))
    end;

    return u0, uf
end