export run, run_with_output

"""
    run(params, print=false)

    run MLM simulation with given parameters
"""
function run(params)
    params.qft0 = calc_qft0(params.RHft, params.Gamma_q, params.sft0, params.Gamma_s);
    
    z0 = 0.0;
    qtM0 = params.RHsurf * q_sat(z0, params.SST0);
    hM0 = Cp * (params.SST0 - 2.0) + L0 * qtM0;
    
    zi0 = 900.0
    u0 = [zi0, hM0, qtM0, params.SST0];

    prob = SteadyStateProblem(mlm, u0, params);
    tspan = 3600.0 * 24.0 * 100.0;
    tol = 1e-6;
    
    sol = solve(prob, DynamicSS(Rosenbrock23(autodiff=false);abstol=tol,reltol=0.0,tspan=tspan));
    uf = sol.u;
    
    return u0, uf, params, sol
end

"""
    run_with_output(params)

    run MLM simulation and print output to file
"""
function run_with_output(params, filename="default.txt")
    u0, uf, par, sol = run(params);
    du = zeros(4);
    mlm(du, uf, par, 0.0);
    zi,hM,qM,SST = uf;
    
    filename = string("experiments/output/",filename)
    open(filename, "a") do io
        write(io, string(sol.retcode,"\n"))
        write(io, string("u0: ",u0,"\n"))
        write(io, string("uf: ",uf,"\n"))
        write(io, string("du/u: ", du ./ uf,"\n"))
        write(io, string("we = ", we(uf,par,par.etype)*1000, " (mm/s)\n"))
        write(io, string("zc = ", zi - calc_LCL(zi,hM,qM)," (m)\n"))
        write(io, string("RH0 = ", RH(0.0,hM,qM)*100.0," (%)\n"))
        write(io, string("LHF = ", calc_LHF(uf,par)," (W/m2)\n"))
        write(io, string("SHF = ", calc_SHF(uf,par)," (W/m2)\n"))
        write(io, string("dR = ", calc_cloudtop_RAD(uf,par,par.rtype)," (W/m^2)\n"))
    end;

    return u0, uf, par
end