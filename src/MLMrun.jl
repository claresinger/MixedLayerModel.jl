module MLMrun

include("Definitions.jl")
using DifferentialEquations
using ..Thermodynamics
using ..SurfaceFluxes
using ..Radiation
using ..Entrainment
using ..MLMODE

function run(params, print=false)    
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
    
    if print
        println(sol.retcode);
    end

    return u0, uf, params
end

function run_with_output(params)
    u0, uf, par = run(params, true);
    
    println("CO2: ",par.CO2);
    println("u0: ",u0);
    println("uf: ",uf);
    
    du = zeros(4);
    mlm(du, uf, par, 0.0);
    println("du/u: ", du ./ uf);
    zi,hM,qM,SST = uf;
    
    println("we = ", we(uf,par,par.etype)*1000, " (mm/s)");
    println("zc = ", zi - calc_LCL(zi,hM,qM)," (m)");
    println("RH0 = ", RH(0.0,hM,qM)*100.0," (%)");
    println("LHF = ", calc_LHF(uf,par)," (W/m2)");
    println("SHF = ", calc_SHF(uf,par)," (W/m2)");
    println("LWP = ", calc_LWP(zi,hM,qM)*1000.0," (g/m^2)");
    println("ΔR = ", calc_cloud_RAD(uf,par)," (W/m^2)");
    println();
    
    return u0, uf, par
end

end