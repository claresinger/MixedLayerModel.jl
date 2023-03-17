using MixedLayerModel
using MixedLayerModel: Cp, g
using Plots

include("mlm_solve_funcs.jl")
include("plot_transient_solution.jl")

par = climatology();
par.etype = enBal();
par.ftype = varFlux();
par.rtype = varRad();
par.fttype = fixEIS();

for regime in ("Sc", "Cu")
    for init in (0.2,1)
        println(regime, " ", init)

        if regime == "Sc"
            par.SST0 = 292; # (K)
            par.D = 5.5e-6; # (1/s)
            par.V = 10
        elseif regime == "Cu"
            par.SST0 = 300; # (K)
            par.D = 3.0e-6; # (1/s)
            par.V = 12
        end
        # par.V = 10 # m/s
        par.RHft = 0.2 # (-)
        par.CO2 = 400 # (ppm)
        par.EIS = 10 # (K)

        # run to equilibrium with fixed SST
        dt = 24.0;
        tmax = 100.0;
        par.stype = fixSST();
        u0, sol = run_mlm(par, init=init, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
        uf = sol.u[end];
        println(uf)
        # OHU = calc_OHU(uf, par, incloud_LWP(uf, calc_LCL(uf)), par.stype);
        # println(OHU)

        # # run to equilibrium with prog SST
        # par.stype = varSST();
        # par.Hw = 0.4; # Hw = 0.1 --> τSST = 5 days
        # par.OHU = OHU;
        # # par.EIS = 10*init;
        # u0, sol = run_mlm(par, init=init, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
        # uf = sol.u[end];
        # println(uf)
        # println()

        local path = "experiments/figures/init_condition/";
        mkpath(path);
        filename = path*regime*"_init"*string(init)*".png";
        plot_sol(sol, filename);
    end
end
