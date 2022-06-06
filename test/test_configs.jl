push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using MixedLayerModel
include("../experiments/mlm_solve_funcs.jl")
include("../experiments/plot_transient_solution.jl")

dt = 12.0;
tmax = 40.0;

ENV["GKSwstype"]="nul"
mkpath("figures/");

# test climatology config
par = climatology();
for entrainment in (enBal(), bflux())
    for fluxes in (fixFlux(), varFlux())
        for radiation in (fixRad(), varRad())
            par.etype = entrainment;
            par.ftype = fluxes;
            par.rtype = radiation;
            par.stype = fixSST();
            par.fttype = fixedFT();

            u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
            plot_sol(sol, "figures/climatology"*string(entrainment)
                *string(fluxes)*string(radiation)*"_sol400.png")
            uf = sol.u[end];
            zb = calc_LCL(uf);
            du = zeros(5);
            mlm(du, uf, par, 0.0);

            println(entrainment, fluxes, radiation)
            println(uf)
            @test all(uf .> 0)
            @test zb <= uf[1]
            @test all(du/uf .< 1e-3)
        end
    end
end

# test co2 config
par = upCO2();
for entrainment in (enBal(), bflux())
    for sst in (fixSST(), varSST())
        for freetrop in (co2dep(), twocol())
            par.etype = entrainment;
            par.ftype = varFlux();
            par.rtype = varRad();
            par.stype = sst;
            par.fttype = freetrop;

            u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
            plot_sol(sol, "figures/upCO2"*string(entrainment)
                *string(sst)*string(freetrop)*"_sol400.png")
            uf = sol.u[end];
            zb = calc_LCL(uf);
            du = zeros(5);
            mlm(du, uf, par, 0.0);

            println(entrainment, sst, freetrop)
            println(uf)
            @test all(uf .> 0)
            @test zb <= uf[1]
            @test all(du/uf .< 1e-3)
        end
    end
end
