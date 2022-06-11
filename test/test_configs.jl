include("../experiments/mlm_solve_funcs.jl")
include("../experiments/plot_transient_solution.jl")

using FileIO

makeplot = true

if makeplot
    dt = 6.0;
    tmax = 40.0;
    ENV["GKSwstype"]="nul"
    mkpath("figures/");
else
    dt = 6.0;
    tmax = 1.0;
end

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
            println(entrainment, fluxes, radiation)

            u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
            pathname = "figures/climatology_";
            filename = pathname*string(entrainment)*string(fluxes)*string(radiation)*"_sol400.png";
            makeplot ? plot_sol(sol, filename) : println("no plots")
            filename = pathname*string(entrainment)*string(fluxes)*string(radiation)*"_sol400.jld2";
            makeplot ? save(filename, Dict("sol" => sol)) : println("no save")
            uf = sol.u[end];
            zb = calc_LCL(uf);
            du = zeros(5);
            mlm(du, uf, par, 0.0);

            @test all(uf .> 0)
            @test zb <= uf[1]
            @test all(du/uf .< 1e-3)

            if makeplot
                pathname = "main_figures/climatology_";
                filename = pathname*string(entrainment)*string(fluxes)*string(radiation)*"_sol400.jld2";
                mainsol = load(filename)["mainsol"]
                println(sol.u[end])
                println(mainsol.u[end])
                @test isapprox(sol.u[end][1], mainsol.u[end][1], rtol = 0.1)
                println()
            end
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
            println(entrainment, sst, freetrop)

            u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
            pathname = "figures/upCO2_";
            filename = pathname*string(entrainment)*string(sst)*string(freetrop)*"_sol400.png";
            makeplot ? plot_sol(sol, filename) : println("no plots")
            filename = pathname*string(entrainment)*string(sst)*string(freetrop)*"_sol400.jld2";
            makeplot ? save(filename, Dict("sol" => sol)) : println("no save")
            uf = sol.u[end];
            zb = calc_LCL(uf);
            du = zeros(5);
            mlm(du, uf, par, 0.0);

            @test all(uf .> 0)
            @test zb <= uf[1]
            @test all(du/uf .< 1e-3)

            if makeplot
                pathname = "main_figures/upCO2_";
                filename = pathname*string(entrainment)*string(sst)*string(freetrop)*"_sol400.jld2";
                mainsol = load(filename)["mainsol"]
                println(sol.u[end])
                println(mainsol.u[end])
                @test isapprox(sol.u[end][1], mainsol.u[end][1], rtol = 0.1)
                println()
            end
        end
    end
end
