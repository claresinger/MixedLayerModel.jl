# push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using MixedLayerModel
using MixedLayerModel: Cp, g
using FileIO
using NCDatasets
using Plots
using Statistics

include("mlm_solve_funcs.jl")

# set up MLM params
par = climatology();
par.rtype = varRad();
par.stype = fixSST();
par.ftype = varFlux();
par.fttype = fixEIS();
par.etype = enBal();

# load boundary conditions from file
file = "experiments/data/box_BCs_JJA_NEP.nc";
ds = Dataset(file, "r");

skip1 = 1
skip2 = 1
lon = ds["lon"][1:skip1:end]
lat = ds["lat"][1:skip2:end]

p = plot(size=(600,300), layout=(2,1), dpi=400, 
        left_margin = 5Plots.mm, bottom_margin = 5Plots.mm);
ct = 0
for (i1,loni) in enumerate(lon)
    for (i2, lati) in enumerate(lat)
        local j1 = (length(lon)+1-i1)*skip1
        local j2 = (length(lat)+1-i2)*skip2
        if (i1 == 12 && i2 == 6) || (i1 == 1 && i2 == 1)
            global ct += 1
            println(i1, ", ", i2)
            println(j1, ", ", j2)

            par.SST0 = ds["sst"][j1,j2];
            par.V = ds["WS"][j1,j2];
            par.D = ds["D500"][j1,j2];
            par.RHft = ds["RH500"][j1,j2];
            par.EIS = ds["EIS"][j1,j2];

            # par.SST0 = ds["sst"][j1,j2];
            # par.V = mean(ds["WS"]);
            # par.D = mean(ds["D500"]);
            # par.RHft = mean(ds["RH500"]);
            # par.EIS = mean(ds["EIS"]);

            # par.SST0 = mean(ds["sst"]);
            # par.V = mean(ds["WS"]);
            # par.D = mean(ds["D500"]);
            # par.RHft = mean(ds["RH500"]);
            # par.EIS = ds["EIS"][j1,j2];

            # par.SST0 = mean(ds["sst"]);
            # par.V = mean(ds["WS"]);
            # par.D = mean(ds["D500"]);
            # par.RHft = ds["RH500"][j1,j2];
            # par.EIS = mean(ds["EIS"]);

            dt, tmax = 24, 50;
            u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
            uf = sol.u[end];
            du = zeros(5);
            mlm(du, uf, par, 0.0);

            println(uf)
            println(du)

            println(lati, " ", loni)
            println(par.SST0, " ", ct)
            t = sol.t / 3600.0 / 24.0;
            zi = getindex.(sol.u,1);
            zb = [calc_LCL(uk) for uk in sol.u];
            cf = getindex.(sol.u,5);
            plot!(t, zi, color=ct, ls=:solid, lw=2, subplot=1, ylabel="zi, zb [m]", ylim=(0,2300));
            plot!(t, zb, color=ct, ls=:dash, lw=2, subplot=1, legend=false);
            plot!(t, cf * 1e2, color=ct, ls=:solid, lw=2, legend=:bottomleft, subplot=2, 
                ylabel="CF [%]", xlabel="Time [days]", ylim=(0,90),
                label="SST="*string(Int(round(par.SST0)))*" (K)");
        end
    end
end

mkpath("experiments/figures/cfmip/")

savefig("experiments/figures/cfmip/JJA_NEP_ScCu_all.png")
# savefig("experiments/figures/cfmip/JJA_NEP_ScCu_onlySST.png")
# savefig("experiments/figures/cfmip/JJA_NEP_ScCu_onlyEIS.png")
# savefig("experiments/figures/cfmip/JJA_NEP_ScCu_onlyRH.png")