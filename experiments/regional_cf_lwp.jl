push!(LOAD_PATH, joinpath(@__DIR__, ".."))

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
par.fttype = fixedFT();
par.etype = enBal();

# load boundary conditions from file
file = "experiments/data/box_BCs_JJA_NEP.nc";
ds = Dataset(file, "r");
# println(ds)

skip1 = 1
skip2 = 1
lon = ds["lon"][1:skip1:end]
lat = ds["lat"][1:skip2:end]
N, M = length(lon), length(lat)
println(N, " ", M)
sst_ss = zeros(N, M)
real_cf = zeros(N, M)
zi_ss = zeros(N, M)
qtM_ss = zeros(N, M)
cf_ss = zeros(N, M)
zb_ss = zeros(N, M)
lwp_ss = zeros(N, M)

for (i1,loni) in enumerate(lon)
    for (i2, lati) in enumerate(lat)
        local j1 = i1*skip1
        local j2 = i2*skip2
        if i1 == 12 && i2 == 6
            println(i1, ", ", i2)
            println(j1, ", ", j2)

            par.SST0 = ds["sst"][j1,j2];
            println(par.SST0)

            # par.OHU = -5;
            # par.LHF = ds["LHF"][j1,j2];
            # par.SHF = 1.0;
            # par.SHF = ds["SHF"][j1,j2];

            # par.D = 6e-6;
            # par.RHft = 0.2;
            
            par.D = ds["D500"][j1,j2];
            par.RHft = ds["RH500"][j1,j2];

            EIS = 10.0; # (K)
            par.sft0 = Cp*(par.SST0 + EIS); # 10 (K) jump
            par.Gamma_s = (Cp*-5e-3) + g; # (K/m)

            dt, tmax = 24, 50;
            u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
            uf = sol.u[end];
            du = zeros(5);
            mlm(du, uf, par, 0.0);

            println(uf)
            println(du)
            zi_ss[i1, i2] = uf[1];
            qtM_ss[i1, i2] = uf[3];
            cf_ss[i1, i2] = uf[5];
            zb_ss[i1, i2] = calc_LCL(uf);
            lwp_ss[i1, i2] = incloud_LWP(uf, zb_ss[i1, i2]);
            sst_ss[i1, i2] = par.SST0;
            real_cf[i1, i2] = ds["allsc"][j1, j2];

            if i1 == 12 && i2 == 6
                t = sol.t / 3600.0 / 24.0;
                zi = getindex.(sol.u,1);
                zb = [calc_LCL(uk) for uk in sol.u];
                cf = getindex.(sol.u,5);
                local p = plot(size=(600,300), layout=(2,1), dpi=200,
                    left_margin = 5Plots.mm, bottom_margin = 5Plots.mm);
                plot!(t, zi, marker="o-", legend=false, subplot=1, ylabel="zi, zb [m]");
                plot!(t, zb, marker="o-", legend=false, subplot=1);
                plot!(t, cf * 1e2, marker="o-", legend=false, subplot=2, 
                    ylabel="CF [%]", xlabel="Time [days]");
                savefig("experiments/figures/box_JJA_NEP_"*string(i1)*"+"*string(i2)*".png")
            end
        end
    end
end

p1 = contourf(lon, lat, cf_ss, c=:Greens, clims=(0,1), 
    title="NEP JJA", xlabel="lon", ylabel="lat",
    colorbar_title="Cloud fraction")
p2 = contourf(lon, lat, lwp_ss*1e3, c=:Blues, clims=(0,1000), 
    title="NEP JJA", xlabel="lon", ylabel="lat",
    colorbar_title="In-cloud LWP [g/m2]")
p3 = contourf(lon, lat, lwp_ss.*cf_ss*1e3, c=:Blues, clims=(0,1000), 
    title="NEP JJA", xlabel="lon", ylabel="lat",
    colorbar_title="All-sky LWP [g/m2]")

plot(p1, p2, p3, layout=(3,1), size=(800,1000), dpi=200)
savefig("experiments/figures/box_JJA_NEP.png")

close(ds)
