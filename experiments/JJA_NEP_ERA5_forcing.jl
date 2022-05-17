push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using MixedLayerModel
using FileIO
using NCDatasets
using Plots
using Statistics

include("mlm_solve_funcs.jl")

# set up MLM params
par = climatology();
par.rtype = varRad();
par.stype = fixSST();
par.fttype = fixedFT();
par.etype = enBal();
par.ftype = varFlux();
# par.Gamma_q = -2e-6; # (kg/kg/m)
# par.Gamma_s = 5e-3; # (K/m)
par.Gamma_q = 0.0; # (kg/kg/m)
par.Gamma_s = 0.0; # (K/m)

# load boundary conditions from file
file = "experiments/JJA_NEP_boundary_conditions_from_reanalysis.nc";
ds = Dataset(file, "r");

#println(ds)
lon = zeros(length(ds["lat"]))
sst_ss = zeros(length(ds["lat"]))
real_cf = zeros(length(ds["lat"]))

zi_ss = zeros(length(ds["lat"]))
qtM_ss = zeros(length(ds["lat"]))
cf_ss = zeros(length(ds["lat"]))
zb_ss = zeros(length(ds["lat"]))
lwp_ss = zeros(length(ds["lat"]))

for (j,lat) in enumerate(ds["lat"])
    if j < 100
        i = 2j
        j = length(ds["lat"])+1 - j
        i = length(ds["lon"])+1 - i

        par.SST0 = ds["sst"][i,j];
        par.V = 10.0;
        par.D = 6e-6;
        # par.D = 6e-6 - 3e-6*((ds["sst"][i,j]-290)/290);
        # par.LHF = -1*ds["lhf"][i,j];
        # par.SHF = -1*ds["shf"][i,j];
        # par.sft0 = ds["s+"][i,j];
        # par.RHft = 1e-2*ds["RH+"][i,j];

        println(ds["s+"][i,j], "\t", 1e-2*ds["RH+"][i,j];)

        par.sft0 = 310; #ds["sst"][i,j] + 15;
        par.RHft = 0.2;

        dt, tmax = 12, 40;
        u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
        uf = sol.u[end];
        du = zeros(5);
        mlm(du, uf, par, 0.0);

        println(j)
        println(uf)
        println()
        lon[j] = ds["lon"][i];
        zi_ss[j] = uf[1];
        qtM_ss[j] = uf[3];
        cf_ss[j] = uf[5];
        zb_ss[j] = calc_LCL(uf);
        lwp_ss[j] = incloud_LWP(uf, zb_ss[j]);
        sst_ss[j] = par.SST0;
        real_cf[j] = ds["low"][i,j] / 100;

        # t = sol.t / 3600.0 / 24.0;
        # zi = getindex.(sol.u,1);
        # hM = getindex.(sol.u,2) * 1e-3;
        # qtM = getindex.(sol.u,3) * 1e3;
        # sst = getindex.(sol.u,4);
        # cf = getindex.(sol.u,5);
        # S = zeros(length(t));
        # LHF = zeros(length(t));
        # zb = zeros(length(t));
        # ΔR = zeros(length(t));
        # LWP = zeros(length(t));
        # for (i,si) in enumerate(S)
        #     zb[i] = calc_LCL(sol.u[i]);
        #     LWP[i] = incloud_LWP(sol.u[i], zb[i]);
        #     S[i] = calc_S(sol.u[i], par, zb[i], LWP[i]);
        #     LHF[i] = calc_LHF(sol.u[i], par);
        #     ΔR[i] = calc_cloudtop_RAD(sol.u[i], par, LWP[i], par.rtype);
        # end 
        # p = plot(size=(1200,800), layout=(5,2), dpi=200, left_margin = 5Plots.mm);
        # plot!(t, zi, marker="o-", legend=false, subplot=1, ylabel="zi, zb [m]");
        # plot!(t, zb, marker="o-", legend=false, subplot=1);
        # plot!(t, hM, marker="o-", legend=false, subplot=2, ylabel="hM [kJ/kg]"); 
        # plot!(t, qtM, marker="o-", legend=false, subplot=3, ylabel="qtM [g/kg]");
        # plot!(t, sst, marker="o-", legend=false, subplot=4, ylabel="SST [K]");
        # plot!(t, cf * 1e2, marker="o-", legend=false, subplot=5, ylabel="CF [%]");
        # plot!(t, LWP * 1e3, marker="o-", legend=false, subplot=6, ylabel="LWP [g/m2]");
        # plot!(t, LHF, marker="o-", legend=false, subplot=7, ylabel="LHF [W/m2]");
        # plot!(t, ΔR, marker="o-", legend=false, subplot=8, ylabel="ΔR [W/m2]");
        # plot!(t, (zi .- zb) ./ zi, marker="o-", legend=false, subplot=9, ylabel="zc/zi [-]", xlabel="time [days]");
        # plot!(t, S, marker="o-", legend=false, subplot=10, ylabel="S [-]", xlabel="time [days]");
        # savefig("experiments/figures/JJA_NEP_MLM_j14_CF01days.png")
        # display(p)
    end
end

p = plot(size=(600,800), layout=(5,1), dpi=200, left_margin = 5Plots.mm, show=true);
plot!(lon, sst_ss, subplot=1, marker=:o, legend=false, ylabel="SST (K)")
plot!(lon, zi_ss, subplot=2, marker=:o)
plot!(lon, zb_ss, subplot=2, marker=:o, legend=false, ylabel="zi, zb (m)")
plot!(lon, lwp_ss*1e3, subplot=3, marker=:o, legend=false, ylabel="LWP (g/m2)")
plot!(lon, cf_ss*100, subplot=4, marker=:o, legend=false, ylabel="CF (%)")
#plot!(lon, real_cf*100, subplot=4, marker=:o)
plot!(lon, qtM_ss*1e3, subplot=5, marker=:o, legend=false, ylabel="qt (g/kg)", xlabel="longitude")
savefig("experiments/figures/JJA_NEP_MLM_transect.png");
display(p)

close(ds)
