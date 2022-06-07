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
file = "experiments/transect_BCs_JJA_NEP.nc";
ds = Dataset(file, "r");
# println(ds)

lon = ds["lon"][1:5:end]
N = length(lon)
sst_ss = zeros(N)
real_cf = zeros(N)
zi_ss = zeros(N)
qtM_ss = zeros(N)
cf_ss = zeros(N)
zb_ss = zeros(N)
lwp_ss = zeros(N)

for (i,loni) in enumerate(lon)
    j = i*5
    if i > 0
        par.SST0 = ds["sst"][j];
        par.V = 10.0;
        # par.LHF = ds["LHF"][j];
        # par.SHF = ds["SHF"][j];
        par.RHft = ds["RH500"][j];
        par.D = ds["D500"][j];
        # par.Tft = par.SST0 + 10;
        # par.Tft = ds["T500"][j] + 30;
        # par.sft0 = Cp*ds["T500"][j] + g*1000;
        # par.RHft = 0.2;
        # par.D = 6e-6;
        par.sft0 = Cp*par.SST0 + g*1000;
        par.Gamma_q = 0.0;
        par.Gamma_s = (Cp*-5e-3) + g; # (K/m)

        # println(loni)
        # println(par)
        
        dt, tmax = 12, 40;
        u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
        uf = sol.u[end];
        du = zeros(5);
        mlm(du, uf, par, 0.0);

        println(j)
        println(uf)
        println()
        zi_ss[i] = uf[1];
        qtM_ss[i] = uf[3];
        cf_ss[i] = uf[5];
        zb_ss[i] = calc_LCL(uf);
        lwp_ss[i] = incloud_LWP(uf, zb_ss[i]);
        sst_ss[i] = par.SST0;
        real_cf[i] = ds["low"][i];

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
        # for (k,si) in enumerate(S)
        #     zb[k] = calc_LCL(sol.u[k]);
        #     LWP[k] = incloud_LWP(sol.u[k], zb[k]);
        #     S[k] = calc_S(sol.u[k], par, zb[k], LWP[k]);
        #     LHF[k] = calc_LHF(sol.u[k], par);
        #     ΔR[k] = calc_cloudtop_RAD(sol.u[k], par, LWP[k], par.rtype);
        # end 
        # p = plot(size=(1200,800), layout=(5,2), dpi=200, left_margin = 5Plots.mm);
        # plot!(t, zi, marker="o-", legend=false, subplot=1, ylabel="zi, zb [m]");
        # plot!(t, zb, marker="o-", legend=false, subplot=1);
        # plot!(t, hM, marker="o-", legend=false, subplot=2, ylabel="hM [kJ/kg]"); 
        # plot!(t, qtM, marker="o-", legend=false, subplot=3, ylabel="qtM [g/kg]");
        # plot!(t, sst, marker="o-", legend=false, subplot=4, ylabel="SST [K]");
        # plot!(t, cf * 1e2, marker="o-", legend=false, subplot=5, ylabel="CF [%]");
        # plot!(t, cf .* LWP * 1e3, marker="o-", legend=false, subplot=6, ylabel="mean LWP [g/m2]");
        # plot!(t, LHF, marker="o-", legend=false, subplot=7, ylabel="LHF [W/m2]");
        # plot!(t, ΔR, marker="o-", legend=false, subplot=8, ylabel="ΔR [W/m2]");
        # plot!(t, (zi .- zb) ./ zi, marker="o-", legend=false, subplot=9, ylabel="zc/zi [-]", xlabel="time [days]");
        # plot!(t, S, marker="o-", legend=false, subplot=10, ylabel="S [-]", xlabel="time [days]");
        # savefig("experiments/figures/transect_JJA_NEP_i2.png")
        # display(p)
    end
end

p = plot(size=(600,600), layout=(3,1), dpi=200, left_margin = 5Plots.mm, show=true);
plot!(lon, sst_ss, subplot=1, marker=:o, color=:black, legend=false, ylabel="SST (K)")
plot!(lon, zi_ss, subplot=2, marker=:o)
plot!(lon, zb_ss, subplot=2, marker=:o, legend=false, ylabel="zb, zi (m)")
plot!(lon, cf_ss*100, subplot=3, marker=:o, legend=false, ylabel="CF (%)", xlabel="longitude")
plot!(lon, real_cf*100, subplot=3, marker=:o, color=:black)
#plot!(lon, lwp_ss*1e3, subplot=4, marker=:o, legend=false, ylabel="LWP (g/m2)")
#plot!(lon, qtM_ss*1e3, subplot=5, marker=:o, legend=false, ylabel="qt (g/kg)", xlabel="longitude")
file = "experiments/transect_BCs_JJA_NEP.nc";
savefig(replace(file, "/"=>"/figures/", "_BCs_"=>"_", ".nc"=>".png"));
# display(p)

close(ds)
