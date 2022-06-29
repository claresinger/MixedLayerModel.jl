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
par.stype = varSST();
par.ftype = fixFlux();
par.fttype = fixedFT();
par.etype = enBal();

# load boundary conditions from file
file = "experiments/data/transect_BCs_JJA_NEP.nc";
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
    local j = i*5
    if i == 1
        par.SST0 = ds["sst"][j];

        par.OHU = -5;
        par.LHF = ds["LHF"][j];
        par.SHF = 1.0;
        # par.SHF = ds["SHF"][j];

        # par.D = 6e-6;
        # par.RHft = 0.2;
        
        par.D = ds["D500"][j];
        par.RHft = ds["RH500"][j];

        EIS = 10.0; # (K)
        par.sft0 = Cp*(par.SST0 + EIS); # 10 (K) jump
        par.Gamma_s = (Cp*-5e-3) + g; # (K/m)

        # println(loni)
        # println(par)
        
        dt, tmax = 24, 100;
        u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
        uf = sol.u[end];
        du = zeros(5);
        mlm(du, uf, par, 0.0);

        println(j)
        println(uf)
        println(calc_OHU(uf, par, incloud_LWP(uf, calc_LCL(uf)), par.stype))
        println()
        zi_ss[i] = uf[1];
        qtM_ss[i] = uf[3];
        cf_ss[i] = uf[5];
        zb_ss[i] = calc_LCL(uf);
        lwp_ss[i] = incloud_LWP(uf, zb_ss[i]);
        sst_ss[i] = par.SST0;
        real_cf[i] = ds["allsc"][j];

        if i in (1, 10, 20)
            t = sol.t / 3600.0 / 24.0;
            zi = getindex.(sol.u,1);
            sM = getindex.(sol.u,2) * 1e-3;
            qtM = getindex.(sol.u,3) * 1e3;
            sst = getindex.(sol.u,4);
            cf = getindex.(sol.u,5);
            S = zeros(length(t));
            LHF = zeros(length(t));
            SHF = zeros(length(t));
            zb = zeros(length(t));
            ΔR = zeros(length(t));
            LWP = zeros(length(t));
            Δsv = zeros(length(t));
            ent = zeros(length(t));
            for (k,si) in enumerate(S)
                zb[k] = calc_LCL(sol.u[k]);
                LWP[k] = incloud_LWP(sol.u[k], zb[k]);
                S[k] = calc_S(sol.u[k], par, zb[k], LWP[k]);
                LHF[k] = calc_LHF(sol.u[k], par);
                SHF[k] = calc_SHF(sol.u[k], par);
                ΔR[k] = calc_cloudtop_RAD(sol.u[k], par, LWP[k], par.rtype);
                Δsv[k] = sv_jump(sol.u[k], par, LWP[k]);
                ent[k] = we(sol.u[k], par, zb[k], LWP[k], par.etype);
            end 
            local p = plot(size=(1200,800), layout=(6,2), dpi=200, left_margin = 5Plots.mm);
            plot!(t, zi, marker="o-", legend=false, subplot=1, ylabel="zi, zb [m]");
            plot!(t, zb, marker="o-", legend=false, subplot=1);
            plot!(t, sM, marker="o-", legend=false, subplot=2, ylabel="sM [kJ/kg]"); 
            plot!(t, qtM, marker="o-", legend=false, subplot=3, ylabel="qtM [g/kg]");
            plot!(t, sst, marker="o-", legend=false, subplot=4, ylabel="SST [K]");
            plot!(t, cf * 1e2, marker="o-", legend=false, subplot=5, ylabel="CF [%]");
            plot!(t, cf .* LWP * 1e3, marker="o-", legend=false, subplot=6, ylabel="mean LWP [g/m2]");
            plot!(t, LHF, marker="o-", legend=false, subplot=7, ylabel="(S/L)HF [W/m2]");
            plot!(t, SHF, marker="o-", legend=false, subplot=7);
            plot!(t, ΔR, marker="o-", legend=false, subplot=8, ylabel="ΔR [W/m2]");
            plot!(t, (zi .- zb) ./ zi, marker="o-", legend=false, subplot=9, ylabel="zc/zi [-]", xlabel="time [days]");
            plot!(t, S, marker="o-", legend=false, subplot=10, ylabel="S [-]", xlabel="time [days]");
            plot!(t, Δsv / Cp, marker="o-", legend=false, subplot=11, ylabel="Δsv/Cp [K]");
            plot!(t, ent*1e3, marker="o-", legend=false, subplot=12, ylabel="we [mm/s]")
            savefig("experiments/figures/transect_JJA_NEP_i"*string(i)*".png")
        end
    end
end

# p = plot(size=(600,600), layout=(3,1), dpi=200, left_margin = 5Plots.mm, show=true);
# plot!(lon, sst_ss, subplot=1, marker=:o, color=:black, legend=false, ylabel="SST (K)")
# plot!(lon, zi_ss, subplot=2, marker=:o)
# plot!(lon, zb_ss, subplot=2, marker=:o, legend=false, ylabel="zb, zi (m)")
# plot!(lon, cf_ss*100, subplot=3, marker=:o, legend=false, ylabel="CF (%)", xlabel="longitude")
# plot!(lon, real_cf*100, subplot=3, marker=:o, color=:black)
# file = "experiments/data/transect_BCs_JJA_NEP.nc";
# savefig(replace(file, "data/"=>"figures/", "_BCs_"=>"_", ".nc"=>".png"));

close(ds)
