push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using MixedLayerModel
using MixedLayerModel: Cp, g
using FileIO
using NCDatasets
using Plots
using Statistics

path = "experiments/figures/remove_cf_radcool/"
mkpath(path)

include("mlm_solve_funcs.jl")

# set up MLM params
par = climatology();
par.rtype = varRad();
par.stype = fixSST();
par.ftype = varFlux();
par.fttype = fixEIS();
par.etype = enBal();

# load boundary conditions from file
file = "experiments/data/transect_BCs_JJA_NEP.nc";
ds = Dataset(file, "r");
# println(ds)

skipi = 1
lon = ds["lon"][1:skipi:end]
N = length(lon)
println(N, " out of ", length(ds["lon"]))
sst_ss = zeros(N)
real_cf = zeros(N)
zi_ss = zeros(N)
qtM_ss = zeros(N)
zb_ss = zeros(N)
lwp_ss = zeros(N)
lhf_ss = zeros(N)
Tsurf_ss = zeros(N)

cf_ss = zeros(N)
cf_ss_min = zeros(N)
cf_ss_max = zeros(N)
cf_ss_onlySST = zeros(N)
cf_ss_onlyWS = zeros(N)
cf_ss_onlyEIS = zeros(N)
cf_ss_onlyD = zeros(N)
cf_ss_onlyRH = zeros(N)

for (i,loni) in enumerate(lon)
    local j = i*skipi
    if i > 0 
        # mean of all parameters, mean cf
        par.SST0 = ds["sst"][j];
        par.V = ds["WS"][j];
        par.D = ds["D500"][j];
        par.RHft = ds["RH500"][j];
        par.EIS = ds["EIS"][j];        
        dt, tmax = 24, 50;
        u0, sol = run_mlm(par, init=1, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
        uf = sol.u[end];
        println(uf);

        # save for file
        zi_ss[i] = uf[1];
        Tsurf_ss[i] = uf[2] / Cp;
        qtM_ss[i] = uf[3];
        cf_ss[i] = uf[5];
        zb_ss[i] = calc_LCL(uf);
        lwp_ss[i] = incloud_LWP(uf, zb_ss[i]);
        lhf_ss[i] = calc_LHF(uf, par);
        sst_ss[i] = par.SST0;
        real_cf[i] = ds["allsc"][j];

        # minimum cf 
        par.SST0 = ds["sst"][j] + ds["sst_std"][j];
        par.V = ds["WS"][j] + ds["WS_std"][j];
        par.D = ds["D500"][j] - ds["D500_std"][j];
        par.RHft = ds["RH500"][j] - ds["RH500_std"][j];
        par.EIS = ds["EIS"][j] - ds["EIS_std"][j];
        u0, sol_min = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
        cf_ss_min[i] = sol_min.u[end][5];

        # maximum cf
        par.SST0 = ds["sst"][j] - ds["sst_std"][j];
        par.V = ds["WS"][j] - ds["WS_std"][j];
        par.D = ds["D500"][j] + ds["D500_std"][j];
        par.RHft = ds["RH500"][j] + ds["RH500_std"][j];
        par.EIS = ds["EIS"][j] + ds["EIS_std"][j];
        u0, sol_max = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
        cf_ss_max[i] = sol_max.u[end][5];

        # set mean values
        ilon = 1:10
        par.SST0 = mean(ds["sst"][ilon]);
        par.V = mean(ds["WS"][ilon]);
        par.D = mean(ds["D500"][ilon]);
        par.RHft = mean(ds["RH500"][ilon]);
        par.EIS = mean(ds["EIS"][ilon]);

        # only SST
        par.SST0 = ds["sst"][j];
        u0, sol_onlySST = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
        cf_ss_onlySST[i] = sol_onlySST.u[end][5];
        par.SST0 = mean(ds["sst"][ilon]);

        # only WS
        par.V = ds["WS"][j];
        u0, sol_onlyWS = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
        cf_ss_onlyWS[i] = sol_onlyWS.u[end][5];
        par.V = mean(ds["WS"][ilon]);

        # only EIS (set and reset)
        par.EIS = ds["EIS"][j];
        u0, sol_onlyEIS = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
        cf_ss_onlyEIS[i] = sol_onlyEIS.u[end][5];
        par.EIS = mean(ds["EIS"][ilon]);

        # only D500
        par.D = ds["D500"][j];
        u0, sol_onlyD = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
        cf_ss_onlyD[i] = sol_onlyD.u[end][5];
        par.D = mean(ds["D500"][ilon]);

        # only RH500
        par.RHft = ds["RH500"][j];
        u0, sol_onlyRH = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
        cf_ss_onlyRH[i] = sol_onlyRH.u[end][5];
        par.RHft = mean(ds["RH500"][ilon]);

        if uf[5] > 0.15 && uf[5] < 0.75
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
            plot!(sol_min.t/3600.0/24.0, getindex.(sol_min.u,5) * 1e2, marker="o-", legend=false, subplot=5);
            plot!(sol_max.t/3600.0/24.0, getindex.(sol_max.u,5) * 1e2, marker="o-", legend=false, subplot=5);
            plot!(t, cf .* LWP * 1e3, marker="o-", legend=false, subplot=6, ylabel="mean LWP [g/m2]");
            plot!(t, LHF, marker="o-", legend=false, subplot=7, ylabel="(S/L)HF [W/m2]");
            plot!(t, SHF, marker="o-", legend=false, subplot=7);
            plot!(t, ΔR, marker="o-", legend=false, subplot=8, ylabel="ΔR [W/m2]");
            plot!(t, (zi .- zb) ./ zi, marker="o-", legend=false, subplot=9, ylabel="zc/zi [-]", xlabel="time [days]");
            plot!(t, S, marker="o-", legend=false, subplot=10, ylabel="S [-]", xlabel="time [days]");
            plot!(t, Δsv / Cp, marker="o-", legend=false, subplot=11, ylabel="Δsv/Cp [K]");
            plot!(t, ent*1e3, marker="o-", legend=false, subplot=12, ylabel="we [mm/s]")
            savefig(path*"transect_JJA_NEP_i"*string(i)*".png")
        end
    end
end

println()
println(cf_ss_min)
println(cf_ss)
println(cf_ss_max)

# p = plot(size=(600,400), layout=(2,1), dpi=300, show=true,
#     left_margin = 5Plots.mm, right_margin = 15Plots.mm);
# # plot observed SST
# plot!(lon, ds["sst"][1:skipi:end], subplot=1, legend=false,
#     lw=2, color=:magenta, yguidefontcolor=:magenta, ytickfontcolor=:magenta,
#     marker=:circle, msw=0,
#     ylabel="SST [K]", ribbon=ds["sst_std"][1:skipi:end], fillalpha=0.3)
# # plot observed EIS
# plot!(twinx(), lon, ds["EIS"][1:skipi:end], lw=2, legend=false,
#         color=:green, yguidefontcolor=:green, ytickfontcolor=:green,
#         marker=:circle, msw=0,
#         ylabel="EIS [K]", ribbon=ds["EIS_std"][1:skipi:end], fillalpha=0.3)

# # plot observed CF
# plot!(lon, ds["allsc"][1:skipi:end]*100, subplot=2, legend=:topleft,
#     lw=2, color=:black, marker=:circle, msw=0, label="Obs",
#     ribbon=ds["allsc_std"][1:skipi:end]*100, fillalpha=0.3)

# # predicted CF with range
# plot!(lon, cf_ss*100, subplot=2, lw=2, ylim=(0,90), 
#     ylabel="CF [%]", color=1, marker=:circle, msw=0, label="ERA5 BCs",
#     ribbon=((cf_ss-cf_ss_min)*100, (cf_ss_max-cf_ss)*100), fillalpha=0.3)

# tkloc, tkstr = xticks(p)[1]
# plot!(xticks=(tkloc, chop.(tkstr,head=1,tail=0).*" °W"))

# savefig(path*"JJA_NEP_transect.png");

# ### only SST and only EIS too
# p = plot(size=(600,400), layout=(2,1), dpi=300, show=true,
#     left_margin = 5Plots.mm, right_margin = 15Plots.mm);
# # plot observed SST
# plot!(lon, ds["sst"][1:skipi:end], subplot=1, legend=false,
#     lw=2, color=:magenta, yguidefontcolor=:magenta, ytickfontcolor=:magenta,
#     marker=:circle, msw=0,
#     ylabel="SST [K]", ribbon=ds["sst_std"][1:skipi:end], fillalpha=0.3)
# # plot observed EIS
# plot!(twinx(), lon, ds["EIS"][1:skipi:end], lw=2, legend=false,
#         color=:green, yguidefontcolor=:green, ytickfontcolor=:green,
#         marker=:circle, msw=0,
#         ylabel="EIS [K]", ribbon=ds["EIS_std"][1:skipi:end], fillalpha=0.3)

# # plot observed CF
# plot!(lon, ds["allsc"][1:skipi:end]*100, subplot=2, legend=:topleft,
#     lw=2, color=:black, marker=:circle, msw=0, label="Obs",
#     ribbon=ds["allsc_std"][1:skipi:end]*100, fillalpha=0.3)

# # predicted CF with range
# plot!(lon, cf_ss*100, subplot=2, lw=2, ylim=(0,90), 
#     ylabel="CF [%]", color=1, marker=:circle, msw=0, label="ERA5 BCs",
#     ribbon=((cf_ss-cf_ss_min)*100, (cf_ss_max-cf_ss)*100), fillalpha=0.3)

# # predicted CF from only varying SST or EIS
# plot!(lon, cf_ss_onlySST*100, subplot=2, lw=2, marker=:circle, msw=0, color=:magenta, label="only SST")
# plot!(lon, cf_ss_onlyEIS*100, subplot=2, lw=2, marker=:circle, msw=0, color=:green, label="only EIS")

# tkloc, tkstr = xticks(p)[1]
# plot!(xticks=(tkloc, chop.(tkstr,head=1,tail=0).*" °W"))

# savefig(path*"JJA_NEP_transect_1var.png");

### only SST and only EIS too
p = plot(size=(600,200), layout=(1,1), dpi=300, show=true,
    left_margin = 5Plots.mm, right_margin = 2Plots.mm, palette = :tab10);
# plot observed CF
plot!(lon, ds["allsc"][1:skipi:end]*100, subplot=1, legend=:topleft,
    lw=2, color=:black, marker=:circle, msw=0, label="Obs",
    ribbon=ds["allsc_std"][1:skipi:end]*100, fillalpha=0.3)
# predicted CF with range
plot!(lon, cf_ss*100, subplot=1, lw=2, ylim=(0,90), 
    ylabel="CF [%]", color=:magenta, marker=:circle, msw=0, label="ERA5 BCs",
    ribbon=((cf_ss-cf_ss_min)*100, (cf_ss_max-cf_ss)*100), fillalpha=0.3)
# predicted CF from only varying on CCF
plot!(lon, cf_ss_onlySST*100, subplot=1, lw=2, marker=:circle, msw=0, color=1, label="SST")
plot!(lon, cf_ss_onlyWS*100, subplot=1, lw=2, marker=:circle, msw=0, color=2, label="WS")
plot!(lon, cf_ss_onlyEIS*100, subplot=1, lw=2, marker=:circle, msw=0, color=3, label="EIS")
plot!(lon, cf_ss_onlyD*100, subplot=1, lw=2, marker=:circle, msw=0, color=4, label="D\$_{500}\$")
plot!(lon, cf_ss_onlyRH*100, subplot=1, lw=2, marker=:circle, msw=0, color=5, label="RH\$_{500}\$")

tkloc, tkstr = xticks(p)[1];
plot!(xticks=(tkloc, chop.(tkstr,head=1,tail=0).*" °W"));
plot!(legend_position=:outertopright);

savefig(path*"JJA_NEP_transect_1var.png");

close(ds)


# create output netcdf file
ds = Dataset(path*"transect_output.nc","c")

# Define the dimension "lon" of size N
defDim(ds,"lon",N)
defVar(ds,"lon",lon,("lon",))

# Define the variables
v = defVar(ds,"zi",zi_ss,("lon",))
v.attrib["units"] = "m"
v.attrib["long_name"] = "cloud top height"

v = defVar(ds,"zb",zb_ss,("lon",))
v.attrib["units"] = "m"
v.attrib["long_name"] = "cloud base height"

v = defVar(ds,"LHF",lhf_ss,("lon",))
v.attrib["units"] = "W/m2"
v.attrib["long_name"] = "surface latent heat flux"

v = defVar(ds,"Tsurf",Tsurf_ss,("lon",))
v.attrib["units"] = "K"
v.attrib["long_name"] = "surface air temperature"

v = defVar(ds,"sst",sst_ss,("lon",))
v.attrib["units"] = "K"
v.attrib["long_name"] = "sea surface temperature"

v = defVar(ds,"cf",cf_ss,("lon",))
v.attrib["units"] = "-"
v.attrib["long_name"] = "cloud fraction"

v = defVar(ds,"cf_min",cf_ss_min,("lon",))
v.attrib["units"] = "-"
v.attrib["long_name"] = "min cloud fraction"

v = defVar(ds,"cf_max",cf_ss_max,("lon",))
v.attrib["units"] = "-"
v.attrib["long_name"] = "max cloud fraction"

v = defVar(ds,"cf_onlySST",cf_ss_onlySST,("lon",))
v.attrib["units"] = "-"
v.attrib["long_name"] = "cloud fraction from only vary SST"

v = defVar(ds,"cf_onlyEIS",cf_ss_onlyEIS,("lon",))
v.attrib["units"] = "-"
v.attrib["long_name"] = "cloud fraction from only vary EIS"

v = defVar(ds,"lwp",lwp_ss,("lon",))
v.attrib["units"] = "kg/m2"
v.attrib["long_name"] = "in-cloud liquid water path"

v = defVar(ds,"obs_cf",real_cf,("lon",))
v.attrib["units"] = "-"
v.attrib["long_name"] = "observed cloud fraction (CASCCAD)"

# print(ds)
close(ds)
