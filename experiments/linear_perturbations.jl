# import necessary packages
push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using MixedLayerModel
using MixedLayerModel: Cp, g
using Plots

include("mlm_solve_funcs.jl")

# # stratocumulus
# regime = "Sc"
# init = 1
# SST0 = 290;
# D0 = 5e-6;
# V0 = 10.0;
# RHft0 = 0.2;
# EIS0 = 10.0;
# CO20 = 400.0;
# ylimCF = (70,80)
# ylimLWP = (75,175)

# cumulus
regime = "Cu"
init = 0
SST0 = 300;
D0 = 2.5e-6;
V0 = 8.0;
RHft0 = 0.2;
EIS0 = 5.0;
CO20 = 400.0;
ylimCF = (5,20)
ylimLWP = (50,150)

N = 11 # 3 or 9
_p2 = collect(range(98,102,N)) ./ 100;
_p10 = collect(range(90,110,N)) ./ 100;
_p20 = collect(range(80,120,N)) ./ 100;
_p50 = collect(range(50,150,N)) ./ 100;

X0_arr = [SST0, V0, D0, RHft0, EIS0, CO20];
Xlabel_arr = ["SST [K]", "V [m s\$^{-1}\$]", "D [10\$^{-6}\$ s\$^{-1}\$]", "RH\$_+\$ [%]", "EIS [K]", "CO\$_2\$ [ppm]"];
perturb = [_p2, _p10, _p10, _p20, _p20, _p50];

CF_arr = zeros(length(X0_arr), N);
LWP_arr = zeros(length(X0_arr), N);
zi_arr = zeros(length(X0_arr), N);
zb_arr = zeros(length(X0_arr), N);

for (j,X0) in enumerate(X0_arr)
    local X_arr = X0 .* perturb[j];

    # run simulation looping over X
    for (i,Xi) in enumerate(X_arr)
        # set up parameters
        local par = climatology();
        par.rtype = varRad();
        par.stype = fixSST();
        par.ftype = varFlux();
        par.fttype = fixEIS();
        par.etype = enBal();
        par.dTdz = -6e-3; # K/km

        par.SST0 = (j == 1) ? Xi : SST0 # (K)
        par.V = (j == 2) ? Xi : V0 # m/s
        par.D = (j == 3) ? Xi : D0 # (1/s)
        par.RHft = (j == 4) ? Xi : RHft0 # (-)
        par.EIS = (j == 5) ? Xi : EIS0 # (K)
        par.CO2 = (j == 6) ? Xi : CO20 # (ppm)

        # run to equilibrium
        dt = 24.0;
        tmax = 50.0;
        println(Xi)
        u0, sol = run_mlm(par, init=init, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));

        # extract end state
        uf = sol.u[end];
        println(uf)
        zi, sM, qM, SST, CFi = uf;
        zb = calc_LCL(uf);
        LWPi = incloud_LWP(uf, zb);
        CF_arr[j,i] = CFi;
        LWP_arr[j,i] = LWPi;
        zi_arr[j,] = zi;
        zb_arr[j,i] = zb;
    end
end

# print results
println(CF_arr)
println(LWP_arr*1e3)

ENV["GKSwstype"]="nul"

# plot and save results
j=1
p1 = plot(X0_arr[j] .* perturb[j], CF_arr[j,:]*100, lw=2, label="",
        color=:black,  yguidefontcolor=:black, guidefontsize=16,
        ylim=ylimCF, ylabel="Cloud fraction [%]",
        left_margin=15Plots.mm, right_margin=0Plots.mm,
        bottom_margin=15Plots.mm, top_margin=5Plots.mm)
plot!(twinx(), X0_arr[j] .* perturb[j], LWP_arr[j,:]*1e3, lw=2, label="",
        color=:blue, yguidefontcolor=:blue, ytickfontcolor=:blue,
        ylim=ylimLWP, yformatter=_->"", guidefontsize=16,
        xlabel=Xlabel_arr[j])
j=2
p2 = plot(X0_arr[j] .* perturb[j], CF_arr[j,:]*100, lw=2, label="",
        color=:black,  yguidefontcolor=:black, guidefontsize=16,
        ylim=ylimCF, yformatter=_->"",
        left_margin=0Plots.mm, right_margin=0Plots.mm,
        bottom_margin=15Plots.mm, top_margin=5Plots.mm)
plot!(twinx(), X0_arr[j] .* perturb[j], LWP_arr[j,:]*1e3, lw=2, label="",
        color=:blue, yguidefontcolor=:blue, ytickfontcolor=:blue,
        ylim=ylimLWP, yformatter=_->"", guidefontsize=16,
        xlabel=Xlabel_arr[j])
j=3
p3 = plot(1e6 .* X0_arr[j] .* perturb[j], CF_arr[j,:]*100, lw=2, label="",
        color=:black,  yguidefontcolor=:black, guidefontsize=16,
        ylim=ylimCF, yformatter=_->"",
        left_margin=0Plots.mm, right_margin=30Plots.mm,
        bottom_margin=15Plots.mm, top_margin=5Plots.mm)
plot!(twinx(), 1e6 .* X0_arr[j] .* perturb[j], LWP_arr[j,:]*1e3, lw=2, label="",
        color=:blue, yguidefontcolor=:blue, ytickfontcolor=:blue,
        ylim=ylimLWP, ylabel="In-cloud LWP [g/m\$^2\$]", guidefontsize=16,
        xlabel=Xlabel_arr[j])
j=4
p4 = plot(1e2 .* X0_arr[j] .* perturb[j], CF_arr[j,:]*100, lw=2, label="",
        color=:black,  yguidefontcolor=:black, guidefontsize=16,
        ylim=ylimCF, ylabel="Cloud fraction [%]",
        left_margin=15Plots.mm, right_margin=0Plots.mm,
        bottom_margin=15Plots.mm, top_margin=5Plots.mm)
plot!(twinx(), 1e2 .* X0_arr[j] .* perturb[j], LWP_arr[j,:]*1e3, lw=2, label="",
        color=:blue, yguidefontcolor=:blue, ytickfontcolor=:blue,
        ylim=ylimLWP, yformatter=_->"", guidefontsize=16,
        xlabel=Xlabel_arr[j])
j=5
p5 = plot(X0_arr[j] .* perturb[j], CF_arr[j,:]*100, lw=2, label="",
        color=:black,  yguidefontcolor=:black, guidefontsize=16,
        ylim=ylimCF, yformatter=_->"",
        left_margin=0Plots.mm, right_margin=0Plots.mm,
        bottom_margin=15Plots.mm, top_margin=5Plots.mm)
plot!(twinx(), X0_arr[j] .* perturb[j], LWP_arr[j,:]*1e3, lw=2, label="",
        color=:blue, yguidefontcolor=:blue, ytickfontcolor=:blue,
        ylim=ylimLWP, yformatter=_->"", guidefontsize=16,
        xlabel=Xlabel_arr[j])
j=6
p6 = plot(X0_arr[j] .* perturb[j], CF_arr[j,:]*100, lw=2, label="",
        color=:black,  yguidefontcolor=:black, guidefontsize=16,
        ylim=ylimCF, yformatter=_->"",
        left_margin=0Plots.mm, right_margin=30Plots.mm,
        bottom_margin=15Plots.mm, top_margin=5Plots.mm)
plot!(twinx(), X0_arr[j] .* perturb[j], LWP_arr[j,:]*1e3, lw=2, label="",
        color=:blue, yguidefontcolor=:blue, ytickfontcolor=:blue,
        ylim=ylimLWP, ylabel="In-cloud LWP [g/m\$^2\$]", guidefontsize=16,
        xlabel=Xlabel_arr[j])
plot(p1, p2, p3, p4, p5, p6, layout=(2,3), link=:y, size=(1200,600), dpi=200);

path = "experiments/figures/linear_perturb/";
mkpath(path);
savefig(path*"linear_perturbations_around_"*regime*"_state.png");
