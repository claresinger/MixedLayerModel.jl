# import necessary packages
push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using MixedLayerModel
using MixedLayerModel: Cp, g
using Plots

include("mlm_solve_funcs.jl")

SST0 = 295;
D0 = 5.5e-6;
V0 = 10.0;
RHft0 = 0.2;
EIS0 = 10.0;
CO20 = 400.0;

perturb_SST = collect(95 : 1.25 : 105) ./ 100;
perturb = collect(80 : 5 : 120) ./ 100;
# perturb_SST = collect(95 : 5 : 105) ./ 100;
# perturb = collect(80 : 20 : 120) ./ 100;

X0_arr = [SST0, V0, D0, RHft0, EIS0, CO20];
Xlabel_arr = ["SST (K)", "V (m/s)", "D x 10\$^6\$ (1/s)", "RH\$^{ft}\$", "EIS (K)", "CO\$_2\$ (ppm)"];

CF_arr = zeros(length(X0_arr), length(perturb));
LWP_arr = zeros(length(X0_arr), length(perturb));
zi_arr = zeros(length(X0_arr), length(perturb));
zb_arr = zeros(length(X0_arr), length(perturb));

for (j,X0) in enumerate(X0_arr)
    if j == 1
        local X_arr = X0 .* perturb_SST;
    else
        local X_arr = X0 .* perturb;
    end
    
    # run simulation looping over X
    for (i,Xi) in enumerate(X_arr)
        # set up parameters
        local par = climatology();
        par.etype = enBal();

        par.SST0 = (j == 1) ? Xi : SST0 # (K)
        par.V = (j == 2) ? Xi : V0 # m/s
        par.D = (j == 3) ? Xi : D0 # (1/s)
        par.RHft = (j == 4) ? Xi : RHft0 # (-)
        EIS = (j == 5) ? Xi : EIS0 # (K)
        par.CO2 = (j == 6) ? Xi : CO20 # (ppm)

        par.sft0 = (par.SST0 + EIS)*Cp; # (K)
        par.Gamma_s = -6e-3*Cp + g; # (K/m)

        # run to equilibrium
        dt = 24.0;
        tmax = 50.0;
        u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));

        # extract end state
        uf = sol.u[end];
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
p1 = plot(X0_arr[j] .* perturb_SST, CF_arr[j,:]*100, lw=2, label="",
        color=:black,  yguidefontcolor=:black,
        ylim=(90,100), ylabel="Cloud fraction [%]",
        left_margin=15Plots.mm, right_margin=30Plots.mm,
        bottom_margin=15Plots.mm, top_margin=5Plots.mm)
plot!(twinx(), X0_arr[j] .* perturb_SST, LWP_arr[j,:]*1e3, lw=2, label="",
        color=:blue, yguidefontcolor=:blue, ytickfontcolor=:blue,
        ylim=(50,150), ylabel="In-cloud LWP [g/m\$^2\$]",
        xlabel=Xlabel_arr[j])
j=2
p2 = plot(X0_arr[j] .* perturb, CF_arr[j,:]*100, lw=2, label="",
        color=:black,  yguidefontcolor=:black,
        ylim=(90,100), ylabel="Cloud fraction [%]",
        left_margin=15Plots.mm, right_margin=30Plots.mm,
        bottom_margin=15Plots.mm, top_margin=5Plots.mm)
plot!(twinx(), X0_arr[j] .* perturb, LWP_arr[j,:]*1e3, lw=2, label="",
        color=:blue, yguidefontcolor=:blue, ytickfontcolor=:blue,
        ylim=(50,150), ylabel="In-cloud LWP [g/m\$^2\$]",
        xlabel=Xlabel_arr[j])
j=3
p3 = plot(1e6 .* X0_arr[j] .* perturb, CF_arr[j,:]*100, lw=2, label="",
        color=:black,  yguidefontcolor=:black,
        ylim=(90,100), ylabel="Cloud fraction [%]",
        left_margin=15Plots.mm, right_margin=30Plots.mm,
        bottom_margin=15Plots.mm, top_margin=5Plots.mm)
plot!(twinx(), 1e6 .* X0_arr[j] .* perturb, LWP_arr[j,:]*1e3, lw=2, label="",
        color=:blue, yguidefontcolor=:blue, ytickfontcolor=:blue,
        ylim=(50,150), ylabel="In-cloud LWP [g/m\$^2\$]",
        xlabel=Xlabel_arr[j])
j=4
p4 = plot(X0_arr[j] .* perturb, CF_arr[j,:]*100, lw=2, label="",
        color=:black,  yguidefontcolor=:black,
        ylim=(90,100), ylabel="Cloud fraction [%]",
        left_margin=15Plots.mm, right_margin=30Plots.mm,
        bottom_margin=15Plots.mm, top_margin=5Plots.mm)
plot!(twinx(), X0_arr[j] .* perturb, LWP_arr[j,:]*1e3, lw=2, label="",
        color=:blue, yguidefontcolor=:blue, ytickfontcolor=:blue,
        ylim=(50,150), ylabel="In-cloud LWP [g/m\$^2\$]",
        xlabel=Xlabel_arr[j])
j=5
p5 = plot(X0_arr[j] .* perturb, CF_arr[j,:]*100, lw=2, label="",
        color=:black,  yguidefontcolor=:black,
        ylim=(90,100), ylabel="Cloud fraction [%]",
        left_margin=15Plots.mm, right_margin=30Plots.mm,
        bottom_margin=15Plots.mm, top_margin=5Plots.mm)
plot!(twinx(), X0_arr[j] .* perturb, LWP_arr[j,:]*1e3, lw=2, label="",
        color=:blue, yguidefontcolor=:blue, ytickfontcolor=:blue,
        ylim=(50,150), ylabel="In-cloud LWP [g/m\$^2\$]",
        xlabel=Xlabel_arr[j])
j=6
p6 = plot(X0_arr[j] .* perturb, CF_arr[j,:]*100, lw=2, label="",
        color=:black,  yguidefontcolor=:black,
        ylim=(90,100), ylabel="Cloud fraction [%]",
        left_margin=15Plots.mm, right_margin=30Plots.mm,
        bottom_margin=15Plots.mm, top_margin=5Plots.mm)
plot!(twinx(), X0_arr[j] .* perturb, LWP_arr[j,:]*1e3, lw=2, label="",
        color=:blue, yguidefontcolor=:blue, ytickfontcolor=:blue,
        ylim=(50,150), ylabel="In-cloud LWP [g/m\$^2\$]",
        xlabel=Xlabel_arr[j])
plot(p1, p2, p3, p4, p5, p6, layout=(2,3), size=(1600,800), dpi=200);

path = "experiments/figures/linear_perturb/";
mkpath(path);
savefig(path*"linear_perturbations_around_Sc_state.png");
