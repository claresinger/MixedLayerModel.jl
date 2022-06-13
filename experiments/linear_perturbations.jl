# import necessary packages
push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using MixedLayerModel
using MixedLayerModel: T0
using Plots

include("mlm_solve_funcs.jl")

SST0 = 25.0;
D0 = 4e-6;
V0 = 8.0;
RHft0 = 0.2;
EIS0 = 10.0;
CO20 = 400.0;

perturb = collect(80 : 5 : 120) ./ 100;

# X_arr = SST0 .* perturb;
# X_label = "SST (K)"; 
# X_arr = D0 .* perturb;
# X_label = "D (1/s)";  
# X_arr = V0 .* perturb;
# X_label = "V (m/s)"; 
# X_arr = RHft0 .* perturb;
# X_label = "RHft";
# X_arr = CO20 .* perturb;
# X_label = "CO2 (ppm)";

X0_arr = [SST0, V0, D0, RHft0, EIS0, CO20];
Xlabel_arr = ["SST (K)", "V (m/s)", "D (1/s)", "RHft", "EIS (K)", "CO2 (ppm)"];

CF_arr = zeros(length(X0_arr), length(perturb));
LWP_arr = zeros(length(X0_arr), length(perturb));
zi_arr = zeros(length(X0_arr), length(perturb));
zb_arr = zeros(length(X0_arr), length(perturb));

for (j,X0) in enumerate(X0_arr)
    local X_arr = X0 .* perturb;
    
    # run simulation looping over X
    for (i,Xi) in enumerate(X_arr)
        # set up parameters
        local par = climatology();
        par.etype = enBal();

        par.SST0 = (j == 1) ? Xi + T0 : SST0 + T0 # (K)
        par.V = (j == 2) ? Xi : V0 # m/s
        par.D = (j == 3) ? Xi : D0 # (1/s)
        par.RHft = (j == 4) ? Xi : RHft0 # (-)
        EIS = (j == 5) ? Xi : EIS0 # (K)
        par.CO2 = (j == 6) ? Xi : CO20 # (ppm)

        par.Gamma_q = -1e-6; # (kg/kg/m)
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
println(X_arr)
println(CF_arr)
println(LWP_arr*1e3)

ENV["GKSwstype"]="nul"

# plot and save results
j=1
p1 = plot(X0_arr[j] .* perturb .+ T0, CF_arr[j,:]*100, lw=2, label="",
        color=:black,  yguidefontcolor=:black,
        ylim=(80,100), ylabel="Cloud fraction [%]",
        left_margin=10Plots.mm, right_margin=25Plots.mm,
        bottom_margin=10Plots.mm)
plot!(twinx(), X0_arr[j] .* perturb .+ T0, LWP_arr[j,:]*1e3, lw=2, label="",
        color=:blue, yguidefontcolor=:blue,
        ylim=(100,250), ylabel="In-cloud LWP [g/m2]",
        xlabel=Xlabel_arr[j])
j=2
p2 = plot(X0_arr[j] .* perturb, CF_arr[j,:]*100, lw=2, label="",
        color=:black,  yguidefontcolor=:black,
        ylim=(80,100), ylabel="Cloud fraction [%]",
        left_margin=10Plots.mm, right_margin=25Plots.mm,
        bottom_margin=10Plots.mm)
plot!(twinx(), X0_arr[j] .* perturb, LWP_arr[j,:]*1e3, lw=2, label="",
        color=:blue, yguidefontcolor=:blue,
        ylim=(100,250), ylabel="In-cloud LWP [g/m2]",
        xlabel=Xlabel_arr[j])
j=3
p3 = plot(X0_arr[j] .* perturb, CF_arr[j,:]*100, lw=2, label="",
        color=:black,  yguidefontcolor=:black,
        ylim=(80,100), ylabel="Cloud fraction [%]",
        left_margin=10Plots.mm, right_margin=25Plots.mm,
        bottom_margin=10Plots.mm)
plot!(twinx(), X0_arr[j] .* perturb, LWP_arr[j,:]*1e3, lw=2, label="",
        color=:blue, yguidefontcolor=:blue,
        ylim=(100,250), ylabel="In-cloud LWP [g/m2]",
        xlabel=Xlabel_arr[j])
j=4
p4 = plot(X0_arr[j] .* perturb, CF_arr[j,:]*100, lw=2, label="",
        color=:black,  yguidefontcolor=:black,
        ylim=(80,100), ylabel="Cloud fraction [%]",
        left_margin=10Plots.mm, right_margin=25Plots.mm,
        bottom_margin=10Plots.mm)
plot!(twinx(), X0_arr[j] .* perturb, LWP_arr[j,:]*1e3, lw=2, label="",
        color=:blue, yguidefontcolor=:blue,
        ylim=(100,250), ylabel="In-cloud LWP [g/m2]",
        xlabel=Xlabel_arr[j])
j=5
p5 = plot(X0_arr[j] .* perturb, CF_arr[j,:]*100, lw=2, label="",
        color=:black,  yguidefontcolor=:black,
        ylim=(80,100), ylabel="Cloud fraction [%]",
        left_margin=10Plots.mm, right_margin=25Plots.mm,
        bottom_margin=10Plots.mm)
plot!(twinx(), X0_arr[j] .* perturb, LWP_arr[j,:]*1e3, lw=2, label="",
        color=:blue, yguidefontcolor=:blue,
        ylim=(100,250), ylabel="In-cloud LWP [g/m2]",
        xlabel=Xlabel_arr[j])
j=6
p6 = plot(X0_arr[j] .* perturb, CF_arr[j,:]*100, lw=2, label="",
        color=:black,  yguidefontcolor=:black,
        ylim=(80,100), ylabel="Cloud fraction [%]",
        left_margin=10Plots.mm, right_margin=25Plots.mm,
        bottom_margin=10Plots.mm)
plot!(twinx(), X0_arr[j] .* perturb, LWP_arr[j,:]*1e3, lw=2, label="",
        color=:blue, yguidefontcolor=:blue,
        ylim=(100,250), ylabel="In-cloud LWP [g/m2]",
        xlabel=Xlabel_arr[j])
plot(p1, p2, p3, p4, p5, p6, layout=(2,3), size=(1600,800), dpi=200);

path = "experiments/figures/linear_perturb/";
mkpath(path);
savefig(path*"all.png");
