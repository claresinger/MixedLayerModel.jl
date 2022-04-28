# import necessary packages
push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using MixedLayerModel
using Plots

include("mlm_solve_funcs.jl")

SST_arr = [288, 289, 290, 291, 292];
D_arr = [5.5, 5.75, 6, 6.25, 6.5]*1e-6;
RHft_arr = [0.18, 0.19, 0.2, 0.21, 0.22];
V_arr = [5,6,7,8,9];

X_arr = SST_arr;
X_label = "SST (K)"; 
# X_arr = D_arr;
# X_label = "D (1/s)";  
# X_arr = V_arr;
# X_label = "V (m/s)"; 
# X_arr = RHft_arr;
# X_label = "RHft";

CF_arr = zeros(length(X_arr));
LWP_arr = zeros(length(X_arr));
# run simulation looping over X
for (i,Xi) in enumerate(X_arr)
    # set up parameters
    par = climatology();
    par.etype = bflux(); #use buoyancy flux entrainment closure

    par.SST0 = Xi # (K), SST_arr[3]
    par.D = D_arr[3] # (1/s), D_arr[3]
    par.V = V_arr[3] # m/s, V_arr[3]
    par.RHft = RHft_arr[3] # (-), RHft_arr[3]

    par.Gamma_q = 0.0; # (kg/kg/m)
    par.sft0 = 305; # (K)
    par.Gamma_s = 0.0; # (K/m)

    # run to equilibrium
    dt = 24.0;
    tmax = 40.0;
    u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));

    # extract end state
    uf = sol.u[end];
    zi, hM, qM, SST, CFi = uf;
    zb = calc_LCL(uf);
    LWPi = incloud_LWP(uf, zb);
    CF_arr[i] = CFi;
    LWP_arr[i] = LWPi;
end

# print results
println(X_arr)
println(CF_arr)
println(LWP_arr*1e3)

# plot and save results
ENV["GKSwstype"]="nul"
plot(size=(800,500), dpi=200, left_margin = 5Plots.mm);
plot!(X_arr, LWP_arr*1e3, marker="o-", xlabel=X_label, ylabel="LWP (g/m2)", legend=false);
path = "experiments/figures/linear_perturb/";
mkpath(path);
savefig(path*"lin_perturb_"*split(X_label, " (")[1]*".png");
