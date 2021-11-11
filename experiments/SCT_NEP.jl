using MixedLayerModel
using FileIO
using Plots

include("mlm_solve_funcs.jl")

# use command line argument to set co2
i = parse(Int,ARGS[1]);
println(i);

# load initial condition from file
path = "experiments/output/SCT_NEP/";
restartfile = path*"mean.jld2"
output = load(restartfile);
u0 = output["uf"];

# set OHU, increase CO2, let SST evolve and check cloud changes
par = climatology();
par.etype = enBal();

# climatology at 4 sites in NEP
SSTs = [290.2, 293.5, 297.2, 299.5];
Ds = [6e-6, 6e-6, 3e-6, 3e-6]; # ??????????????
Vs = [2.7, 5.2, 8.0, 9.0];
RHfts = [0.25, 0.21, 0.24, 0.35];
sft0s = [308.6, 309.6, 311.4, 312.9];

par.SST0 = SSTs[i]; # (K)
par.D = Ds[i]; # (1/s)
par.V = Vs[i] # m/s 

par.RHft = RHfts[i]; # (-) 
par.Gamma_q = 0.0; # (kg/kg/m)
par.sft0 = sft0s[i]; # (K) 
par.Gamma_s = 0.0; # (K/m)

dt, tmax = 4.0, 30;

# plot time series
ENV["GKSwstype"]="nul"
u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
#u0, sol = run_mlm_from_init(u0, par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
t = sol.t / 3600.0 / 24.0;
zi = getindex.(sol.u,1);
hM = getindex.(sol.u,2) * 1e-3;
qtM = getindex.(sol.u,3) * 1e3;
sst = getindex.(sol.u,4);
cf = getindex.(sol.u,5);
S = zeros(length(t));
LHF = zeros(length(t));
zb = zeros(length(t));
ΔR = zeros(length(t));
LWP = zeros(length(t));
for (i,si) in enumerate(S)
    zb[i] = calc_LCL(sol.u[i]);
    S[i] = calc_S(sol.u[i], par);
    LHF[i] = calc_LHF(sol.u[i], par);
    ΔR[i] = calc_cloudtop_RAD(sol.u[i], par, par.rtype);
    LWP[i] = incloud_LWP(sol.u[i]);
end 
plot(size=(1200,800), layout=(5,2), dpi=200, left_margin = 5Plots.mm);
plot!(t, zi, marker="o-", legend=false, subplot=1, ylabel="zi, zb [m]");
plot!(t, zb, marker="o-", legend=false, subplot=1);
plot!(t, hM, marker="o-", legend=false, subplot=2, ylabel="hM [kJ/kg]"); 
plot!(t, qtM, marker="o-", legend=false, subplot=3, ylabel="qtM [g/kg]");
plot!(t, sst, marker="o-", legend=false, subplot=4, ylabel="SST [K]");
plot!(t, cf * 1e2, marker="o-", legend=false, subplot=5, ylabel="CF [%]");
plot!(t, LWP .* cf * 1e3, marker="o-", legend=false, subplot=6, ylabel="LWP [g/m2]");
plot!(t, LHF, marker="o-", legend=false, subplot=7, ylabel="LHF [W/m2]");
plot!(t, ΔR, marker="o-", legend=false, subplot=8, ylabel="ΔR [W/m2]");
plot!(t, (zi .- zb) ./ zi, marker="o-", legend=false, subplot=9, ylabel="zc/zi [-]", xlabel="time [days]");
plot!(t, S, marker="o-", legend=false, subplot=10, ylabel="S [-]", xlabel="time [days]");
mkpath(replace(path, "output"=>"figures"));
savefig(replace(path, "output"=>"figures")*"sol"*string(i)*"_t.png");

## save steady-state solution
uf = sol.u[end];
du = zeros(5);
mlm(du, uf, par, 0.0);
zi,hM,qM,SST = uf;
zb = calc_LCL(uf);
println(uf);
println(du);

output = Dict("p" => par, "u0" => u0, "uf" => uf, "du/u" => du./uf, 
"we" => we(uf,par,par.etype), "zb" => zb, "zc" => zi-zb,
"RHsurf" => RH(0.0, hM, qM), "LHF" => calc_LHF(uf,par), "SHF" => calc_SHF(uf,par),
"ΔR" => calc_cloudtop_RAD(uf,par,par.rtype), "OHU" => calc_OHU(uf,par,par.stype))

save(path*"sol_"*string(i)*".jld2", output)