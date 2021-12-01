push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using MixedLayerModel
using FileIO
using Plots

include("mlm_solve_funcs.jl")

# use command line argument to set co2
# newCO2 = parse(Float64,ARGS[1]);
newCO2 = 600.0;
println(newCO2);

# load initial condition from file
path = "experiments/output/zbout/";
restarttry1 = path*"co2_upstep_"*string(Int(newCO2-100))*".jld2";
restarttry2 = path*"co2_upstep_"*string(Int(newCO2-200))*".jld2";
restarttry3 = path*"co2_upstep_"*string(Int(newCO2-400))*".jld2";
if isfile(restarttry1)
    output = load(restarttry1);
elseif isfile(restarttry2)
    output = load(restarttry2);
elseif isfile(restarttry3)
    output = load(restarttry3);
else
    output = load(path*"co2_400.jld2");
end
u0 = output["uf"];
OHU = output["OHU"];
println("restarting from CO2 = "*string(output["p"].CO2));

# set OHU, increase CO2, let SST evolve and check cloud changes
par = upCO2();
par.Hw = 0.1;
par.OHU = OHU;
par.CO2 = newCO2;
par.etype = enBal();
par.fttype = co2dep();
par.rtype = varRad();
par.stype = varSST();
dt, tmax = 2.0, 30;

# u0, sol = run_mlm_ss_from_init(u0, par, dt=3600.0*dt, tspan=3600.0*24.0*tmax);
# code = sol.retcode;
# println(code);

# uf = sol.u;
# du = zeros(5);
# mlm(du, uf, par, 0.0);
# zi,hM,qM,SST = uf;
# zb = calc_LCL(uf);
# println(uf);
# println(du);

# output = Dict("code" => code, "p" => par, "u0" => u0, "uf" => uf, "du/u" => du./uf, 
# "we" => we(uf,par,par.etype), "zb" => zb, "zc" => zi-zb,
# "RHsurf" => RH(0.0, hM, qM), "LHF" => calc_LHF(uf,par), "SHF" => calc_SHF(uf,par),
# "ΔR" => calc_cloudtop_RAD(uf,par,par.rtype), "OHU" => calc_OHU(uf,par,par.stype))

# save(path*"co2_upstep_"*string(Int(newCO2))*".jld2", output)

# plot time series
ENV["GKSwstype"]="nul"
u0, sol = run_mlm_from_init(u0, par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));

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
    S[i] = calc_S(sol.u[i], par, zb[i]);
    LHF[i] = calc_LHF(sol.u[i], par);
    ΔR[i] = calc_cloudtop_RAD(sol.u[i], par, zb[i], par.rtype);
    LWP[i] = incloud_LWP(sol.u[i], zb[i]);
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
savefig(replace(path, "output"=>"figures")*"sol"*string(Int(newCO2))*"_t.png");

## save steady-state solution
uf = sol.u[end];
du = zeros(5);
mlm(du, uf, par, 0.0);
zi,hM,qM,SST = uf;
zb = calc_LCL(uf);
println(uf);
println(du);

output = Dict("p" => par, "u0" => u0, "uf" => uf, "du/u" => du./uf, 
"we" => we(uf,par,zb,par.etype), "zb" => zb, "zc" => zi-zb,
"RHsurf" => RH(0.0, hM, qM), "LHF" => calc_LHF(uf,par), "SHF" => calc_SHF(uf,par),
"ΔR" => calc_cloudtop_RAD(uf,par,zb,par.rtype), "OHU" => calc_OHU(uf,par,zb,par.stype))

save(path*"co2_upstep_"*string(Int(newCO2))*".jld2", output)