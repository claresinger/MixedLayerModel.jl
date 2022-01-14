push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using MixedLayerModel
using MixedLayerModel: μ, L0
using FileIO
using Plots

include("mlm_solve_funcs.jl")

# use command line argument to set co2
newCO2 = parse(Float64,ARGS[1]);
# newCO2 = 200.0;
println(newCO2);

# load initial condition from file
path = "experiments/output/modCF/";
restarttry1 = path*"co2_downstep_"*string(Int(newCO2+100))*".jld2";
restarttry2 = path*"co2_downstep_"*string(Int(newCO2+200))*".jld2";
restarttry3 = path*"co2_downstep_"*string(Int(newCO2+400))*".jld2";
if isfile(restarttry1)
    output = load(restarttry1);
elseif isfile(restarttry2)
    output = load(restarttry2);
elseif isfile(restarttry3)
    output = load(restarttry3);
else
    output = load(path*"co2_upstep_1600.jld2");
end
u0 = output["uf"];
OHU = output["OHU"];
println("restarting from CO2 = "*string(output["p"].CO2));

# get toa net rad @ 400 ppm
output = load(path*"co2_400.jld2");
u400 = output["uf"];
zb = output["zb"];
R_s_400 = toa_net_rad(u400, zb);

# set OHU, increase CO2, let SST evolve and check cloud changes
par = upCO2();
par.Hw = 0.1;
par.OHU = OHU;
par.R_s_400 = R_s_400;
par.CO2 = newCO2;
par.etype = enBal();
par.fttype = co2dep();
par.rtype = varRad();
par.stype = varSST();
dt, tmax = 24.0, 60.0;

println(par.OHU, "\t", par.R_s_400);

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

# plot time series
ENV["GKSwstype"]="nul"
u0, sol = run_mlm_from_init(u0, par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));

t = sol.t / 3600.0 / 24.0;
zi = getindex.(sol.u,1);
hM = getindex.(sol.u,2);
qtM = getindex.(sol.u,3);
sst = getindex.(sol.u,4);
cf = getindex.(sol.u,5);
S = zeros(length(t));
LHF = zeros(length(t));
zb = zeros(length(t));
ΔR = zeros(length(t));
LWP = zeros(length(t));
Δsvl = zeros(length(t));
ent = zeros(length(t));
trop_SST = zeros(length(t));
for (i,si) in enumerate(S)
    zb[i] = calc_LCL(sol.u[i]);
    LWP[i] = incloud_LWP(sol.u[i], zb[i]);
    S[i] = calc_S(sol.u[i], par, zb[i], LWP[i]);
    LHF[i] = calc_LHF(sol.u[i], par);
    ΔR[i] = calc_cloudtop_RAD(sol.u[i], par, LWP[i], par.rtype);
    Δsvl[i] = Δs(sol.u[i], par, LWP[i]);
    ent[i] = we(sol.u[i], par, zb[i], LWP[i], par.etype);
    trop_SST[i] = trop_sst(sol.u[i], par, LWP[i]);
end 
plot(size=(1200,800), layout=(6,2), dpi=200, left_margin = 5Plots.mm);
plot!(t, zi, marker="o-", legend=false, subplot=1, ylabel="zi, zb [m]");
plot!(t, zb, marker="o-", legend=false, subplot=1);
plot!(t, hM * 1e-3, marker="o-", legend=false, subplot=2, ylabel="hM [kJ/kg]"); 
plot!(t, qtM * 1e3, marker="o-", legend=false, subplot=3, ylabel="qtM [g/kg]");
plot!(t, sst, marker="o-", legend=false, subplot=4, ylabel="SST [K]");
#plot!(t, trop_SST, marker="o-", legend=false, subplot=4);
plot!(t, cf * 1e2, marker="o-", legend=false, subplot=5, ylabel="CF [%]");
plot!(t, LWP .* cf * 1e3, marker="o-", legend=false, subplot=6, ylabel="LWP [g/m2]");
plot!(t, Δsvl * 1e-3, marker="o-", legend=false, subplot=7, ylabel="Δs (kJ/kg)");
plot!(t, ent*1e3, marker="o-", legend=false, subplot=8, ylabel="we (mm/s)")
plot!(t, LHF, marker="o-", legend=false, subplot=9, ylabel="LHF [W/m2]");
plot!(t, ΔR, marker="o-", legend=false, subplot=10, ylabel="ΔR [W/m2]");
plot!(t, (zi .- zb) ./ zi, marker="o-", legend=false, subplot=11, ylabel="zc/zi [-]", xlabel="time [days]");
plot!(t, S, marker="o-", legend=false, subplot=12, ylabel="S [-]", xlabel="time [days]");
mkpath(replace(path, "output"=>"figures"));
savefig(replace(path, "output"=>"figures")*"down"*string(Int(newCO2))*"_t.png");

# ### AGU plots
# Plots.scalefontsizes(2)
# plot(size=(1000,500), layout=(2,2), dpi=200, left_margin = 5Plots.mm, bottom_margin=5Plots.mm);
# plot!(t, ΔR, marker="o-", legend=false, subplot=1, ylabel="ΔR [W/m\$^2\$]");
# plot!(t, S, marker="o-", legend=false, subplot=2, ylabel="Stability, \$S\$",
#         yscale=:log10, yticks=([0.2,0.5,2,5], ["0.2","0.5","2","5"]), ylim=[0.2,10]);
# plot!(t, cf * 1e2, marker="o-", legend=false, subplot=3, ylabel="CF [%]", xlabel="Time [days]");
# plot!(t, sst, marker="o-", legend=false, subplot=4, ylabel="SST [K]", xlabel="Time [days]");
# mkpath(replace(path, "output"=>"figures"));
# savefig(replace(path, "output"=>"figures")*"AGU-down"*string(Int(newCO2))*"_t.png");
# Plots.scalefontsizes(1/2)
# ### 

## save steady-state solution
uf = sol.u[end];
du = zeros(5);
mlm(du, uf, par, 0.0);
zi,hM,qM,SST,CF = uf;
zb = calc_LCL(uf);
LWP = incloud_LWP(uf, zb);
RH = min(qM / q_sat(0.0, temp(0.0, hM, qM)), 1.0);
println(uf);
println(du);
println("cloud base: ",zb)
println("LWP: ", LWP);
println("tropical sst: ", trop_sst(uf, par, LWP));
println("ft qt: ", qjump(uf, par, zb, par.fttype) + qM);

output = Dict("p" => par, "u0" => u0, "uf" => uf, "du/u" => du./uf, 
"we" => we(uf,par,zb,LWP,par.etype), "zb" => zb, "zc" => zi-zb,
"RHsurf" => RH, "LHF" => calc_LHF(uf,par), "SHF" => calc_SHF(uf,par),
"ΔR" => calc_cloudtop_RAD(uf,par,LWP,par.rtype), 
"OHU" => calc_OHU(uf,par,LWP,par.stype))

save(path*"co2_downstep_"*string(Int(newCO2))*".jld2", output)