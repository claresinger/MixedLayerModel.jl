# MLM behavior with slab ocean

Instead if we want to allow SSTs to change interactively, then we can set the ocean heat uptake (OHU) to the value defined by the 400 ppm simulation. Also be sure to switch to `varSST()` mode.
```@example
using MixedLayerModel
using OrdinaryDiffEq
using SteadyStateDiffEq
using Plots
include("../../experiments/mlm_solve_funcs.jl")

# run simulation, 400 ppm (steady-state)
par = upCO2();
par.etype = enBal();
par.fttype = co2dep();
par.rtype = varRad();
par.stype = fixSST();
u0, sol_ss = run_mlm_ss(par);
uf = sol_ss.u;
OHU = calc_OHU(uf,par,par.stype);

# run simulation, 800 ppm
newCO2 = 800.0;
par.CO2 = newCO2;
par.OHU = OHU;
par.stype = varSST();
u0, sol = run_mlm_from_init(uf, par, dt=3600.0*2.0, tspan=(0.0,3600.0*24.0*15.0));

# run simulation, 800 ppm (steady-state)
u0, sol_ss = run_mlm_ss_from_init(uf, par, dt=3600.0*2.0, tspan=3600.0*24.0*15.0);
# get output #hide
uf = sol_ss.u; #hide
du = zeros(5); #hide
mlm(du, uf, par, 0.0); #hide
zi,sM,qM,SST,CF = uf; #hide
zb = calc_LCL(uf); #hide
LWP = incloud_LWP(uf); #hide

# plot #hide
t = sol.t / 3600.0 / 24.0 #hide
zi = getindex.(sol.u,1) #hide
sM = getindex.(sol.u,2) * 1e-3 #hide
qtM = getindex.(sol.u,3) * 1e3 #hide
LWP = incloud_LWP.(sol.u) * 1e3; #hide
plot(size=(600,500), layout=(4,1)) #hide
plot!(t, zi, marker="o-", legend=:topleft, subplot=1, label="zi(t) [m]") #hide
plot!(t, sM, marker="o-", legend=:topleft, subplot=2, label="sM(t) [kJ/kg]") #hide
plot!(t, qtM, marker="o-", legend=:topleft, subplot=3, label="qtM(t) [g/kg]") #hide
plot!(t, LWP, marker="o-", legend=:topleft, subplot=4, label="LWP(t) [g/m2]") #hide

uf = sol_ss.u; #hide
zi = uf[1]; #hide
sM = uf[2] * 1e-3; #hide
qtM = uf[3] * 1e3; #hide
LWP = incloud_LWP(uf) * 1e3; #hide
hline!([zi], subplot=1, label=string(round(zi,digits=1))) #hide
hline!([sM], subplot=2, label=string(round(sM,digits=1))) #hide
hline!([qtM], subplot=3, label=string(round(qtM,digits=1))) #hide
hline!([LWP], subplot=4, label=string(round(LWP,digits=1))) #hide
title!(string(Int(newCO2))*" (ppm)", subplot=1) #hide
xaxis!("t (days)", subplot=3) #hide
```