# MLM behavior with fixed SST

Here we initialize the MLM from the steady-state condition with 400 ppm CO``_2`` and then calculate the new equilibrium state when we perturb the CO``_2`` to 800 ppm. In this example we still have fixed SSTs. 
```@example
using MixedLayerModel
using OrdinaryDiffEq
using SteadyStateDiffEq
using Plots
include("../../experiments/mlm_solve_funcs.jl")

# run simulation, 400 ppm (steady-state)
par = basic_params();
par.etype = enBal();
par.fttype = co2dep();
par.rtype = varRad();
par.stype = fixSST();
u0, sol_ss = run_mlm_ss(par);
uf = sol_ss.u;

# run simulation, 800 ppm
newCO2 = 800.0;
par.CO2 = newCO2;
u0, sol = run_mlm_from_init(uf, par, dt=3600.0*3.0, tspan=(0.0,3600.0*24.0*10.0));

# run simulation, 800 ppm (steady-state)
u0, sol_ss = run_mlm_ss_from_init(uf, par, dt=3600.0*3.0, tspan=3600.0*24.0*10.0);
# get output #hide
uf = sol_ss.u; #hide
du = zeros(4); #hide
mlm(du, uf, par, 0.0); #hide
zi,hM,qM,SST = uf; #hide
zb = calc_LCL(zi,hM,qM); #hide
LWP = calc_LWP(zi,hM,qM); #hide

# plot #hide
t = sol.t / 3600.0 / 24.0 #hide
zi = getindex.(sol.u,1) #hide
hM = getindex.(sol.u,2) * 1e-3 #hide
qtM = getindex.(sol.u,3) * 1e3 #hide
LWP = calc_LWP.(zi, hM*1e3, qtM*1e-3) * 1e3; #hide
plot(size=(600,500), layout=(4,1)) #hide
plot!(t, zi, marker="o-", legend=:topleft, subplot=1, label="zi(t) [m]") #hide
plot!(t, hM, marker="o-", legend=:topleft, subplot=2, label="hM(t) [kJ/kg]") #hide
plot!(t, qtM, marker="o-", legend=:topleft, subplot=3, label="qtM(t) [g/kg]") #hide
plot!(t, LWP, marker="o-", legend=:topleft, subplot=4, label="LWP(t) [g/m2]") #hide

uf = sol_ss.u; #hide
zi = uf[1]; #hide
hM = uf[2] * 1e-3; #hide
qtM = uf[3] * 1e3; #hide
LWP = calc_LWP(zi, hM*1e3, qtM*1e-3) * 1e3; #hide
hline!([zi], subplot=1, label=string(round(zi,digits=1))) #hide
hline!([hM], subplot=2, label=string(round(hM,digits=1))) #hide
hline!([qtM], subplot=3, label=string(round(qtM,digits=1))) #hide
hline!([LWP], subplot=4, label=string(round(LWP,digits=1))) #hide
title!(string(Int(newCO2))*" (ppm)", subplot=1) #hide
xaxis!("t (days)", subplot=3) #hide
```