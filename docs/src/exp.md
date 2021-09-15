# Running an experiment

First use the MLM to calculate the steady-state clouds given an atmosphere with 400 ppm CO``_2`` and fix the SST at 290 K.
```@example
using MixedLayerModel
using OrdinaryDiffEq
using SteadyStateDiffEq
using Plots
include("../../experiments/mlm_solve_funcs.jl")

# run simulation
par = basic_params();
par.etype = enBal();
par.stype = fixSST();
u0, sol = run_mlm(par);

# run simulation (steady-state)
u0, sol_ss = run_mlm_ss(par);
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
plot!(t, zi, marker="o-", legend=:topright, subplot=1, label="zi(t) [m]") #hide
plot!(t, hM, marker="o-", legend=:topright, subplot=2, label="hM(t) [kJ/kg]") #hide
plot!(t, qtM, marker="o-", legend=:topright, subplot=3, label="qtM(t) [g/kg]") #hide
plot!(t, LWP, marker="o-", legend=:topright, subplot=4, label="LWP(t) [g/m2]") #hide

uf = sol_ss.u; #hide
zi = uf[1]; #hide
hM = uf[2] * 1e-3; #hide
qtM = uf[3] * 1e3; #hide
LWP = calc_LWP(zi, hM*1e3, qtM*1e-3) * 1e3; #hide
hline!([zi], subplot=1, label=string(round(zi,digits=1))) #hide
hline!([hM], subplot=2, label=string(round(hM,digits=1))) #hide
hline!([qtM], subplot=3, label=string(round(qtM,digits=1))) #hide
hline!([LWP], subplot=4, label=string(round(LWP,digits=1))) #hide
title!("400 (ppm)", subplot=1) #hide
xaxis!("t (days)", subplot=3) #hide
```