# Running an experiment

First use the MLM to calculate the steady-state clouds given an atmosphere with 400 ppm CO``_2`` and fix the SST at 290 K.
```@example
using MixedLayerModel
using OrdinaryDiffEq
using Plots
include("../../experiments/mlm_solve_funcs.jl")

# run simulation
par = upCO2();
par.etype = enBal();
par.fttype = co2dep();
par.rtype = varRad();
par.stype = fixSST();
dt = 2.0;
tmax = 40.0;
ENV["GKSwstype"]="nul"
u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));

# plot #hide
t = sol.t / 3600.0 / 24.0 #hide
zi = getindex.(sol.u,1) #hide
zb = calc_LCL.(sol.u) #hide
hM = getindex.(sol.u,2) * 1e-3 #hide
qtM = getindex.(sol.u,3) * 1e3 #hide
SST = getindex.(sol.u,4) #hide
CF = getindex.(sol.u,5) #hide
LWP = incloud_LWP.(sol.u, zb) #hide
plot(size=(600,500), layout=(5,1)) #hide
plot!(t, zi, marker="o-", label="", subplot=1, ylabel="zi, zb [m]") #hide
plot!(t, zb, marker="o-", subplot=1, label="") #hide
plot!(t, hM, marker="o-", label="", subplot=2, ylabel="hM [kJ/kg]") #hide
plot!(t, qtM, marker="o-", label="", subplot=3, ylabel="qtM [g/kg]") #hide
plot!(t, SST, marker="o-", label="", subplot=4, ylabel="SST [K]") #hide
plot!(t, CF*1e2, marker="o-", label="", subplot=5, ylabel="CF [%]") #hide
```