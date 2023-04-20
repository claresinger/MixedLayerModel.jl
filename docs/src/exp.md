# Running an experiment

## Generate initial equilibrium condition
First use the MLM to calculate the steady-state clouds given an atmosphere with 400 ppm CO``_2`` and fix the SST at 290 K.
```@example
using MixedLayerModel
using OrdinaryDiffEq
using Plots
include("../../experiments/mlm_solve_funcs.jl")

# run simulation
par = upCO2();
par.etype = enBal();
par.fttype = co2EIS();
par.rtype = varRad();
par.stype = fixSST();
dt = 12.0;
tmax = 50.0;
ENV["GKSwstype"]="nul"
u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
uf = sol.u[end];

# plot #hide
t = sol.t / 3600.0 / 24.0 #hide
zi = getindex.(sol.u,1) #hide
zb = calc_LCL.(sol.u) #hide
sM = getindex.(sol.u,2) * 1e-3 #hide
qtM = getindex.(sol.u,3) * 1e3 #hide
SST = getindex.(sol.u,4) #hide
CF = getindex.(sol.u,5) #hide
LWP = incloud_LWP.(sol.u, zb) #hide
plot(size=(600,500), layout=(5,1)) #hide
plot!(t, zi, marker="o-", label="", subplot=1, ylabel="zi, zb [m]") #hide
plot!(t, zb, marker="o-", subplot=1, label="") #hide
plot!(t, sM, marker="o-", label="", subplot=2, ylabel="sM [kJ/kg]") #hide
plot!(t, qtM, marker="o-", label="", subplot=3, ylabel="qtM [g/kg]") #hide
plot!(t, SST, marker="o-", label="", subplot=4, ylabel="SST [K]") #hide
plot!(t, CF*1e2, marker="o-", label="", subplot=5, ylabel="CF [%]") #hide
xaxis!("t (days)", subplot=5) #hide
```

## Add CO``_2`` with fixed SSTs
Here we initialize the MLM from the steady-state condition with 400 ppm CO``_2`` and then calculate the new equilibrium state when we perturb the CO``_2`` to 800 ppm. In this example we still have fixed SSTs. 

```@example
using MixedLayerModel #hide
using OrdinaryDiffEq #hide
using Plots #hide
include("../../experiments/mlm_solve_funcs.jl") #hide
par = upCO2(); #hide
par.etype = enBal(); #hide
par.fttype = co2EIS(); #hide
par.rtype = varRad(); #hide
par.stype = fixSST(); #hide
dt = 12.0; #hide
tmax = 50.0; #hide
ENV["GKSwstype"]="nul" #hide
u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax)); #hide
uf = sol.u[end]; #hide

# run simulation, 800 ppm
newCO2 = 800.0;
par.CO2 = newCO2;
u0, sol = run_mlm_from_init(uf, par, dt=3600.0*12.0, tspan=(0.0,3600.0*24.0*50.0));

# plot #hide
t = sol.t / 3600.0 / 24.0 #hide
zi = getindex.(sol.u,1) #hide
sM = getindex.(sol.u,2) * 1e-3 #hide
qtM = getindex.(sol.u,3) * 1e3 #hide
SST = getindex.(sol.u,4) #hide
CF = getindex.(sol.u,5) #hide
zb = calc_LCL.(sol.u) #hide
LWP = incloud_LWP.(sol.u, zb) * 1e3 #hide
plot(size=(600,900), layout=(6,1)) #hide
plot!(t, zi, marker="o-", legend=:topleft, subplot=1, label="zi(t) [m]") #hide
plot!(t, zb, marker="o-", legend=:topleft, subplot=1, label="zb(t) [m]") #hide
plot!(t, sM, marker="o-", legend=:topleft, subplot=2, label="sM(t) [kJ/kg]") #hide
plot!(t, qtM, marker="o-", legend=:topleft, subplot=3, label="qtM(t) [g/kg]") #hide
plot!(t, LWP, marker="o-", legend=:topleft, subplot=4, label="LWP(t) [g/m2]") #hide
plot!(t, SST, marker="o-", legend=:topleft, subplot=5, label="SST(t) [K]") #hide
plot!(t, CF*100, marker="o-", legend=:topleft, subplot=6, label="CF(t) [%]") #hide

title!(string(Int(newCO2))*" (ppm)", subplot=1) #hide
xaxis!("t (days)", subplot=6) #hide
```

## Add CO``_2`` and let SSTs dynamically adjust
Instead if we want to allow SSTs to change interactively, then we can set the ocean heat uptake (OHU) to the value defined by the 400 ppm simulation. Also be sure to switch to `varSST()` mode.

```@example
using MixedLayerModel #hide
using OrdinaryDiffEq #hide
using Plots #hide
include("../../experiments/mlm_solve_funcs.jl") #hide
par = upCO2(); #hide
par.etype = enBal(); #hide
par.fttype = co2EIS(); #hide
par.rtype = varRad(); #hide
par.stype = fixSST(); #hide
dt = 12.0; #hide
tmax = 50.0; #hide
ENV["GKSwstype"]="nul" #hide
u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax)); #hide
uf = sol.u[end]; #hide

zb = calc_LCL(uf);
LWP = incloud_LWP(uf, zb);
OHU = calc_OHU(uf, par, LWP, par.stype);

# run simulation, 800 ppm
newCO2 = 800.0;
par.CO2 = newCO2;
par.OHU = OHU;
par.stype = varSST();
u0, sol = run_mlm_from_init(uf, par, dt=3600.0*12.0, tspan=(0.0,3600.0*24.0*50.0));

# plot #hide
t = sol.t / 3600.0 / 24.0 #hide
zi = getindex.(sol.u,1) #hide
sM = getindex.(sol.u,2) * 1e-3 #hide
qtM = getindex.(sol.u,3) * 1e3 #hide
SST = getindex.(sol.u,4) #hide
CF = getindex.(sol.u,5) #hide
zb = calc_LCL.(sol.u) #hide
LWP = incloud_LWP.(sol.u, zb) * 1e3 #hide
plot(size=(600,900), layout=(6,1)) #hide
plot!(t, zi, marker="o-", legend=:topleft, subplot=1, label="zi(t) [m]") #hide
plot!(t, zb, marker="o-", legend=:topleft, subplot=1, label="zb(t) [m]") #hide
plot!(t, sM, marker="o-", legend=:topleft, subplot=2, label="sM(t) [kJ/kg]") #hide
plot!(t, qtM, marker="o-", legend=:topleft, subplot=3, label="qtM(t) [g/kg]") #hide
plot!(t, LWP, marker="o-", legend=:topleft, subplot=4, label="LWP(t) [g/m2]") #hide
plot!(t, SST, marker="o-", legend=:topleft, subplot=5, label="SST(t) [K]") #hide
plot!(t, CF*100, marker="o-", legend=:topleft, subplot=6, label="CF(t) [%]") #hide

title!(string(Int(newCO2))*" (ppm)", subplot=1) #hide
xaxis!("t (days)", subplot=6) #hide
```
