# MLM behavior with fixed SST

Here we initialize the MLM from the steady-state condition with 400 ppm CO``_2`` and then calculate the new equilibrium state when we perturb the CO``_2`` to 800 ppm. In this example we still have fixed SSTs. 
```@example
using MixedLayerModel
using OrdinaryDiffEq
using Plots
include("../../experiments/mlm_solve_funcs.jl")

# run simulation, 400 ppm (steady-state)
par = upCO2();
par.etype = enBal();
par.fttype = co2EIS();
par.rtype = varRad();
par.stype = fixSST();
u0, sol = run_mlm(par, dt=3600.0*3.0, tspan=(0.0,3600.0*24.0*10.0));
uf = sol.u[end];

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
xaxis!("t (days)", subplot=3) #hide
```
