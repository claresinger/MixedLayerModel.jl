exp_path = "20230419_co2_slow_ramp/";
path = "experiments/output/"*exp_path;

using MixedLayerModel
using Plots
include("mlm_solve_funcs.jl")

# create parameters
par = upCO2();
par.etype = enBal();
par.fttype = co2EIS();
par.rtype = varRad();
par.stype = fixSST();
par.ctype = fixCO2();

# adjust tunable parameters
par.Cd = 7.9e-4;
par.α_vent = 1.69e-3;
par.SW_b = 140;

# 400 ppm, fixed SST to get OHU
dt, tmax = 10.0, 100.0; # days
u0, sol = run_mlm(par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
uf = sol.u[end];
zb = calc_LCL(uf);
LWP = incloud_LWP(uf, zb);
OHU_400 = calc_OHU(uf,par,LWP,par.stype);
println(OHU_400)

# varying CO2, with OHU
dt, tmax = 10.0, 365.0*20; # days
par.stype = varSST();
par.ctype = varCO2();
par.rate = 1.1; # 10% increase per year
par.OHU = OHU_400;
u0, sol = run_mlm_from_init(uf, par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
time = sol.t / 3600 / 24 / 365 # years

# plot
plot(size=(800,400), layout=(1,1), dpi=200, 
    left_margin = 5Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm);
plot!(time, getindex.(sol.u,5) * 100, legend=false, xlabel="Time [years]", ylabel="Cloud fraction [%]")
plot!(twinx(), time, co2_rate.(sol.t, Ref(par), Ref(par.ctype)), color="red", legend=false, ylabel="CO2 [ppmv]")
mkpath(replace(path, "output"=>"figures"));
savefig(replace(path, "output"=>"figures")*"sol.png")
