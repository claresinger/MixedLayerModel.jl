exp_path = "20221211_calibparam_fixWV/";
path = "experiments/output/"*exp_path;

using MixedLayerModel
using JLD2
include("mlm_solve_funcs.jl")
include("plot_transient_solution.jl")

# create parameters
par = upCO2();
par.etype = enBal();
par.fttype = co2EIS();
par.rtype = varRad();
par.stype = fixSST();
par.wvtype = wvRADOFF();
dt, tmax = 10.0, 100.0; # days

# adjust tunable parameters
par.Cd = 7.9e-4;
par.α_vent = 1.69e-3;
par.SW_b = 140;

# # 400 ppm
# u0, sol = run_mlm(par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
# uf = sol.u[end];
# zb = calc_LCL(uf);
# LWP = incloud_LWP(uf, zb);
# OHU_400 = calc_OHU(uf,par,LWP,par.stype);

# # upsteps/downsteps
# CO2updn_list = [200, 400, 600, 800, 1000, 1200, 1600, 2000, 2400, 2800, 3200, 3600,
#                 3600, 3200, 2800, 2400, 2000, 1600, 1200, 1000, 800, 600, 400, 200];
# I = 12;
# par.stype = varSST();
# for (i,newCO2) in enumerate(CO2updn_list)
#     par.CO2 = newCO2;
#     par.OHU = OHU_400;
    
#     global u0, sol = run_mlm_from_init(uf, par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
    
#     # plot
#     mkpath(replace(path, "output"=>"figures"));
#     if i > I
#         filename = replace(path, "output"=>"figures")*"down"*string(Int(newCO2))*"_t.png";
#     else
#         filename = replace(path, "output"=>"figures")*"up"*string(Int(newCO2))*"_t.png";
#     end
#     plot_sol(sol, filename);

#     # print
#     global uf = sol.u[end];
#     zi, sM, qM, SST, CF = uf;
#     println(newCO2, ": ", CF)

#     # save
#     du = zeros(5);
#     mlm(du, uf, par, 0.0);
#     zb = calc_LCL(uf);
#     LWP = incloud_LWP(uf, zb);
#     RH = min(qM / q_sat(0.0, temp(0.0, sM, qM)), 1.0);
#     output = Dict("p" => par, "u0" => u0, "uf" => uf, "du/u" => du./uf, 
#     "we" => we(uf,par,zb,LWP,par.etype), "zb" => zb, "zc" => zi-zb,
#     "RHsurf" => RH, "LHF" => calc_LHF(uf,par), "SHF" => calc_SHF(uf,par),
#     "ΔR" => calc_cloudtop_RAD(uf,par,LWP,par.rtype), 
#     "OHU" => calc_OHU(uf,par,LWP,par.stype))
#     if i > I
#         save(path*"co2_downstep_"*string(Int(newCO2))*".jld2", output)
#     else
#         save(path*"co2_upstep_"*string(Int(newCO2))*".jld2", output)
#     end
# end

# plot
co2u = "[200, 400, 600, 800, 1000, 1200, 1600, 2000, 2400, 2800, 3200, 3600]";
co2d = "[3600, 3200, 2800, 2400, 2000, 1600, 1200, 1000, 800, 600, 400, 200]";
ARGS = [exp_path, co2u, co2d];
xtks = ([0, 1000, 2000, 3000])
xrange = [0,3800]
SSTrange = [285, 310]
LHFrange = [50, 250]

orig_co2u = "[200, 300, 400, 600, 800, 1000, 1200, 1300, 1400, 1600]";
orig_co2d = "[1600, 1400, 1300, 1200, 1000, 800, 600, 400, 300, 200]";
orig_path = "20221211_calibparam/";
title = "fixWV"
include("plot_hysteresis_noLES.jl")
