using MixedLayerModel
using FileIO
using Plots

include("mlm_solve_funcs.jl")

# use command line argument to set co2
newCO2 = parse(Float64,ARGS[1]);
println(newCO2);

# load initial condition from file
path = "experiments/output/new_alpha_enBal_invco2/";
output = load(path*"co2_400.jld2");
u0 = output["uf"];
OHU = output["OHU"];

# set OHU, increase CO2, let SST evolve and check cloud changes
par = basic_params();
par.Hw = 0.1;
par.OHU = OHU;
par.CO2 = newCO2;
par.etype = enBal();
par.fttype = co2dep();
par.rtype = varRad();
par.stype = varSST();
dt, tmax = 2.0, 30.0;
u0, sol = run_mlm_ss_from_init(u0, par, dt=3600.0*dt, tspan=3600.0*24.0*tmax);
code = sol.retcode;
println(code);

uf = sol.u;
du = zeros(4);
mlm(du, uf, par, 0.0);
zi,hM,qM,SST = uf;
zb = calc_LCL(zi,hM,qM);
println(uf);
println(du);

output = Dict("code" => code, "p" => par, "u0" => u0, "uf" => uf, "du/u" => du./uf, 
"we" => we(uf,par,par.etype), "zb" => zb, "zc" => zi-zb,
"RHsurf" => RH(0.0, hM, qM), "LHF" => calc_LHF(uf,par), "SHF" => calc_SHF(uf,par),
"ΔR" => calc_cloudtop_RAD(uf,par,par.rtype), "OHU" => calc_OHU(uf,par,par.stype))

save(path*"co2_upstep_"*string(Int(newCO2))*".jld2", output)

ENV["GKSwstype"]="nul"
u0, sol = run_mlm_from_init(uf, par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
t = sol.t / 3600.0 / 24.0;
zi = getindex.(sol.u,1);
hM = getindex.(sol.u,2) * 1e-3;
qtM = getindex.(sol.u,3) * 1e3;
sst = getindex.(sol.u,4);
plot(size=(600,500), layout=(4,1), dpi=200);
plot!(t, zi, legend=:topleft, subplot=1, label="zi(t) [m]");
plot!(t, hM, legend=:topleft, subplot=2, label="hM(t) [kJ/kg]"); 
plot!(t, qtM, legend=:topleft, subplot=3, label="qtM(t) [g/kg]");
plot!(t, sst, legend=:topleft, subplot=4, label="SST(t) [K]");
mkpath(replace(path, "output"=>"figures"));;
savefig(replace(path, "output"=>"figures")*"sol"*string(Int(newCO2))*"_t.png");