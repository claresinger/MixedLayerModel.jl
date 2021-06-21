using MixedLayerModel
using FileIO

# use command line argument to set co2
newCO2 = parse(Float64,ARGS[1]);
println(newCO2);

# load initial condition from file
path = "experiments/output/test_deep_ocean/";
output = load(path*"co2_400.jld2");
u0 = output["uf"];
OHU = output["OHU"];

# set OHU, increase CO2, let SST evolve and check cloud changes
par = basic_params();
par.Hw = 10;
par.OHU = OHU;
par.CO2 = newCO2;
par.etype = enBal();
par.rtype = varRad();
par.stype = varSST();
u0, sol = run_mlm_from_init(u0, par);

code = sol.retcode;
println(code);

uf = sol.u;
du = zeros(4);
mlm(du, uf, par, 0.0);
zi,hM,qM,SST = uf;
zb = calc_LCL(zi,hM,qM);

output = Dict("code" => code, "p" => par, "u0" => u0, "uf" => uf, "du/u" => du./uf, 
"we" => we(uf,par,par.etype), "zb" => zb, "zc" => zi-zb,
"RHsurf" => RH(0.0, hM, qM), "LHF" => calc_LHF(uf,par), "SHF" => calc_SHF(uf,par),
"ΔR" => calc_cloudtop_RAD(uf,par,par.rtype), "OHU" => calc_OHU(uf,par,par.stype))

save(path*"co2_upstep_"*string(Int(newCO2))*".jld2", output)