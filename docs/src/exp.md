# Running an experiment

First use the MLM to calculate the steady-state clouds given an atmosphere with 400 ppm CO_2.
Save this result to a file. In this example we use the energy balance entrainment parameterization and have interactive radiation, but fixed sea-surface temperatures.
```julia
using MixedLayerModel
using FileIO

# run simulation with 400 ppm CO2
par = basic_params();
par.etype = enBal();
par.rtype = varRad();
u0, sol = run_mlm(par);

# get output
code = sol.retcode;
println(code);
uf = sol.u;
du = zeros(4);
mlm(du, uf, par, 0.0);
zi,hM,qM,SST = uf;
zb = calc_LCL(zi,hM,qM);

# save output
output = Dict("code" => code, "p"=>par, "u0" => u0, "uf" => uf, "du/u" => du./uf, 
"we" => we(uf,par,par.etype), "zb" => zb, "zc" => zi-zb,
"RHsurf" => RH(0.0, hM, qM), "LHF" => calc_LHF(uf,par), "SHF" => calc_SHF(uf,par),
"ΔR" => calc_cloudtop_RAD(uf,par,par.rtype), "OHU" => calc_OHU(uf,par,par.stype))
path = "experiments/output/exp_name/";
save(path*"co2_400.jld2", output)
```

Now we can load in the results from the 400 ppm simulation and use the calculated steady-state to initialize a new simulation with 800 ppm CO_2.
In this example we still have fixed SSTs. 
```julia
using MixedLayerModel
using FileIO

# load initial condition from file
path = "experiments/output/exp_name/";
output = load(path*"co2_400.jld2");
u0 = output["uf"];

# run experiment now with 800 ppm CO2
par = basic_params();
newCO2 = 800.0;
par.CO2 = newCO2;
par.etype = enBal();
par.rtype = varRad();
u0, sol = run_mlm_from_init(u0, par);

# get output
code = sol.retcode;
println(code);
uf = sol.u;
du = zeros(4);
mlm(du, uf, par, 0.0);
zi,hM,qM,SST = uf;
zb = calc_LCL(zi,hM,qM);

# save output
output = Dict("code" => code, "p" => par, "u0" => u0, "uf" => uf, "du/u" => du./uf, 
"we" => we(uf,par,par.etype), "zb" => zb, "zc" => zi-zb,
"RHsurf" => RH(0.0, hM, qM), "LHF" => calc_LHF(uf,par), "SHF" => calc_SHF(uf,par),
"ΔR" => calc_cloudtop_RAD(uf,par,par.rtype), "OHU" => calc_OHU(uf,par,par.stype))
save(path*"co2_upstep_fixSST_"*string(Int(newCO2))*".jld2", output)
```

Instead if we want to allow SSTs to change interactively, then we can set the ocean heat uptake (OHU) to the value defined by the 400 ppm simulation.
Also be sure to switch to ``varSST()`` mode.
```julia
using MixedLayerModel
using FileIO

# load initial condition from file
path = "experiments/output/exp_name/";
output = load(path*"co2_400.jld2");
u0 = output["uf"];
OHU = output["OHU"];

# run experiment now with 800 ppm CO2
par = basic_params();
newCO2 = 800.0;
par.CO2 = newCO2;
par.OHU = OHU;
par.etype = enBal();
par.rtype = varRad();
par.stype = varSST();
u0, sol = run_mlm_from_init(u0, par);

# get output
code = sol.retcode;
println(code);
uf = sol.u;
du = zeros(4);
mlm(du, uf, par, 0.0);
zi,hM,qM,SST = uf;
zb = calc_LCL(zi,hM,qM);

# save output
output = Dict("code" => code, "p" => par, "u0" => u0, "uf" => uf, "du/u" => du./uf, 
"we" => we(uf,par,par.etype), "zb" => zb, "zc" => zi-zb,
"RHsurf" => RH(0.0, hM, qM), "LHF" => calc_LHF(uf,par), "SHF" => calc_SHF(uf,par),
"ΔR" => calc_cloudtop_RAD(uf,par,par.rtype), "OHU" => calc_OHU(uf,par,par.stype))
save(path*"co2_upstep_"*string(Int(newCO2))*".jld2", output)
```