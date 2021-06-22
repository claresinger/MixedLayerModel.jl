using MixedLayerModel
using FileIO

# define OHU from 400 ppm simulation
par = basic_params();
par.etype = Sally();
par.a = 1.1;
u0, sol = run_mlm(par);
code = sol.retcode;
println(code);

uf = sol.u;
du = zeros(4);
mlm(du, uf, par, 0.0);
zi,hM,qM,SST = uf;
zb = calc_LCL(zi,hM,qM);
println(uf);
println(du);

output = Dict("code" => code, "p"=>par, "u0" => u0, "uf" => uf, "du/u" => du./uf, 
"we" => we(uf,par,par.etype), "zb" => zb, "zc" => zi-zb,
"RHsurf" => RH(0.0, hM, qM), "LHF" => calc_LHF(uf,par), "SHF" => calc_SHF(uf,par),
"ΔR" => calc_cloudtop_RAD(uf,par,par.rtype), "OHU" => calc_OHU(uf,par,par.stype));
path = "experiments/output/Sally1.1/";
save(path*"co2_400.jld2", output);