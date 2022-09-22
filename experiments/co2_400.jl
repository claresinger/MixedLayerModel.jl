push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using MixedLayerModel
using FileIO

include("mlm_solve_funcs.jl")
include("plot_transient_solution.jl")

# define path to save file (which experiment are you running?)
exp_path = ARGS[1];
path = "experiments/output/"*exp_path;

# par = upCO2();
# par.etype = enBal();
# par.fttype = co2EIS();
# par.rtype = varRad();
# par.stype = fixSST();
# dt, tmax = 48.0, 50.0;

## solve and plot
ENV["GKSwstype"]="nul"
u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
mkpath(replace(path, "output"=>"figures"));
filename = replace(path, "output"=>"figures")*"sol400_t.png";
plot_sol(sol, filename);

## save steady-state solution
uf = sol.u[end];
du = zeros(5);
mlm(du, uf, par, 0.0);
zi,sM,qM,SST,CF = uf;
zb = calc_LCL(uf);
LWP = incloud_LWP(uf, zb);
RH = min(qM / q_sat(0.0, temp(0.0, sM, qM)), 1.0);
println(uf);
println(du);
println("cloud base: ",zb)
println("LWP: ", LWP);
println("tropical sst: ", trop_sst(uf, par, LWP));
println("ft T: ", temp_ft(trop_sst(uf, par, LWP), zi, par));
println("cloud top T: ", temp(zi, sM, qM));
println("ft qt: ", qjump(uf, par, LWP, par.fttype) + qM);

output = Dict("p" => par, "u0" => u0, "uf" => uf, "du/u" => du./uf, 
"we" => we(uf,par,zb,LWP,par.etype), "zb" => zb, "zc" => zi-zb,
"RHsurf" => RH, "LHF" => calc_LHF(uf,par), "SHF" => calc_SHF(uf,par),
"ΔR" => calc_cloudtop_RAD(uf,par,LWP,par.rtype), 
"OHU" => calc_OHU(uf,par,LWP,par.stype))

save(path*"co2_400.jld2", output);