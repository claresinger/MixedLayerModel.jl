using MixedLayerModel
using FileIO
using Plots

include("mlm_solve_funcs.jl")

# use command line argument to set co2
newCO2 = parse(Float64,ARGS[1]);
println(newCO2);

# load initial condition from file
path = "experiments/output/test_cloudfrac/";
restarttry = path*"co2_upstep_"*string(Int(newCO2-200))*".jld2";
restarttry2 = path*"co2_upstep_"*string(Int(newCO2-400))*".jld2";
if isfile(restarttry)
    output = load(restarttry);
    println("found file with CO2 = "*string(Int(newCO2-200)));
elseif isfile(restarttry2)
    output = load(restarttry2);
    println("found file with CO2 = "*string(Int(newCO2-400)));
else
    output = load(path*"co2_400.jld2");
    println("no file with CO2 = "*string(Int(newCO2-200)));
end
u0 = output["uf"];
OHU = output["OHU"];
println("restarting from CO2 = "*string(output["p"].CO2));

# set OHU, increase CO2, let SST evolve and check cloud changes
par = upCO2();
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
du = zeros(5);
mlm(du, uf, par, 0.0);
zi,hM,qM,SST = uf;
zb = calc_LCL(uf);
println(uf);
println(du);

output = Dict("code" => code, "p" => par, "u0" => u0, "uf" => uf, "du/u" => du./uf, 
"we" => we(uf,par,par.etype), "zb" => zb, "zc" => zi-zb,
"RHsurf" => RH(0.0, hM, qM), "LHF" => calc_LHF(uf,par), "SHF" => calc_SHF(uf,par),
"ΔR" => calc_cloudtop_RAD(uf,par,par.rtype), "OHU" => calc_OHU(uf,par,par.stype))

save(path*"co2_upstep_"*string(Int(newCO2))*".jld2", output)

ENV["GKSwstype"]="nul"
u0, sol = run_mlm_from_init(u0, par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
t = sol.t / 3600.0 / 24.0;
zi = getindex.(sol.u,1);
hM = getindex.(sol.u,2) * 1e-3;
qtM = getindex.(sol.u,3) * 1e3;
sst = getindex.(sol.u,4);
cf = getindex.(sol.u,5) * 1e2;
S = zeros(length(t));
for (i,si) in enumerate(S)
    S[i] = calc_S(sol.u[i], par);
end
plot(size=(600,800), layout=(6,1), dpi=200, left_margin = 5Plots.mm);
plot!(t, zi, legend=false, subplot=1, ylabel="zi [m]");
plot!(t, hM, legend=false, subplot=2, ylabel="hM [kJ/kg]"); 
plot!(t, qtM, legend=false, subplot=3, ylabel="qtM [g/kg]");
plot!(t, sst, legend=false, subplot=4, ylabel="SST [K]");
plot!(t, cf, legend=false, subplot=5, ylabel="CF [%]");
plot!(t, S, legend=false, subplot=6, ylabel="S [-]", xlabel="time [days]");
mkpath(replace(path, "output"=>"figures"));
savefig(replace(path, "output"=>"figures")*"sol"*string(Int(newCO2))*"_t.png");
