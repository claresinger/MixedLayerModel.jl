using MixedLayerModel
using FileIO

include("mlm_solve_funcs.jl")
include("plot_transient_solution.jl")

println()
# use command line argument to set co2
newCO2 = parse(Float64,ARGS[1]);
println(newCO2);
exp_path = ARGS[2];
path = "experiments/output/"*exp_path;

# par = upCO2();
par.CO2 = newCO2;
# par.etype = enBal();
# par.fttype = co2EIS();
# par.rtype = varRad();
# par.stype = varSST();
# dt, tmax = 48.0, 50.0;

# load initial condition from file
restarttry1 = path*"co2_downstep_"*string(Int(newCO2+100))*".jld2";
restarttry2 = path*"co2_downstep_"*string(Int(newCO2+200))*".jld2";
restarttry3 = path*"co2_downstep_"*string(Int(newCO2+400))*".jld2";
restarttry4 = path*"co2_downstep_"*string(Int(newCO2+800))*".jld2";
restarttry5 = path*"co2_upstep_"*string(Int(newCO2))*".jld2";
if isfile(restarttry1)
    output = load(restarttry1);
elseif isfile(restarttry2)
    output = load(restarttry2);
elseif isfile(restarttry3)
    output = load(restarttry3);
elseif isfile(restarttry4)
    output = load(restarttry4);
elseif isfile(restarttry5)
    output = load(restarttry5);
else
    println("no restart file")
end
u0 = output["uf"];
OHU = output["OHU"];
println("restarting from CO2 = "*string(output["p"].CO2));

# set OHU
par.OHU = OHU;

# solve and plot
ENV["GKSwstype"]="nul"
u0, sol = run_mlm_from_init(u0, par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
mkpath(replace(path, "output"=>"figures"));
filename = replace(path, "output"=>"figures")*"down"*string(Int(newCO2))*"_t.png";
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

output = Dict("p" => par, "u0" => u0, "uf" => uf, "du/u" => du./uf, 
"we" => we(uf,par,zb,LWP,par.etype), "zb" => zb, "zc" => zi-zb,
"RHsurf" => RH, "LHF" => calc_LHF(uf,par), "SHF" => calc_SHF(uf,par),
"ΔR" => calc_cloudtop_RAD(uf,par,LWP,par.rtype), 
"OHU" => calc_OHU(uf,par,LWP,par.stype))

save(path*"co2_downstep_"*string(Int(newCO2))*".jld2", output)
