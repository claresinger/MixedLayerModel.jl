push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using MixedLayerModel
using FileIO
using Plots

include("mlm_solve_funcs.jl")
include("plot_transient_solution.jl")

# use command line argument to set co2
newCO2 = parse(Float64,ARGS[1]);
# newCO2 = 200.0;
println(newCO2);

par = upCO2();
par.CO2 = newCO2;
par.etype = enBal();
par.fttype = co2dep();
par.rtype = varRad();
par.stype = varSST();
dt, tmax = 48.0, 50.0;

# load initial condition from file
path = "experiments/output/cfmip_modCF_newCF/";
restarttry1 = path*"co2_downstep_"*string(Int(newCO2+50))*".jld2";
restarttry2 = path*"co2_downstep_"*string(Int(newCO2+100))*".jld2";
restarttry3 = path*"co2_downstep_"*string(Int(newCO2+200))*".jld2";
if isfile(restarttry1)
    output = load(restarttry1);
elseif isfile(restarttry2)
    output = load(restarttry2);
elseif isfile(restarttry3)
    output = load(restarttry3);
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
println("cloud base: ",zb)
println("LWP: ", LWP);
println("tropical sst: ", trop_sst(uf, par, LWP));
println("ft qt: ", qjump(uf, par, LWP, par.fttype) + qM);

output = Dict("p" => par, "u0" => u0, "uf" => uf, "du/u" => du./uf, 
"we" => we(uf,par,zb,LWP,par.etype), "zb" => zb, "zc" => zi-zb,
"RHsurf" => RH, "LHF" => calc_LHF(uf,par), "SHF" => calc_SHF(uf,par),
"ΔR" => calc_cloudtop_RAD(uf,par,LWP,par.rtype), 
"OHU" => calc_OHU(uf,par,LWP,par.stype))

save(path*"co2_downstep_"*string(Int(newCO2))*".jld2", output)

# ### AGU plots
# Plots.scalefontsizes(2)
# plot(size=(1000,500), layout=(2,2), dpi=200, left_margin = 5Plots.mm, bottom_margin=5Plots.mm);
# plot!(t, ΔR, marker="o-", legend=false, subplot=1, ylabel="ΔR [W/m\$^2\$]");
# plot!(t, S, marker="o-", legend=false, subplot=2, ylabel="Stability, \$S\$",
#         yscale=:log10, yticks=([0.2,0.5,2,5], ["0.2","0.5","2","5"]), ylim=[0.2,10]);
# plot!(t, cf * 1e2, marker="o-", legend=false, subplot=3, ylabel="CF [%]", xlabel="Time [days]");
# plot!(t, sst, marker="o-", legend=false, subplot=4, ylabel="SST [K]", xlabel="Time [days]");
# mkpath(replace(path, "output"=>"figures"));
# savefig(replace(path, "output"=>"figures")*"AGU-down"*string(Int(newCO2))*"_t.png");
# Plots.scalefontsizes(1/2)
# ### 