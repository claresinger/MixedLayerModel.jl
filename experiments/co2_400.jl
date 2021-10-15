using MixedLayerModel
using FileIO
using Plots

include("mlm_solve_funcs.jl")

# define path to save file (which experiment are you running?)
path = "experiments/output/test_cloudfrac/";

# define OHU from 400 ppm simulation
par = upCO2();
par.etype = enBal();
par.fttype = co2dep();
dt = 4.0;
tmax = 25.0;

u0, sol = run_mlm_ss(par, dt=3600.0*dt, tspan=3600.0*24.0*tmax);
code = sol.retcode;
println(code);
println(u0);

uf = sol.u;
du = zeros(5);
mlm(du, uf, par, 0.0);
zi,hM,qM,SST,CF = uf;
zb = calc_LCL(uf);
println(uf);
println(du);

output = Dict("code" => code, "p"=>par, "u0" => u0, "uf" => uf, "du/u" => du./uf, 
"we" => we(uf,par,par.etype), "zb" => zb, "zc" => zi-zb,
"RHsurf" => RH(0.0, hM, qM), "LHF" => calc_LHF(uf,par), "SHF" => calc_SHF(uf,par),
"ΔR" => calc_cloudtop_RAD(uf,par,par.rtype), "OHU" => calc_OHU(uf,par,par.stype));
save(path*"co2_400.jld2", output);

ENV["GKSwstype"]="nul"
u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));
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
savefig(replace(path, "output"=>"figures")*"sol400_t.png");