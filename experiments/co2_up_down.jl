exp_path = "20230407_OU-ensemble10_noise3_10days/";
path = "experiments/output/"*exp_path;

using MixedLayerModel
using JLD2, Statistics
using SciMLBase.EnsembleAnalysis
include("mlm_solve_funcs.jl")
include("plot_transient_solution.jl")

# create parameters
par = upCO2();
par.etype = enBal();
par.fttype = co2EIS();
par.rtype = varRad();
par.stype = fixSST();
dt, tmax = 1/6, 100.0; # days
nhalf = Int(tmax/dt/2);

# adjust tunable parameters
par.Cd = 8e-4; #7.9e-4;
par.α_vent = 1.5e-3; #1.69e-3;
par.SW_b = 120;
par.Dcrit = 1.1;
par.noise = 3e-3;

# 400 ppm
# u0, sol = run_mlm(par, dt=3600.0*24.0*dt*100, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
# uf = sol.u[end];
u0, sim, sol = stochastic_run_mlm(par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
u = timeseries_steps_mean(sim);
uf = dropdims(mean(u[end-nhalf:end], dims=2), dims=2);
println("init 400: ", uf[5])
zb = calc_LCL(uf);
LWP = incloud_LWP(uf, zb);
OHU_400 = calc_OHU(uf,par,LWP,par.stype);
println("OHU = ", OHU_400)

mkpath(replace(path, "output"=>"figures"));
plot_sol(sol, replace(path, "output"=>"figures")*"init400.png", ndims(u))

# upsteps/downsteps
CO2updn_list = [200,300,400,600,800,1000,1200,1300,1400,1600,1600,1400,1300,1200,1000,800,600,400,300,200];
I = 10;
par.stype = varSST();
for (i,newCO2) in enumerate(CO2updn_list)
    par.CO2 = newCO2;
    par.OHU = OHU_400;

    # local u0, sol = run_mlm_from_init(uf, par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
    local u0, sim, sol = stochastic_run_mlm_from_init(uf, par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
    local u = timeseries_steps_mean(sim);

    # plot
    if i > I
        filename = replace(path, "output"=>"figures")*"down"*string(Int(newCO2))*"_t.png";
    else
        filename = replace(path, "output"=>"figures")*"up"*string(Int(newCO2))*"_t.png";
    end
    plot_sol(sol, filename, ndims(u));

    # print
    global uf = dropdims(mean(u[end-nhalf:end], dims=2), dims=2);
    local zi, sM, qM, SST, CF = uf;
    local zb = calc_LCL(uf);
    local lwp = incloud_LWP(uf, zb);
    println(newCO2, ": ", CF)

    # save
    du = zeros(5);
    mlm(du, uf, par, 0.0);
    output = Dict(
        "p" => par, 
        "u0" => u0, 
        "uf" => uf, 
        "du/u" => du./uf, 
        "LWP" => lwp,
        "zb" => zb, 
        "we" => we(uf, par, zb, lwp, par.etype), 
        "LHF" => calc_LHF(uf, par), 
        "SHF" => calc_SHF(uf, par),
        "ΔR" => calc_cloudtop_RAD(uf, par, lwp, par.rtype), 
        "OHU" => calc_OHU(uf, par, lwp, par.stype),
        "De" => calc_decoupling(uf, par, zb, lwp),
    )
    if i > I
        save(path*"co2_downstep_"*string(Int(newCO2))*".jld2", output)
    else
        save(path*"co2_upstep_"*string(Int(newCO2))*".jld2", output)
    end
end

# plot
co2u = "[200, 300, 400, 600, 800, 1000, 1200, 1300, 1400, 1600]";
co2d = "[1600, 1400, 1300, 1200, 1000, 800, 600, 400, 300, 200]";
ARGS = [exp_path, co2u, co2d];
include("plot_hysteresis.jl")
