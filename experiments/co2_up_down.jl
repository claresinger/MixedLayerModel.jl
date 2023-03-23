exp_path = "20230322_newD/";
path = "experiments/output/"*exp_path;

using MixedLayerModel
using JLD2
include("mlm_solve_funcs.jl")
include("plot_transient_solution.jl")

# create parameters
par = upCO2();
par.etype = enBal();
par.fttype = co2EIS();
par.rtype = varRad();
par.stype = fixSST();
dt, tmax = 10.0, 100.0; # days

# adjust tunable parameters
par.Cd = 8e-4; #7.9e-4;
par.α_vent = 1e-3; #1.69e-3;
par.SW_b = 120;
par.Dcrit = 1;

# 400 ppm
u0, sol = run_mlm(par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
uf = sol.u[end];
zb = calc_LCL(uf);
LWP = incloud_LWP(uf, zb);
OHU_400 = calc_OHU(uf,par,LWP,par.stype);
println(OHU_400)
plot_sol(sol,replace(path, "output"=>"figures")*"init400.png";)

# upsteps/downsteps
CO2updn_list = [200,300,400,600,800,1000,1200,1300,1400,1600,1600,1400,1300,1200,1000,800,600,400,300,200];
I = 10;
par.stype = varSST();
for (i,newCO2) in enumerate(CO2updn_list)
    par.CO2 = newCO2;
    par.OHU = OHU_400;
    
    local u0, sol = run_mlm_from_init(uf, par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
    
    # plot
    mkpath(replace(path, "output"=>"figures"));
    if i > I
        filename = replace(path, "output"=>"figures")*"down"*string(Int(newCO2))*"_t.png";
    else
        filename = replace(path, "output"=>"figures")*"up"*string(Int(newCO2))*"_t.png";
    end
    plot_sol(sol, filename);

    # print
    uarr = sol.u[end-100:end];
    global uf = mean(uarr);
    local zi, sM, qM, SST, CF = uf;
    local zb = calc_LCL.(uarr);
    local lwp = incloud_LWP.(uarr, zb);
    println(newCO2, ": ", CF)

    # save
    du = zeros(5);
    mlm(du, uf, par, 0.0);
    output = Dict(
        "p" => par, 
        "u0" => u0, 
        "uf" => uf, 
        "du/u" => du./uf, 
        "LWP" => mean(lwp),
        "zb" => mean(zb), 
        "zc" => zi-mean(zb),
        "we" => mean(we.(uarr,Ref(par), zb, lwp, Ref(par.etype))), 
        # "RHsurf" => min(qM / q_sat(0.0, temp(0.0, sM, qM)), 1.0), 
        "LHF" => mean(calc_LHF.(uarr, Ref(par))), 
        "SHF" => mean(calc_SHF.(uarr, Ref(par))),
        "ΔR" => mean(calc_cloudtop_RAD.(uarr, Ref(par), lwp, Ref(par.rtype))), 
        "OHU" => mean(calc_OHU.(uarr, Ref(par), lwp, Ref(par.stype))),
        "De" => mean(calc_decoupling.(uarr, Ref(par), zb, lwp)),
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
