exp_path = "20260208/";
path = "experiments/output/"*exp_path;

using MixedLayerModel
using JLD2
include("mlm_solve_funcs.jl")

do_calc = false
do_plot = true

expN = 50;
CO2updn_list = 400:100:2000;

if do_calc
    # create parameters
    par = upCO2();
    par.etype = enBal();
    par.fttype = co2EIS();
    par.rtype = varRad();
    dt, tmax = 10.0, 100.0; # days

    # adjust tunable parameters
    par.Cd = 7.9e-4;
    par.α_vent = 1.69e-3;
    par.SW_b = 140;

    for expi in 1:expN
        println(expi)
        try
            # reset initial conditions
            par.SST0 = 290.0;
            par.EIS0 = 8.0;
            par.V = 10.0;
            par.D = 6.0e-6;
            par.RHft = 0.2;

            # perturb initial conditions
            par.SST0 += 0.3*randn()
            par.EIS0 += 0.1*randn()
            par.V += 0.01*randn()
            par.D += 0.05e-6*randn()
            par.RHft += 0.01*randn()

            # 400 ppm
            par.stype = fixSST();
            par.CO2 = 400;
            u0, sol = run_mlm(par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
            uf = sol.u[end];
            zb = calc_LCL(uf);
            LWP = incloud_LWP(uf, zb);
            OHU_400 = calc_OHU(uf,par,LWP,par.stype);
            par.stype = varSST();
            par.OHU = OHU_400;

            # println(par)
            # println()

            # upsteps
            for (i,newCO2) in enumerate(CO2updn_list)
                par.CO2 = newCO2;
                u0, sol = run_mlm_from_init(uf, par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
                uf = sol.u[end];
                zi, sM, qM, SST, CF = uf;

                # # print
                # println(newCO2)
                # println("u0: ", uf)
                # println("uf: ", uf)
                # println()

                # save
                du = zeros(5);
                mlm(du, uf, par, 0.0);
                zb = calc_LCL(uf);
                LWP = incloud_LWP(uf, zb);
                RH = min(qM / q_sat(0.0, temp(0.0, sM, qM)), 1.0);
                output = Dict("p" => par, "u0" => u0, "uf" => uf, "du/u" => du./uf, 
                "we" => we(uf,par,zb,LWP,par.etype), "zb" => zb, "zc" => zi-zb,
                "RHsurf" => RH, "LHF" => calc_LHF(uf,par), "SHF" => calc_SHF(uf,par),
                "ΔR" => calc_cloudtop_RAD(uf,par,LWP,par.rtype), 
                "OHU" => calc_OHU(uf,par,LWP,par.stype))
                save(path*string(expi)*"co2_upstep_"*string(Int(newCO2))*".jld2", output)
            end
        catch
            println("fail")
        end
    end
end

if do_plot
    # plot
    ARGS = [exp_path, CO2updn_list, expN];
    include("plot_breakup_varICs.jl")
end