@everywhere module GModel

using Distributed
using MixedLayerModel
include("mlm_solve_funcs.jl")
include("normalize.jl")

export run_ensembles
export run_forward

function run_ensembles(params, N_ens)
    # g_ens = zeros(nd*2, N_ens)
    # for i in 1:N_ens
    #     g_ens[:, i] = run_forward(params[:,i]) # map θ to G(θ)
    # end

    g_ens = zeros(nd*2, N_ens)
    params = [params[:,i] for i in 1:size(params,2)]
    g_ens[:, :] = hcat(pmap(x -> run_forward(x), params)...) # map θ to G(θ)

    return g_ens
end

function run_forward(params)
    # create parameters
    par = upCO2();
    par.etype = enBal();
    par.fttype = co2EIS();
    par.rtype = varRad();
    par.stype = fixSST();
    dt, tmax = 10.0, 100.0; # days

    # adjust tunable parameters
    # par.decoup_slope = params[1];
    par.Cd = params[1];
    par.α_vent = params[2];
    # par.EIS0 = params[3];
    # par.ECS = params[4];
    # par.Eexport = params[5];
    # par.SW_b = params[6];

    # 400 ppm
    u0, sol = run_mlm(par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
    uf = sol.u[end];
    zb = calc_LCL(uf);
    LWP = incloud_LWP(uf, zb);
    OHU_400 = calc_OHU(uf,par,LWP,par.stype);

    # upsteps
    par.stype = varSST();
    Gstacked = zeros(nd*2, 1);
    try
        for (i,newCO2) in enumerate(CO2updn_list)
            par.CO2 = newCO2;
            par.OHU = OHU_400;
            
            u0, sol = run_mlm_from_init(uf, par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
            uf = sol.u[end];
            
            zi, sM, qM, SST, CF = uf;
            Gstacked[i,:] .= normalize_data(SST, "SST");
            Gstacked[nd+i,:] .= normalize_data(calc_LHF(uf, par), "LHF");
        end
    catch
        println("catching error in MLM")
        println(Gstacked)
        fill!(Gstacked, NaN)
    end
    return Gstacked
end

end