@everywhere module GModel

using Distributed

using MixedLayerModel
include("mlm_solve_funcs.jl")

export run_ensembles
export run_forward

CO2updn_list = [200,300,400,800,1000,1200,1300,1400,1600,1400,1300,1200,1000,800,400,300,200];
nd = length(CO2updn_list)

function run_ensembles(params, N_ens)
    # g_ens = zeros(nd, N_ens)
    # for i in 1:N_ens
    #     # run the model with the current parameters, i.e., map θ to G(θ)
    #     g_ens[:, i] = run_forward(params[:,i])
    # end

    g_ens = zeros(nd, N_ens)
    params = [params[:,i] for i in 1:size(params,2)]
    g_ens[:, :] = hcat(pmap(x -> run_forward(x), params)...)
    return g_ens
end

function run_forward(params)
    # create parameters
    par = upCO2();
    par.etype = enBal();
    par.fttype = co2EIS();
    par.rtype = varRad();
    par.stype = fixSST();
    dt, tmax = 48.0, 50.0;

    # adjust tunable parameters
    par.decoup_slope = params[1]; #10;
    par.α_vent = params[2]; #1.08e-3;
    par.EIS0 = params[3]; #10.0;
    par.ECS = params[4]; #3.0;
    par.Eexport = params[5]; #15.0;
    par.SW_b = params[6]; #150;

    # 400 ppm
    u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
    uf = sol.u[end];
    zb = calc_LCL(uf);
    LWP = incloud_LWP(uf, zb);
    OHU_400 = calc_OHU(uf,par,LWP,par.stype);

    # upsteps
    par.stype = varSST();
    Gstacked = zeros(nd, 1);

    for (i,newCO2) in enumerate(CO2updn_list)
        par.CO2 = newCO2;
        par.OHU = OHU_400;
        u0, sol = run_mlm_from_init(uf, par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
        uf = sol.u[end];
        zi, sM, qM, SST, CF = uf;
        Gstacked[i,:] .= SST;
    end

    return Gstacked
end

end