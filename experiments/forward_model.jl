@everywhere module GModel

using Distributed
using Statistics

using MixedLayerModel
include("mlm_solve_funcs.jl")

export run_ensembles
export run_forward
export normalize_data, unnormalize_data

CO2updn_list = [200,300,400,800,1000,1200,1300,1400,1600,1400,1300,1200,1000,800,400,300,200];
nd = length(CO2updn_list)

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
    par.EIS0 = params[3];
    par.ECS = params[4];
    par.Eexport = params[5];
    par.SW_b = params[6];

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
            
            u0, sol = run_mlm_from_init(uf, par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
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
    # if (maximum(Gstacked) > 1e5) | (minimum(Gstacked) < -1e5)
    #     println(Gstacked)
    #     println("NaN fill, fit too terrible: ", maximum(Gstacked), "\t", minimum(Gstacked))
    #     fill!(Gstacked, NaN)
    # end
    return Gstacked
end

SSTupdn_list = [287.7,289.1,290.0,292.2,293.2,294.3,304.5,305.8,308.0,306.2,304.7,303.7,302.0,300.9,297.9,296.8,287.6];
LHFupdn_list = [97.8,103.9,107.1,112.8,115.3,120.5,208.7,213.5,220.9,214.4,209.4,206.0,199.7,195.2,183.7,179.1,96.1];
ziupdn_list = [1442,1349,1266,1078,1011,972,834,781,703,766,826,866,949,1013,1230,1340,1416];
μSST, σSST = mean(SSTupdn_list), std(SSTupdn_list);
μLHF, σLHF = mean(LHFupdn_list), std(LHFupdn_list);

function normalize_data(x, name="SST")
    if name == "SST"
        x = (x .- μSST) ./ σSST;
    end
    if name == "LHF"
        x = (x .- μLHF) ./ σLHF;
    end
    return x
end

function unnormalize_data(x, name="SST")
    if name == "SST"
        x = x .* σSST .+ μSST;
    end
    if name == "LHF"
        x = x .* σLHF .+ μLHF;
    end
    return x
end

end