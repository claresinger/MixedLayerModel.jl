push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using MixedLayerModel
using Plots

function plot_sol(sol, filename)
    ## plot and save
    t = sol.t / 3600.0 / 24.0;
    zi = getindex.(sol.u,1);
    sM = getindex.(sol.u,2);
    qtM = getindex.(sol.u,3);
    sst = getindex.(sol.u,4);
    cf = getindex.(sol.u,5);
    S = zeros(length(t));
    LHF = zeros(length(t));
    zb = zeros(length(t));
    ΔR = zeros(length(t));
    LWP = zeros(length(t));
    Δsv = zeros(length(t));
    ent = zeros(length(t));
    tropSST = zeros(length(t));
    for (i,si) in enumerate(S)
        zb[i] = calc_LCL(sol.u[i]);
        LWP[i] = incloud_LWP(sol.u[i], zb[i]);
        S[i] = calc_S(sol.u[i], par, zb[i], LWP[i]);
        LHF[i] = calc_LHF(sol.u[i], par);
        ΔR[i] = calc_cloudtop_RAD(sol.u[i], par, LWP[i], par.rtype);
        Δsv[i] = sv_jump(sol.u[i], par, LWP[i]);
        ent[i] = we(sol.u[i], par, zb[i], LWP[i], par.etype);
        tropSST[i] = trop_sst(sol.u[i], par, LWP[i]);
    end 
    plot(size=(1200,800), layout=(6,2), dpi=200, left_margin = 5Plots.mm);
    plot!(t, zi, line=2, marker=:circle, legend=false, subplot=1, ylabel="zi, zb [m]");
    plot!(t, zb, line=2, marker=:circle, legend=false, subplot=1);
    plot!(t, sM * 1e-3, line=2, marker=:circle, legend=false, subplot=2, ylabel="sM [kJ/kg]"); 
    plot!(t, qtM * 1e3, line=2, marker=:circle, legend=false, subplot=3, ylabel="qtM [g/kg]");
    plot!(t, sst, line=2, marker=:circle, legend=false, subplot=4, ylabel="SST [K]");
    plot!(t, tropSST, line=2, marker=:circle, legend=false, subplot=4);
    plot!(t, cf * 1e2, line=2, marker=:circle, legend=false, subplot=5, ylabel="CF [%]");
    plot!(t, LWP .* cf * 1e3, line=2, marker=:circle, legend=false, subplot=6, ylabel="LWP [g/m2]");
    plot!(t, Δsv * 1e-3, line=2, marker=:circle, legend=false, subplot=7, ylabel="Δsv (kJ/kg)");
    plot!(t, ent*1e3, line=2, marker=:circle, legend=false, subplot=8, ylabel="we (mm/s)")
    plot!(t, LHF, line=2, marker=:circle, legend=false, subplot=9, ylabel="LHF [W/m2]");
    plot!(t, ΔR, line=2, marker=:circle, legend=false, subplot=10, ylabel="ΔR [W/m2]");
    plot!(t, (zi .- zb) ./ zi, line=2, marker=:circle, legend=false, subplot=11, ylabel="zc/zi [-]", xlabel="time [days]");
    plot!(t, S, line=2, marker=:circle, legend=false, subplot=12, ylabel="S [-]", xlabel="time [days]");
    savefig(filename);

    ### AGU plots
    if (cf[1] > 0.5) && (cf[end] < 0.5)
        Plots.scalefontsizes(1.5)
        plot(size=(1000,500), layout=(2,2), dpi=200, left_margin = 5Plots.mm, bottom_margin=10Plots.mm);
        plot!(t, ΔR, marker="o-", legend=false, subplot=1, 
                ylabel="ΔR [W/m\$^2\$]", ylim=[0,60]);
        plot!(t, sst, marker="o-", legend=false, subplot=2, 
                ylabel="SST [K]", ylim=[290,315]);
        plot!(t, S, marker="o-", legend=false, subplot=3, xlabel="Time [days]",
                ylabel="Decoupling, \$\\mathcal{D}\$", yscale=:log10,
                yticks=([0.3,1,3], ["0.3","1","3"]), ylim=[0.25,3]);
        plot!(t, cf * 1e2, marker="o-", legend=false, subplot=4, 
                ylabel="CF [%]", xlabel="Time [days]", ylim=[0,100]);
        savefig(replace(filename, "up"=>"AGU-up"))
        Plots.scalefontsizes(1/1.5)
    end
    ### 
end