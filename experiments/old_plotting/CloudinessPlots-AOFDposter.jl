module CloudinessPlots

###########
# import functions and constants from other files
###########
using PyPlot
using LaTeXStrings

include("AtmosParams.jl")
include("TDFunctions.jl")
using .AtmosParams
using .TDFunctions

###########
# export functions for other files
###########
export plot_mult_lwp_var

###########
# functions
###########
function plot_mult_lwp_var(u_arr, par_arr, var_arr,save)
    colors = ["#FF6C0C","#73A950","#00A1DF"];
    symbols = ["o","s","^"];

    figure(figsize=(8,5));
    ax1 = subplot(111);
    ax2 = ax1.twiny();
    ax3 = ax1.twiny();
    p1,p2,p3 = 0,0,0;

    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
    rcParams["font.size"] = 15;

    for (i,var) in enumerate(var_arr)
        us = u_arr[i];
        pars = par_arr[i];
        ax = ax1;

        for (j,u) in enumerate(us)
            p = pars[j];
            dSST, D, DFR, CT, V, qft0, Gamma_q, sft0, Gamma_s, A = p;
            zi, hM, qM, SST = u;
            zb = LCL(zi,hM,qM);
            p_avg = pres((zi+zb)/2.0,hM,qM);
            T_avg = temp((zi+zb)/2.0,hM,qM);
            rho_avg = rho(p_avg, T_avg);
            Γl = 1.5e-6;
            
            lwp = 0.5*rho_avg*Γl*(zi-zb)^2;

            if i == 1
                ax = ax1;
                p1, = ax.plot(DFR,lwp*1e3,marker=symbols[i],color=colors[i],markerfacecolor="none",markersize=10,label=L"$f_r$ (W/m$^2$)");
            end
            if i == 2
                ax = ax2;
                p2, = ax.plot(sft0,lwp*1e3,marker=symbols[i],color=colors[i],markerfacecolor="none",markersize=10,label=L"$s^{ft}_0/C_p$ (K)");
            end
            if i == 3
                ax = ax3;
                p3, = ax.plot(D*1e6,lwp*1e3,marker=symbols[i],color=colors[i],markerfacecolor="none",markersize=10,label=L"$D$ ($10^{-6}$ s$^{-1}$)");
            end
        end
        ax.tick_params(axis="x",labelcolor=colors[i]);
    end
    ax3.spines["top"].set_position(("axes", 1.15));
    ax1.set_ylabel(L"LWP (g/m$^2$)");
    ax1.grid(linestyle=":");
    ax2.set_xticks(collect(295:0.5:297));
    dots = [p3,p2,p1];
    ax1.legend(dots, [d.get_label() for d in dots], loc=4);
    ax1.text(31,10,"MLM")
    ylim([0,210]);
    tight_layout();
    savefig(string("./figs/cloudiness/",save,"_lwp_all.png"),dpi=300);
    show();
end

end