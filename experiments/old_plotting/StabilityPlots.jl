###########
# import functions and constants from other files
###########
using PyPlot
using LaTeXStrings
include("EntrainmentDefs.jl")

###########
# functions
###########
function plot_S_multvar(u0, par0, us_list, pars_list, var_list)
    fig = figure(figsize=(9,6));
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
    rcParams["font.size"] = 15;

    ax1 = subplot(221);
    ax2 = subplot(222);
    ax3 = subplot(223);
    ax4 = subplot(224);
    
    color_list = ["#E2CC00","#FF6C0C","#73A950","#00A1DF"]
    for (j,var) in enumerate(var_list)
        us = us_list[j];
        pars = pars_list[j];

        x0 = getX(var, u0, par0);
        ΔL0, LHF0, h0, S0 = calc_S(u0, par0);

        for (i,u) in enumerate(us)
            xi = getX(var, u, pars[i]);
            ΔL, LHF, h, S = calc_S(u, pars[i]);

            #x = log(xi/x0);
            x = xi/x0;
            y1, y2, y3, y4 = ΔL, LHF, h, S;

            ax1.plot(x, y1, "o", color = color_list[j])
            ax2.plot(x, y2, "o", color = color_list[j])
            ax3.plot(x, y3, "o", color = color_list[j])
            if i == 1
                name,unit = getXlabel(var);
                ax4.plot(x, y4, "o", color = color_list[j], label=name)
            else
                ax4.plot(x, y4, "o", color = color_list[j])
            end
        end
    end

    # ax4.set_xscale(:log)
    # ax4.set_yscale(:log)

    #xlab = string(L"$\Delta \log$(X)")
    xlab = L"$X/X_0$"
    ax1.set_xlabel(xlab);
    ax2.set_xlabel(xlab);
    ax3.set_xlabel(xlab);
    ax4.set_xlabel(xlab);

    ax1.set_ylabel(L"$\Delta R$ (W/m$^2$)");
    ax2.set_ylabel(L"LHF (W/m$^2$)");
    ax3.set_ylabel(L"z_c/z_i");
    ax4.set_ylabel(L"instability, $S$");
        
    ax3.set_ylim([0,1]);
    ax4.set_ylim([0,1]);

    ax1.grid(linestyle=":");
    ax2.grid(linestyle=":");
    ax3.grid(linestyle=":");
    ax4.grid(linestyle=":");

    ax4.legend(loc=1,fontsize=12);
    
    figtitle, savename = get_saves(var, par0, "all_var");
    title(figtitle);
    tight_layout();
    savefig(savename,dpi=300);

    show();

end
    
function plot_S2_multvar(u0, par0, us_list, pars_list, var_list)
    fig = figure(figsize=(9,6));
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
    rcParams["font.size"] = 15;

    ax1 = subplot(221);
    ax2 = subplot(222);
    ax3 = subplot(223);
    ax4 = subplot(224);
    
    color_list = ["#E2CC00","#FF6C0C","#73A950","#00A1DF"]
    for (j,var) in enumerate(var_list)
        us = us_list[j];
        pars = pars_list[j];

        x0 = getX(var, u0, par0);
        ΔL0, LHF0, h0, S0 = calc_S(u0, par0);
        zi0, hM0, qM0, SST0 = u0;
        zb0 = LCL(zi0, hM0, qM0);
            
        len = length(us);
        x,y1,y2,y3,y4 = zeros(len),zeros(len),zeros(len),zeros(len),zeros(len);

        for (i,u) in enumerate(us)
            xi = getX(var, u, pars[i]);
            ΔL, LHF, h, S = calc_S(u, pars[i]);
            zi, hM, qM, SST = u;
            zb = LCL(zi, hM, qM);
        
            x[i] = xi/x0;
            y1[i], y2[i], y3[i], y4[i] = zi, zb, h, S;

        end
        name,unit = getXlabel(var);
        ax1.plot(x, y1, "-", color = color_list[j])
        ax2.plot(x, y2, "-", color = color_list[j])
        ax3.plot(x, y3, "-", color = color_list[j])
        ax4.plot(x, y4, "-", color = color_list[j], label=name)
    end

    xlab = L"$X/X_0$"
    ax1.set_xlabel(xlab);
    ax2.set_xlabel(xlab);
    ax3.set_xlabel(xlab);
    ax4.set_xlabel(xlab);
        
    ax1.set_ylabel(L"$z_i$ (m)");
    ax2.set_ylabel(L"$z_b$ (m)");
    ax3.set_ylabel(L"$z_c/z_i = (z_i-z_b)/z_i$");
    ax4.set_ylabel(L"$S = (LHF/\Delta R)*(z_c/z_i)$");
        
    ax3.set_ylim([0,1]);
    ax4.set_ylim([0,1]);

    ax1.grid(linestyle=":");
    ax2.grid(linestyle=":");
    ax3.grid(linestyle=":");
    ax4.grid(linestyle=":");

    ax4.legend(loc=1,fontsize=12);
    
    figtitle, savename = get_saves(var, par0, "all_var2");
    title(figtitle);
    tight_layout();
    savefig(savename,dpi=300);

    show();

end
    
function plot_S3(u0, par0, us, pars, var)
    # s fig
    fig = figure(figsize=(9,6));
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
    rcParams["font.size"] = 15;

    ax1 = subplot(221);
    ax2 = subplot(222);
    ax3 = subplot(223);
    ax4 = subplot(224);

    x0 = getX(var, u0, par0);
    ΔR0, LHF0, h0, S0 = calc_S(u0, par0);
    f0 = LHF0/ΔR0;
    zi0, hM0, qM0, SST0 = u0;
    zb0 = LCL(zi0, hM0, qM0);
        
    len = length(us);
    x,y1a,y1b,y2,y3,y4 = zeros(len),zeros(len),zeros(len),zeros(len),zeros(len),zeros(len);

    for (i,u) in enumerate(us)
        x[i] = getX(var, u, pars[i]);
        ΔR, LHF, h, S = calc_S(u, pars[i]);
        f = LHF/ΔR;
        zi, hM, qM, SST = u;
        zb = LCL(zi, hM, qM);
        
        y1a[i], y1b[i], y2[i], y3[i], y4[i] = zi, zb, h, f, S;
    end
        
    ax1.plot(x, y1a, "k-", label=L"$z_i$")
    ax1.plot(x, y1b, "k:", label=L"$z_b$")
    ax2.plot(x, y2, "k-")
    ax3.plot(x, y3, "k-")
    ax4.plot(x, y4, "k-")

    y1a, y1b, y2, y3, y4 = zi0, zb0, h0, f0, S0;
    ax1.plot(x0, y1a, "ro")
    ax1.plot(x0, y1b, "ro")
    ax2.plot(x0, y2, "ro")
    ax3.plot(x0, y3, "ro")
    ax4.plot(x0, y4, "ro")

    name, unit = getXlabel(var);
    xlab = string(name," ",unit);

    ax1.set_xlabel(xlab);
    ax2.set_xlabel(xlab);
    ax3.set_xlabel(xlab);
    ax4.set_xlabel(xlab);

    ax1.set_ylabel("z (m)");
    ax1.legend()
    ax2.set_ylabel(L"$z_c/z_i$");
    ax3.set_ylabel(L"$LHF/\Delta R$");
    ax4.set_ylabel(L"$S$");
        
    ax2.set_ylim([0,1]);
    ax4.set_ylim([0,1]);

    ax1.grid(linestyle=":");
    ax2.grid(linestyle=":");
    ax3.grid(linestyle=":");
    ax4.grid(linestyle=":");

    figtitle, savename = get_saves(var, par0, "inst3");
    title(figtitle);
    tight_layout();
    savefig(savename,dpi=300);

    close();
end

function plot_S2(u0, par0, us, pars, var)
    # s fig
    fig = figure(figsize=(9,6));
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
    rcParams["font.size"] = 15;

    ax1 = subplot(221);
    ax2 = subplot(222);
    ax3 = subplot(223);
    ax4 = subplot(224);

    x0 = getX(var, u0, par0);
    ΔL0, LHF0, h0, S0 = calc_S(u0, par0);
    zi0, hM0, qM0, SST0 = u0;
    zb0 = LCL(zi0, hM0, qM0);
        
    len = length(us);
    x,y1,y2,y3,y4 = zeros(len),zeros(len),zeros(len),zeros(len),zeros(len);

    for (i,u) in enumerate(us)
        x[i] = getX(var, u, pars[i]);
        ΔL, LHF, h, S = calc_S(u, pars[i]);
        zi, hM, qM, SST = u;
        zb = LCL(zi, hM, qM);
        
        #y1, y2, y3, y4 = zi, zb, h, S;
        y1[i], y2[i], y3[i], y4[i] = zi, zb, h, S;
    end
        
    ax1.plot(x, y1, "k-")
    ax2.plot(x, y2, "k-")
    ax3.plot(x, y3, "k-")
    ax4.plot(x, y4, "k-")

    y1, y2, y3, y4 = zi0, zb0, h0, S0;
    ax1.plot(x0, y1, "ro")
    ax2.plot(x0, y2, "ro")
    ax3.plot(x0, y3, "ro")
    ax4.plot(x0, y4, "ro")

    name, unit = getXlabel(var);
    xlab = string(name," ",unit);

    ax1.set_xlabel(xlab);
    ax2.set_xlabel(xlab);
    ax3.set_xlabel(xlab);
    ax4.set_xlabel(xlab);

    ax1.set_ylabel(L"$z_i$ (m)");
    ax2.set_ylabel(L"$z_b$ (m)");
    ax3.set_ylabel(L"$z_c/z_i = (z_i-z_b)/z_i$");
    ax4.set_ylabel(L"$S = (LHF/\Delta R)*(z_c/z_i)$");
        
    ax3.set_ylim([0,1]);
    ax4.set_ylim([0,1]);

    ax1.grid(linestyle=":");
    ax2.grid(linestyle=":");
    ax3.grid(linestyle=":");
    ax4.grid(linestyle=":");

    figtitle, savename = get_saves(var, par0, "inst2");
    title(figtitle);
    tight_layout();
    savefig(savename,dpi=300);

    close();
end

function plot_S(u0, par0, us, pars, var)
    # s fig
    fig = figure(figsize=(9,6));
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
    rcParams["font.size"] = 15;

    ax1 = subplot(221);
    ax2 = subplot(222);
    ax3 = subplot(223);
    ax4 = subplot(224);

    x0 = getX(var, u0, par0);
    ΔL0, LHF0, h0, S0 = calc_S(u0, par0);

    for (i,u) in enumerate(us)
        x = getX(var, u, pars[i]);
        ΔL, LHF, h, S = calc_S(u, pars[i]);
        
        y1, y2, y3, y4 = ΔL, LHF, h, S;
        
        ax1.plot(x, y1, "ko")
        ax2.plot(x, y2, "ko")
        ax3.plot(x, y3, "ko")
        ax4.plot(x, y4, "ko")
    end

    y1, y2, y3, y4 = ΔL0, LHF0, h0, S0;
    ax1.plot(x0, y1, "ro")
    ax2.plot(x0, y2, "ro")
    ax3.plot(x0, y3, "ro")
    ax4.plot(x0, y4, "ro")

    name, unit = getXlabel(var);
    xlab = string(name," ",unit);
    #xlab = string(L"$\Delta \log$(",name,")");

    ax1.set_xlabel(xlab);
    ax2.set_xlabel(xlab);
    ax3.set_xlabel(xlab);
    ax4.set_xlabel(xlab);

    ax1.set_ylabel(L"$\Delta R$ (W/m$^2$)");
    ax2.set_ylabel(L"LHF (W/m$^2$)");
    ax3.set_ylabel(L"z_c/z_i");
    ax4.set_ylabel(L"instability, $S$");
        
    ax3.set_ylim([0,1]);
    ax4.set_ylim([0,1]);

    ax1.grid(linestyle=":");
    ax2.grid(linestyle=":");
    ax3.grid(linestyle=":");
    ax4.grid(linestyle=":");

    figtitle, savename = get_saves(var, par0, "inst");
    title(figtitle);
    tight_layout();
    savefig(savename,dpi=300);

    close();
end

function plot_Z(u0, par0, us, pars, var)
    # z fig
    fig = figure(figsize=(9,6));
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
    rcParams["font.size"] = 15;

    ax1 = subplot(221);
    ax2 = subplot(222);
    ax3 = subplot(223);
    ax4 = subplot(224);

    x0 = getX(var, u0, par0);
    zi0, hM0, qM0, SST0 = u0;
    zb0 = LCL(zi0, hM0, qM0);

    for (i,u) in enumerate(us)
        x = getX(var, u, pars[i]);

        zi, hM, qM, SST = u;
        zb = LCL(zi, hM, qM);
        
        y1, y2, y3, y4 = zi, zb, (zi-zb), (zi-zb)/zi;
        
        ax1.plot(x, y1, "ko")
        ax2.plot(x, y2, "ko")
        ax3.plot(x, y3, "ko")
        ax4.plot(x, y4, "ko")
    end

    y1, y2, y3, y4 = zi0, zb0, (zi0-zb0), (zi0-zb0)/zi0;
    ax1.plot(x0, y1, "ro")
    ax2.plot(x0, y2, "ro")
    ax3.plot(x0, y3, "ro")
    ax4.plot(x0, y4, "ro")

    name, unit = getXlabel(var);
    xlab = string(name," ",unit);

    ax1.set_xlabel(xlab);
    ax2.set_xlabel(xlab);
    ax3.set_xlabel(xlab);
    ax4.set_xlabel(xlab);

    ax1.set_ylabel(L"$z_i$ (m)");
    ax2.set_ylabel(L"$z_b$ (m)");
    ax3.set_ylabel(L"$z_c$ (m)");
    ax4.set_ylabel(L"$z_c/z_i$");
        
    ax4.set_ylim([0,1]);

    ax1.grid(linestyle=":");
    ax2.grid(linestyle=":");
    ax3.grid(linestyle=":");
    ax4.grid(linestyle=":");

    figtitle, savename = get_saves(var, par0, "z");
    title(figtitle);
    tight_layout();
    savefig(savename,dpi=300);

    close();

end

function calc_S(u,par)
    zi, hM, qM, SST = u;
    zb = LCL(zi, hM, qM);
    ΔL = calc_DFR(zi, hM, qM, SST, par, par.rtype);
    LHF = calc_LHF(qM, SST, par);

    h = (zi - zb) / zi;
    if zi < 0.1
        h = 0.0;
        println("zi = 0")
        println(zi,"\t",zb,"\t",h)
    end
        
    S = (LHF / ΔL) * h;
    return ΔL, LHF, h, S
end

function getX(var, u, par)
    zi, hM, qM, SST = u;
    if var == "DFR"
        x = calc_DFR(zi, hM, qM, SST, par, par.rtype);
    elseif var == "D"
        x = par.D * 1e6;
    elseif var == "inv"
        Tft = temp(zi, h_ft(zi, par), q_ft(zi, par));
        Tzi = temp(zi, hM, qM);
        println(Tft,"\t",Tzi);
        x = Tft - Tzi;
    elseif var == "CTq"
        x = par.CTq;
    elseif var == "LHF"
        x = par.LHF;
    elseif var == "SST"
        x = SST;
    end

    return x
end

function getXlabel(var)
    if var == "DFR"
        name = L"$\Delta R$";
        unit = L"(W/m$^2$)";
    elseif var == "D"
        name = L"$D$";
        unit = L"($10^{-6}$ s$^{-1}$)";
    elseif var == "inv"
        name = L"$\Delta T$";
        unit = "(K)";
    elseif var == "CTq"
        name = L"$C_{Tq}$";
        unit = "";
    elseif var == "LHF"
        name = "LHF";
        unit = L"(W/m$^2$)";
    elseif var == "SST"
        name = "SST";
        unit = "(K)";
    end

    return name, unit
end
                                            
## helper function for legend titles
function get_saves(var, par, plttype)
    ent_str = string(typeof(par.etype));
    rad_str = string(typeof(par.rtype));
    flux_str = string(typeof(par.ftype));

    dirname = string("./figs/instability/",ent_str,"_",rad_str,"_",flux_str,"/");
    if !isdir(dirname)
        mkdir(dirname);
    end
    savename = string(dirname,plttype,"-",var,".png");
    figtitle = string(ent_str,", ",rad_str,", ",flux_str);

    return [figtitle,savename]
end