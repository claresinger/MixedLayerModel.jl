###########
# import functions and constants from other files
###########
using PyPlot
using LaTeXStrings
include("EntrainmentDefs.jl")

###########
# functions
###########

stepz = 0.1;
maxz = 1500.0;
colors = ["red","orange","lime","green","blue","purple"];
colors = ["red","lime","purple"];
fontsize = 20;

#####
# make rainbow profile plots of generic variables "phi1" and "phi2"
# for several values of climate variable "var"
#####
function plot_three_profs_generic(us, pars, phi1, phi2, phi3, var)
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
    rcParams["font.size"] = fontsize;

    fig = figure(figsize=(15,8));
    ax1 = subplot(131);
    ax2 = subplot(132);
    ax3 = subplot(133);

    for i in 1:length(us)
        u = us[i];
        p = pars[i];
        zi, hM, qM, SST = u;

        z = [collect(0:stepz:zi);collect(zi:stepz:maxz)];

        x1 = get_phi(u,p,phi1,stepz,maxz);
        x2 = get_phi(u,p,phi2,stepz,maxz);
        x3 = get_phi(u,p,phi3,stepz,maxz);
        legname = get_leg(u,p,var);

        ax1.plot(x1,z,linewidth=2,linestyle="-",color=colors[i])
        ax2.plot(x2,z,linewidth=2,linestyle="-",color=colors[i])
        ax3.plot(x3,z,linewidth=2,linestyle="-",color=colors[i],label=legname)
    end

    phi123 = string(phi1,"_AND_",phi2,"_AND_",phi3)
    figtitle, legtitle, savename = get_saves_prof(var,phi123,pars[1]);

    ax1.set_xlabel(get_xlabel(phi1)); ax1.set_ylabel(L"$z$ (m)");
    ax1.set_xlim(get_xlim(phi1));
    ax1.set_ylim([0,maxz]);
    ax1.grid(linestyle=":");

    ax2.set_xlabel(get_xlabel(phi2));
    ax2.set_ylim([0,maxz])
    ax2.set_yticklabels([]);
    ax2.grid(linestyle=":");
    
    ax3.set_xlabel(get_xlabel(phi3));
    ax3.set_ylim([0,maxz])
    ax3.set_yticklabels([]);
    ax3.legend(title=legtitle,loc=4);
    ax3.grid(linestyle=":");

    suptitle(figtitle);
    tight_layout();
    subplots_adjust(top=0.93);
    savefig(savename,dpi=300);
    #close();
end

#####
# make rainbow profile plots of generic variables "phi1" and "phi2"
# for several values of climate variable "var"
#####
function plot_two_profs_generic(us, pars, phi1, phi2, var)
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
    rcParams["font.size"] = fontsize;

    fig1 = figure(figsize=(10,8));
    ax1 = subplot(121);
    ax2 = subplot(122);

    for i in 1:length(us)
        u = us[i];
        p = pars[i];
        zi, hM, qM, SST = u;

        z = [collect(0:stepz:zi);collect(zi:stepz:maxz)];

        x1 = get_phi(u,p,phi1,stepz,maxz);
        x2 = get_phi(u,p,phi2,stepz,maxz);
        legname = get_leg(u,p,var);

        ax1.plot(x1,z,linewidth=2,linestyle="-",color=colors[i])
        ax2.plot(x2,z,linewidth=2,linestyle="-",color=colors[i],label=legname)
    end

    phi12 = string(phi1,"_AND_",phi2)
    figtitle, legtitle, savename = get_saves_prof(var,phi12,pars[1]);

    ax1.set_xlabel(get_xlabel(phi1)); ax1.set_ylabel(L"$z$ (m)");
    ax1.grid(linestyle=":");

    ax2.set_xlabel(get_xlabel(phi2));
    ax2.set_yticklabels([]);
    ax2.legend(title=legtitle,loc=1);
    ax2.grid(linestyle=":");

    suptitle(figtitle);
    tight_layout();
    subplots_adjust(top=0.93);
    savefig(savename,dpi=300);
    close();
end

#####
# make rainbow profile plots of generic variable "phi"
# for several values of climate variable "var"
#####
function plot_prof_generic(us, pars, phi, var)
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
    rcParams["font.size"] = fontsize;

    fig1 = figure(figsize=(6,8));
    for i in 1:length(us)
        u = us[i];
        p = pars[i];
        zi, hM, qM, SST = u;

        z = [collect(0:stepz:zi);collect(zi:stepz:maxz)];

        x = get_phi(u, p, phi, stepz, maxz);
        legname = get_leg(u,p,var);
        
        plot(x,z,linewidth=2,linestyle="-",color=colors[i],label=legname)
    end

    figtitle, legtitle, savename = get_saves_prof(var,phi,pars[1]);

    xlabel(get_xlabel(phi)); ylabel(L"$z$ (m)");
    legend(title=legtitle,loc=4);
    grid(linestyle=":");
#     title(figtitle);
    tight_layout();
    savefig(savename,dpi=300);
    close();
end

## helper function compute phi
function get_phi(u, p, phi, stepz, maxz)
    zi, hM, qM, SST = u;
    DFR = calc_DFR(zi,hM,qM,SST,p,p.rtype);

    z1 = collect(0:stepz:zi); z2 = collect(zi:stepz:maxz);
    z = [z1;z2];
    h = [hM .* ones(length(z1)); h_ft(z2, p)];
    qt = [qM .* ones(length(z1)); q_ft(z2, p)];
    T = temp.(z,h,qt);
    
    if phi == "theta"       
        x = theta(z,h,qt);
    elseif phi == "T"
        x = T;
    elseif phi == "p"
        x = pres(z,T) / 100.0; # Pa -> hPa
    elseif phi == "h"        
        x = h / 1000.0; # J/kg -> kJ/kg
    elseif phi == "qt"
        x = qt * 1000.0; # kg/kg -> g/kg
    elseif phi == "RH"
        x = RH(z, h, qt) * 100.0; # -> %
    elseif phi == "ql"
        x = q_l(z,T,qt) * 1000.0; # kg/kg -> g/kg
    elseif phi == "cf"
        x = ceil.(q_l(z,T,qt)) * 100.0; # -> %
    elseif phi == "bflux"                                                                           
        zb = LCL(zi,hM,qM);
        z = collect(0:stepz:maxz);
        z1 = z[z .< zb];
        z2 = intersect(z[z .>= zb], z[z .< zi];)
        z3 = z[z .>= zi];
        
        bflux = calc_bflux(u,p,z1,z2,z3,p.etype);                            
        x = bflux * 1e4; # m^2/s^3 -> 10^-4 m^2/s^3
    else
        print("error: no option for phi == "+phi)
    end

    return x
end

## helper function for x labels
function get_xlabel(phi)
    if phi == "theta"
        str = L"$\theta$ (K)"
    elseif phi == "T"
        str = L"$T$ (K)"
    elseif phi == "p"
        str = L"$p$ (hPa)"
    elseif phi == "h"
        str = L"$h$ (kJ/kg)"
    elseif phi == "qt"
        str = L"$q_t$ (g/kg)"
    elseif phi == "RH"
        str = "Relative humidity (%)"
    elseif phi == "ql"
        str = L"$q_\ell$ (g/kg)"
    elseif phi == "cf"
        str = "Cloud fraction (%)"
    elseif phi == "bflux"
        str = L"Buoyancy flux, $\langle w' b' \rangle$ ($10^{-4}$ m$^2$/s$^3$)"
    else
        print("error: no option for phi == "+phi)
    end

    return str
end
                                                                        
## helper function for x limits for plots
function get_xlim(phi)
    if phi == "theta"
        lims = [290,310]
    elseif phi == "T"
        lims = [280,300]
    elseif phi == "p"
        lims = [850,1050]
    elseif phi == "h"
        lims = [305,321]
    elseif phi == "qt"
        lims = [1,12]
    elseif phi == "RH"
        lims = [0,100]
    elseif phi == "ql"
        lims = [0,0.7]
    elseif phi == "cf"
        lims = [0,100]
    elseif phi == "bflux"
        lims = [0,20]
    else
        print("error: no option for phi == "+phi)
    end

    return lims
end

## helper function for legend labels
function get_leg(u,p,var)
    zi, hM, qM, SST = u;
    if var == "DFR"
        legname = string(round(p.DFR*100)/100);
    elseif var == "D"
        legname = string(round(p.D*1e6*100)/100);
    elseif var == "inv"
        legname = string(round(p.sft0*100)/100);
    elseif var == "CTq"
        legname = string(round(p.CTq*1e4*100)/100);
    elseif var == "LHF"
        legname = string(round(p.LHF));
    elseif var == "SST"
        legname = string(round(SST*100)/100);
    elseif var == "SHF"
        legname = string(round(p.SHF*10)/10);
    end

    return legname
end

## helper function for legend titles
function get_saves_prof(var, phi, par)
    ent_str = string(typeof(par.etype));
    rad_str = string(typeof(par.rtype));
    flux_str = string(typeof(par.ftype));

    if var == "DFR"
        if rad_str == "direct"
            legtitle = L"$\Delta R$ (W/m$^2$)";
        elseif rad_str == "LWPcalc"
            legtitle = L"$f_r$ (W/m$^2$)";
        end
    elseif var == "D"
        legtitle = L"$D$ ($10^{-6}$ s$^{-1}$)";
    elseif var == "inv"
        legtitle = L"$\tilde{s}^{ft}_0$ (K)";
    elseif var == "CTq"
        legtitle = L"$C_{Tq}$ ($10^{-4}$)";
    elseif var == "LHF"
        legtitle = L"LHF (W/m$^2$)";
    elseif var == "SHF"
        legtitle = L"SHF (W/m$^2$)";
    elseif var == "SST"
        legtitle = "SST (K)";
    end

    dirname = string("./figs/profiles/",ent_str,"_",rad_str,"_",flux_str,"/");
    if !isdir(dirname)
        mkdir(dirname);
    end
    savename = string(dirname,var,"_",phi,".png");
    figtitle = string(ent_str,", ",rad_str,", ",flux_str);

    return [figtitle, legtitle, savename]
end