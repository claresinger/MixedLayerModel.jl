###########
# import functions and constants from other files
###########
using PyPlot
using LaTeXStrings
using Polynomials
include("EntrainmentDefs.jl")

###########
# functions
###########
function plot_lwp(ubase, pbase, us, pars, var)
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
    rcParams["font.size"] = 15;
    fig = figure(figsize=(7.5,6));
    
    len = length(us);
    x,lwp = zeros(len), zeros(len);
    
    x0 = getX(var, ubase, pbase);
    zi0, hM0, qtM0, SST0 = ubase;
    lwp0 = calc_LWP(zi0,hM0,qtM0) * 1000.0;
    
    for (i,u) in enumerate(us)
        par = pars[i];
        zi, hM, qtM, SST = u;
        lwp[i] = calc_LWP(zi,hM,qtM) * 1000.0; #kg/m^2 -> g/m^2
        x[i] = getX(var, u, par);
    end
    
    a, b = (x .- x0) ./ x0 * 100.0, (lwp .- lwp0) ./ lwp0 * 100.0;
    plot(a, b, "ko");
    fit = polyfit(a, b, 1);
    c1,c2 = coeffs(fit);
    println("slope = ", c2);
    plot(a, c1 .+ c2 .* a, "r-");
    text(-2,5,string("y=",round(c2*10)/10," x + ",round(c1*10)/10));
    
    title(getTitle(var));
    xlabel(string("% change in ",getXlabel(var)));
    ylabel("% change in LWP");
    
    ent_str = string(typeof(pbase.etype));
    rad_str = string(typeof(pbase.rtype));
    flux_str = string(typeof(pbase.ftype));
    dirname = string("./figs/linAnalysis/",ent_str,"_",rad_str,"_",flux_str,"/");
    if !isdir(dirname)
        mkpath(dirname);
    end
    savename = string(dirname,var,".png");
    savefig(savename,dpi=300);
    
    show()
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
        x = Tft - Tzi;
    elseif var == "CTq"
        x = par.CTq;
    elseif var == "LHF"
        x = par.LHF;
    elseif var == "SHF"
        x = par.SHF;
    elseif var == "SST"
        x = SST;
    end
    return x
end
                    
function getXlabel(var)
    if var == "DFR"
        name = L"$\Delta R$";
    elseif var == "D"
        name = L"$D$";
    elseif var == "inv"
        name = L"$\Delta_i T$";
    elseif var == "LHF"
        name = "LHF";
    elseif var == "SHF"
        name = "SHF";
    elseif var == "SST"
        name = "SST";
    end

    return name
end
                                            
function getTitle(var)
    if var == "DFR"
        name = "Radiative cooling";
    elseif var == "D"
        name = "Subsidence";
    elseif var == "inv"
        name = "Inversion strength";
    elseif var == "LHF"
        name = "Latent heat flux";
    elseif var == "SHF"
        name = "Sensible heat flux";
    elseif var == "SST"
        name = "Sea surface temperature";
    end

    return name
end