###########
# import functions and constants from other files
###########
using PyPlot
using LaTeXStrings
include("EntrainmentDefs.jl")

###########
# functions
###########
function plot_generic_contour(us, pars, var1, var2, phi, exp_type)
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
    rcParams["font.size"] = 15;
    fig = figure(figsize=(7.5,6));
    
    n2 = length(us);
    n = convert(Int,round(sqrt(n2)));
    
    x, y, z = zeros(n2), zeros(n2), zeros(n2);
    for (i,u) in enumerate(us)
        x[i] = getX(var1, u, pars[i]);
        y[i] = getX(var2, u, pars[i]);
        z[i] = getZ(phi, u);
    end
    xgrid = reshape(x,(n,n));
    ygrid = reshape(y,(n,n));
    zgrid = reshape(z,(n,n));
    
    contourf(xgrid, ygrid, zgrid, levels=get_Zlevs(phi), cmap=get_Zcols(phi), extend="both");
    colorbar();

#     if var1 == "SHF" && var2 == "LHF"
#         # plot constant bowen ratio lines
#         bgrid = xgrid ./ ygrid;
#         cs = contour(xgrid,ygrid,bgrid, levels=[0.01,0.05,0.1], colors="k", linestyles=":");
#         clabel(cs,fmt="%1.2f");
#     end
    
    # plot center cross
    xcent,ycent = (maximum(x)+minimum(x))/2.0, (maximum(y)+minimum(y))/2.0;
    plot([xcent,xcent],[minimum(y),maximum(y)],"k:");
    plot([minimum(x),maximum(x)],[ycent,ycent],"k:");

    xlabel(getXlabel(var1));
    ylabel(getXlabel(var2));
    title(getZlabel(phi));
    
    savename = get_saves(var1,var2,pars[1],phi,exp_type);
    tight_layout();
    savefig(savename,dpi=300);

    show()
end

function getZ(phi, u)
    zi, hM, qtM, SST = u;
    if phi == "LWP"
        z = calc_LWP(zi,hM,qtM) * 1000.0; #kg/m^2 -> g/m^2
    elseif phi == "zi"
        z = zi;
    elseif phi == "zb"
        z = LCL(zi,hM,qtM);
    elseif phi == "zc"
        z = zi - LCL(zi,hM,qtM);            
    end
    return z
end
                    
function getZlabel(phi)
    if phi == "LWP"
        name = "LWP";
        unit = L"(g/m$^2$)";
    elseif phi == "zi"
        name = L"$z_{top}$";
        unit = "(m)";
    elseif phi == "zb"
        name = L"$z_{base}$";
        unit = "(m)";
    elseif phi == "zc"
        name = L"$z_c$";
        unit = "(m)";
    end

    return string(name," ",unit)
end
                        
function get_Zlevs(phi)
    n = 10
    if phi == "LWP"
        levs = collect(range(50, length=n, stop=150));
    elseif phi == "zi"
        levs = collect(range(500, length=n, stop=1500));
    elseif phi == "zb"
        levs = collect(range(300, length=n, stop=1300));
    elseif phi == "zc"
        levs = collect(range(0, length=n, stop=400));
    end

    return levs
end
                                    
function get_Zcols(phi)
    if phi == "LWP"
        cmap = "YlGnBu";
    elseif phi == "zi"
        cmap = "YlOrBr";
    elseif phi == "zb"
        cmap = "YlOrBr";
    elseif phi == "zc"
        cmap = "YlOrBr";
    end

    return cmap
end

function getX(var, u, par)
    zi, hM, qM, SST = u;
    if var == "DFR"
        x = calc_DFR(zi, hM, qM, SST, par, par.rtype);
    elseif var == "D"
        x = par.D * 1e6;
    elseif var == "inv"
#         Tft = temp(zi, h_ft(zi, par), q_ft(zi, par));
#         Tzi = temp(zi, hM, qM);
#         x = Tft - Tzi;
#         x = theta(zi, h_ft(zi, par), q_ft(zi, par)) - theta(zi,hM,qM);
        x = par.sft0;
    elseif var == "CTq"
        x = par.CTq;
    elseif var == "LHF"
        x = par.LHF;
    elseif var == "SHF"
        x = par.SHF;
    elseif var == "Bo"
        x = par.SHF / par.LHF;
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
        name = L"$\tilde{s}^{ft}_0$";
#         name = L"$\Delta T$";
#         name = L"$\Delta \theta$";
        unit = "(K)";
    elseif var == "CTq"
        name = L"$C_{Tq}$";
        unit = "";
    elseif var == "LHF"
        name = "LHF";
        unit = L"(W/m$^2$)";
    elseif var == "SHF"
        name = "SHF";
        unit = L"(W/m$^2$)";
    elseif var == "Bo"
        name = L"$\beta$";
        unit = "";
    elseif var == "SST"
        name = "SST";
        unit = "(K)";
    end

    return string(name," ",unit)
end
                                                                
## helper function for legend titles
function get_saves(var1, var2, par, phi, exp_type)
    ent_str = string(typeof(par.etype));
    rad_str = string(typeof(par.rtype));
    flux_str = string(typeof(par.ftype));

    dirname = string("./figs/contours/",ent_str,"_",rad_str,"_",flux_str,"/",exp_type,"/");
    if !isdir(dirname)
        mkpath(dirname);
    end
    savename = string(dirname,var1,"vs",var2,"-",phi,".png");
    #figtitle = string(ent_str,", ",rad_str,", ",flux_str);

    return savename
end