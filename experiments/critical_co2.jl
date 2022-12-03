using MixedLayerModel
using Plots
using NCDatasets
include("mlm_solve_funcs.jl")

exp_path = "experiments/figures/20221118_critical_co2/";
mkpath(exp_path)

function co2_loop(par; co2try=400:100:2500)
     # 400 ppm
     par.stype = fixSST();
     par.CO2 = 400;
     u0, sol = run_mlm(par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
     uf = sol.u[end];
     zb = calc_LCL(uf);
     LWP = incloud_LWP(uf, zb);
     OHU_400 = calc_OHU(uf,par,LWP,par.stype);
 
     # upsteps
     par.stype = varSST();
     for co2i in co2try
        par.CO2 = co2i;
        par.OHU = OHU_400;
        
        u0, sol = run_mlm_from_init(uf, par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
        uf = sol.u[end];
        # println(co2i, "\t", uf[5])
        if uf[5] < 0.3
            return co2i
        end
    end
    return NaN
end

par = upCO2();
par.etype = enBal();
par.fttype = co2EIS();
par.rtype = varRad();
par.stype = fixSST();
dt, tmax = 10.0, 100.0; # days

# adjust tunable parameters
par.Cd = 8e-4;
par.α_vent = 1.22e-3;
par.EIS0 = 10.0;
par.ECS = 1.5;
par.Eexport = 10.0;
par.SW_b = 150;

N = 10
x_list = zeros(3,N)
co2c_list = zeros(3,N)

# create output netcdf file
isfile(exp_path*"critical_co2.nc") ? rm(exp_path*"critical_co2.nc") : "no file"
ds = Dataset(exp_path*"critical_co2.nc","c")
defDim(ds,"var",3)
defDim(ds,"x",N)
defVar(ds,"var",["α_vent","Cd","SW_b"],("var",))

# α_vent
x_list[1,:] = range(0.9e-3,2e-3,N)
for (i,x) in enumerate(x_list[1,:])
    par.α_vent = x
    co2c_list[1,i] = co2_loop(par, co2try=400:50:2000)
end
par.α_vent = 1.22e-3

# Cd
x_list[2,:] = range(7.8e-4,9.2e-4,N)
for (i,x) in enumerate(x_list[2,:])
    par.Cd = x
    co2c_list[2,i] = co2_loop(par, co2try=400:50:2000)
end
par.Cd = 8e-4

# SWb
x_list[3,:] = range(100,200,N)
for (i,x) in enumerate(x_list[3,:])
    par.SW_b = x
    co2c_list[3,i] = co2_loop(par, co2try=unique([400:100:1100; 1100:50:1800; 1800:200:2000]))
end
par.SW_b = 150

# plot
pα = plot(x_list[1,:]*1e3, co2c_list[1,:], marker=:circle, label=false,
    xlabel="\$\\alpha_{\\mathrm{vent}} \\times 10^{3}\$ [m s⁻¹]", ylabel="Critical CO₂ [ppm]", ylims=[500,2000])
pCd = plot(x_list[2,:]*1e4, co2c_list[2,:], marker=:circle, label=false, 
    xlabel="\$C_d \\times 10^4\$", ylims=[500,2500], ytick=([500,1000,1500,2000],[]))
pSWb = plot(x_list[3,:], co2c_list[3,:], marker=:circle, label=false,
    xlabel="\$b_{\\mathrm{SW}}\$ [W m⁻²]", ylims=[500,2500], ytick=([500,1000,1500,2000],[]))
plot(pα, pCd, pSWb, layout=(1,3), size=(1000,300), dpi=300, 
    left_margin=10Plots.mm, bottom_margin=10Plots.mm)
savefig(exp_path*"co2_crit_params.png")

defVar(ds,"Xval",x_list,("var","x"))
defVar(ds,"critCO2",co2c_list,("var","x"))
print(ds)
close(ds)