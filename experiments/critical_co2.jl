using MixedLayerModel
using Plots
using NCDatasets
include("mlm_solve_funcs.jl")

exp_path = "experiments/figures/20221206_critical_co2/";
mkpath(exp_path)

function co2_loop(par; co2try=400:100:2500)
    dt, tmax = 10.0, 100.0; # days

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
        if uf[5] < 0.3
            return co2i
        end
    end
    return NaN
end

# par = upCO2();
# par.etype = enBal();
# par.fttype = co2EIS();
# par.rtype = varRad();
# par.stype = fixSST();

# # adjust tunable parameters
# par.Cd = 7.9e-4;
# par.α_vent = 1.69e-3;
# par.SW_b = 140;
# # par.EIS0 = 10.0;
# # par.ECS = 1.5;
# # par.Eexport = 10.0;

# N = 11
# x_list = zeros(3,N)
# co2c_list = zeros(3,N)

# # create output netcdf file
# isfile(exp_path*"critical_co2.nc") ? rm(exp_path*"critical_co2.nc") : "no file"
# ds = Dataset(exp_path*"critical_co2.nc","c")
# defDim(ds,"var",3)
# defDim(ds,"x",N)
# defVar(ds,"var",["α_vent","Cd","SW_b"],("var",))

# # Cd
# x_list[2,:] = range(0.95,1.05,N) .* 7.9e-4
# for (i,x) in enumerate(x_list[2,:])
#     par.Cd = x
#     co2c_list[2,i] = co2_loop(par, co2try=400:50:2000)
# end
# par.Cd = 7.9e-4

# # α_vent
# x_list[1,:] = range(0.85,1.15,N) .* 1.69e-3
# for (i,x) in enumerate(x_list[1,:])
#     par.α_vent = x
#     co2c_list[1,i] = co2_loop(par, co2try=400:50:2000)
# end
# par.α_vent = 1.69e-3

# # SWb
# x_list[3,:] = range(0.75,1.25,N) .* 140
# for (i,x) in enumerate(x_list[3,:])
#     par.SW_b = x
#     co2c_list[3,i] = co2_loop(par, co2try=unique([400:100:1100; 1100:50:1800; 1800:200:2000]))
# end
# par.SW_b = 140

# defVar(ds,"Xval",x_list,("var","x"))
# defVar(ds,"critCO2",co2c_list,("var","x"))
# print(ds)
# close(ds)

ds = Dataset(exp_path*"critical_co2.nc","r")
println(ds)
x_list = ds["Xval"]
co2c_list = ds["critCO2"]
N = length(co2c_list[1,:])
println(N)

# plot
default = Int(ceil(N/2))
pCd = plot(x_list[2,:]*1e4, co2c_list[2,:], marker=:circle, color=:black, label=false, 
    xlabel="\$V\$ [mm s⁻¹]", ylabel="Critical CO₂ [ppmv]", ylims=[500,2000],
    title="a)", titleloc=:left, titlefont = font(10))
plot!([x_list[2,default]*1e4], [co2c_list[2,default]], marker=:circle, color=:red3, label=false)
pα = plot(x_list[1,:]*1e3, co2c_list[1,:], marker=:circle, color=:black, label=false,
    xlabel="\$\\alpha_{\\mathrm{vent}}\$ [mm s⁻¹]", ylims=[500,2000], ytick=([500,1000,1500,2000],[]),
    title="b)", titleloc=:left, titlefont = font(10))
plot!([x_list[1,default]*1e3], [co2c_list[1,default]], marker=:circle, color=:red3, label=false)
pSWb = plot(x_list[3,:], co2c_list[3,:], marker=:circle, color=:black, label=false,
    xlabel="\$b_{\\mathrm{SW}}\$ [W m⁻²]", ylims=[500,2000], ytick=([500,1000,1500,2000],[]),
    title="c)", titleloc=:left, titlefont = font(10))
plot!([x_list[3,default]], [co2c_list[3,default]], marker=:circle, color=:red3, label=false)

plot(pCd, pα, pSWb, layout=(1,3), size=(900,300), dpi=300, 
    left_margin=5Plots.mm, bottom_margin=5Plots.mm)
savefig(exp_path*"co2_crit_params.png")

close(ds)