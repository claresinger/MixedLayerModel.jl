using MixedLayerModel
using Plots
using NCDatasets
include("mlm_solve_funcs.jl")

exp_path = "experiments/figures/20250407_critical_co2_BCs/";
mkpath(exp_path)

# rate = rate of change in D per warming [%/K]
function co2_loop_modD(par, rate; co2try=400:200:3000)
    dt, tmax = 10.0, 100.0; # days

    # 400 ppm
    par.stype = fixSST();
    par.CO2 = 400;
    u0, sol = run_mlm(par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
    uf = sol.u[end];
    zb = calc_LCL(uf);
    LWP = incloud_LWP(uf, zb);
    OHU_400 = calc_OHU(uf,par,LWP,par.stype);
    D_400 = par.D;

    # upsteps
    par.stype = varSST();
    for co2i in co2try
        par.CO2 = co2i;
        par.OHU = OHU_400;
        par.D = D_400 * (1 + rate * par.ECS * log(par.CO2 / 400) / log(2));
        
        u0, sol = run_mlm_from_init(uf, par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
        uf = sol.u[end];
        if uf[5] < 0.3
            par.D = D_400;
            return co2i
        end
    end
    par.D = D_400;
    return NaN
end

# rate = rate of change in RH+ per warming [%/K]
function co2_loop_modRH(par, rate; co2try=400:200:3000)
    dt, tmax = 10.0, 100.0; # days

    # 400 ppm
    par.stype = fixSST();
    par.CO2 = 400;
    u0, sol = run_mlm(par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
    uf = sol.u[end];
    zb = calc_LCL(uf);
    LWP = incloud_LWP(uf, zb);
    OHU_400 = calc_OHU(uf,par,LWP,par.stype);
    RHft_400 = par.RHft;

    # upsteps
    par.stype = varSST();
    for co2i in co2try
        par.CO2 = co2i;
        par.OHU = OHU_400;
        par.RHft = RHft_400 * (1 + rate * par.ECS * log(par.CO2 / 400) / log(2));
        
        u0, sol = run_mlm_from_init(uf, par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
        uf = sol.u[end];
        if uf[5] < 0.3
            par.RHft = RHft_400;
            return co2i
        end
    end
    par.RHft = RHft_400;
    return NaN
end

# rate = rate of change in OHU per warming [%/K]
function co2_loop_modOHU(par, rate; co2try=400:200:3000)
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
        par.OHU = OHU_400 * (1 + rate * par.ECS * log(par.CO2 / 400) / log(2));
        
        u0, sol = run_mlm_from_init(uf, par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
        uf = sol.u[end];
        if uf[5] < 0.3
            par.OHU = OHU_400;
            return co2i
        end
    end
    par.OHU = OHU_400;
    return NaN
end

par = upCO2();
par.etype = enBal();
par.fttype = co2EIS();
par.rtype = varRad();

# adjust tunable parameters
par.Cd = 7.9e-4;
par.α_vent = 1.69e-3;
par.SW_b = 140;

N = 11
x_list = zeros(3,N)
co2c_list = zeros(3,N)

# create output netcdf file
isfile(exp_path*"critical_co2.nc") ? rm(exp_path*"critical_co2.nc") : "no file"
ds = Dataset(exp_path*"critical_co2.nc","c")
defDim(ds,"var",3)
defDim(ds,"x",N)
defVar(ds,"var",["D","RHft","OHU"],("var",))

# D, subsidence
x_list[1,:] = range(-0.05, 0, N) # from 0 to 5%/K decrease
for (i,x) in enumerate(x_list[1,:])
    co2c_list[1,i] = co2_loop_modD(par, x, co2try=500:50:2000)
end
println(co2c_list)

# RHft, free-tropospheric relative humidity
x_list[2,:] = range(-0.01, 0, N) # from 0 to 1%/K decrease
for (i,x) in enumerate(x_list[2,:])
    co2c_list[2,i] = co2_loop_modRH(par, x, co2try=500:50:2000)
end
println(co2c_list)

# OHU
x_list[3,:] = range(-0.2, 0.2, N) # from -20% to +20%/K change in OHU
for (i,x) in enumerate(x_list[3,:])
    co2c_list[3,i] = co2_loop_modOHU(par, x, co2try=500:50:2000)
end

defVar(ds,"Xval",x_list,("var","x"))
defVar(ds,"critCO2",co2c_list,("var","x"))
print(ds)
close(ds)

# read output to plot
ds = Dataset(exp_path*"critical_co2.nc","r")
# println(ds)
x_list = ds["Xval"]
co2c_list = ds["critCO2"]
println(co2c_list[:])
N = length(co2c_list[1,:])
# println(N)

# plot
isfile(exp_path*"co2_crit_BCs.png") ? rm(exp_path*"co2_crit_BCs.png") : "no file"
center = Int(ceil(N/2))
ms = 6
pD = plot(x_list[1,:]*100, co2c_list[1,:], 
    marker=:circle, markersize=ms, color=:black, markerstrokewidth=0, label=false, 
    xlabel="\$\\Delta D\$ [% K⁻¹]", ylabel="Critical CO₂ [ppmv]", ylims=[500,2000],
    title="a)", titleloc=:left, titlefont = font(10))
plot!([x_list[1,end]*100], [co2c_list[1,end]], 
    marker=:circle, markersize=ms, color=:red3, markerstrokewidth=0, label=false)

pRHft = plot(x_list[2,:]*100, co2c_list[2,:], 
    marker=:circle, markersize=ms, color=:black, markerstrokewidth=0, label=false,
    xlabel="\$\\Delta\$ RH₊ [% K⁻¹]", ylims=[500,2000], ytick=([500,1000,1500,2000],[]),
    title="b)", titleloc=:left, titlefont = font(10))
plot!([x_list[2,end]*100], [co2c_list[2,end]], 
    marker=:circle, markersize=ms, color=:red3, markerstrokewidth=0, label=false)

pOHU = plot(x_list[3,:]*100, co2c_list[3,:], 
    marker=:circle, markersize=ms, color=:black, markerstrokewidth=0, label=false,
    xlabel="\$\\Delta\$ OHU [% K⁻¹]", ylims=[500,2000], ytick=([500,1000,1500,2000],[]),
    title="c)", titleloc=:left, titlefont = font(10))
plot!([x_list[3,center]*100], [co2c_list[3,center]], 
    marker=:circle, markersize=ms, color=:red3, markerstrokewidth=0, label=false)

plot(pD, pRHft, pOHU, layout=(1,3), size=(900,300), dpi=300, 
    left_margin=5Plots.mm, bottom_margin=5Plots.mm)
savefig(exp_path*"co2_crit_BCs.png")

close(ds)