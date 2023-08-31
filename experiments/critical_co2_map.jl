using MixedLayerModel
using Plots
using NCDatasets
include("mlm_solve_funcs.jl")

exp_path = "experiments/figures/20230830_critical_co2_map/";
mkpath(exp_path)

function co2_loop(par; co2try=200:100:2000)
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

par = upCO2();
par.etype = enBal();
par.fttype = co2EIS();
par.rtype = varRad();

# adjust tunable parameters
par.Cd = 7.9e-4;
par.α_vent = 1.69e-3;
par.SW_b = 140;

file = "experiments/data/crit_co2_map_input.nc";
ds = Dataset(file, "r");
SST_list = ds["sst"][:];
EIS_list = ds["EIS"][:];
lon = ds["lon"][:]
lat = ds["lat"][:]
N = length(SST_list);
co2c_list = fill(NaN, N);

for (i,SSTi) in enumerate(SST_list)
    if (~ismissing(SSTi))
        print(i," (", SSTi, ", ", EIS_list[i],") : ")
        par.SST0 = SSTi
        par.EIS0 = 8.0 #EIS_list[i]
        co2c_list[i] = co2_loop(par, co2try=200:100:2000)
        println(co2c_list[i])
    end
end

co2c_list = reshape(co2c_list, size(ds["sst"]))
close(ds)

isfile(exp_path*"critical_co2_map.nc") ? rm(exp_path*"critical_co2_map.nc") : "no file"
ds = Dataset(exp_path*"critical_co2_map.nc","c")
nx,ny = size(co2c_list)
defDim(ds,"lon",nx)
defDim(ds,"lat",ny)
defVar(ds,"lon",lon,("lon",))
defVar(ds,"lat",lat,("lat",))
defVar(ds,"critCO2",co2c_list,("lon","lat"))
close(ds)

# N = 15
# SST_list = range(288, 294, N)
# EIS_list = range(9.5, 5, N)
# co2c_list = zeros(N)

# # create output netcdf file
# isfile(exp_path*"critical_co2.nc") ? rm(exp_path*"critical_co2.nc") : "no file"
# ds = Dataset(exp_path*"critical_co2.nc","c")
# defDim(ds,"x",N)

# for (i,SSTi) in enumerate(SST_list)
#     par.SST0 = SSTi
#     par.EIS0 = EIS_list[i]
#     co2c_list[i] = co2_loop(par, co2try=100:50:2400)
#     println(SSTi, "\t", co2c_list[i])
# end

# defVar(ds,"SST",SST_list,("x",))
# defVar(ds,"EIS",EIS_list,("x",))
# defVar(ds,"critCO2",co2c_list,("x",))
# print(ds)
# close(ds)

# ds = Dataset(exp_path*"critical_co2.nc","r")
# SST_list = ds["SST"]
# EIS_list = ds["EIS"]
# co2c_list = ds["critCO2"]
# println(SST_list[:])
# println(EIS_list[:])
# println(co2c_list[:])
# N = length(co2c_list)

# # plot
# isfile(exp_path*"co2_crit_ICs.pdf") ? rm(exp_path*"co2_crit_ICs.pdf") : "no file"
# default = 5
# ms = 6
# p = plot(SST_list, co2c_list, 
#     marker=:circle, markersize=ms, color=:black, markerstrokewidth=0, label=false, 
#     xlabel="Initial SST [K]", ylabel="Critical CO₂ [ppmv]", ylims=[100,2200],
#     size=(500,400), dpi=300, left_margin=5Plots.mm, bottom_margin=5Plots.mm, right_margin=2Plots.mm)

# plot!([SST_list[default]], [co2c_list[default]], 
#     marker=:circle, markersize=ms, color=:red3, markerstrokewidth=0, label=false,
#     z_order=:front)

# plot!(twiny(), EIS_list, co2c_list, 
#     marker=:circle, markersize=ms, color=:black, markerstrokewidth=0, label=false,
#     xflip=true, xlabel="Intitial EIS [K]", ylims=[100,2200],
#     z_order=:back)

# close(ds)
# savefig(p, exp_path*"co2_crit_ICs.pdf")
