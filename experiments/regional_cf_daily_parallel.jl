using MixedLayerModel
using FileIO
using NCDatasets
using Plots
using Statistics
using StatsBase
using Random
include("mlm_solve_funcs.jl")

path = "experiments/figures/20230126_CFdaily_parallel_10days/";
mkpath(path);

# load boundary conditions from file
# import xarray as xr
# ds = xr.open_dataset("data/box_BCs_daily_JJA_NEP_subonly.nc")
# x = ds.where((ds.D500 > 0) & (ds.D500 < 10e-6) 
#                 & (ds.EIS > 0) & (ds.EIS < 20) 
#                 & (ds.RH500 < 1) 
#                 & (ds.WS > 0) & (ds.WS < 15)
#                 & (ds.sst == ds.sst))
# x.to_netcdf("data/regional_daily_good_BCs_JJA_NEP_subonly.nc")
file = "experiments/data/regional_daily_good_BCs_JJA_NEP_subonly.nc";
ds = Dataset(file, "r");
ds = NCDatasets.@select(ds, 10 <= lat <= 40 && -160 <= lon <= -110)
time = ds["time"]

# set up MLM params
par = climatology();
par.rtype = varRad();
par.stype = fixSST();
par.ftype = varFlux();
par.fttype = fixEIS();
par.etype = enBal();

par.decoup_slope = 8;
par.α_vent = 1.69e-3;
par.Cd = 6e-4; #7.9e-4;

#####
skip1 = 50 #20, 5, 200
skip2 = 50 #10, 5, 200
lon = ds["lon"][1:skip1:end]
lat = ds["lat"][1:skip2:end]
N, M = length(lon), length(lat)

Ndays = 10 # about 55% success rate
Random.seed!(1234)
days_indices = sample(1:length(time), Ndays, replace=false, ordered=true)
CFsave = fill(NaN, (N, M, Ndays))
println(size(CFsave))

@time begin
Threads.@threads for i in 1:Ndays
    dayi = days_indices[i]
    println(i, ": ", dayi)

    for i1 in 1:N
        for i2 in 1:M
            # println(i1,", ", i2)
            local j1 = 1 + (i1-1)*skip1
            local j2 = 1 + (i2-1)*skip2
            
            if typeof(ds["sst"][j1,j2,dayi]) == Missing
                println("skip")
            else
                par.SST0 = ds["sst"][j1,j2,dayi] # (K)
                par.V = ds["WS"][j1,j2,dayi] # m/s
                par.D = ds["D500"][j1,j2,dayi] # (1/s)
                par.RHft = ds["RH500"][j1,j2,dayi] # (-)
                par.EIS0 = ds["EIS"][j1,j2,dayi] # (K)
                par.CO2 = 400 # (ppm)
                try
                    _, sol = run_mlm(par, dt=3600.0*24.0*5, tspan=(0.0,3600.0*24.0*100), quiet=true);
                    CFsave[i1,i2,i] = sol.u[end][5];
                catch
                    println("fail")
                end
            end
        end
    end
end
end

# create output netcdf file
filename = path*"CF_daily.nc"
isfile(filename) ? rm(filename) : "no file"
ds_save = Dataset(filename,"c")

# Define the dimension "lon" and "lat" with the size N and M resp.
defDim(ds_save,"lon",N)
defDim(ds_save,"lat",M)
defDim(ds_save,"day",Ndays)
defVar(ds_save,"lon",lon,("lon",))
defVar(ds_save,"lat",lat,("lat",))
defVar(ds_save,"day",time[days_indices],("day",))

# Define the variables
v = defVar(ds_save,"CF",CFsave,("lon","lat","day"))
v.attrib["long_name"] = "baseline cloud fraction"

close(ds) # close data file
close(ds_save) # close output file
