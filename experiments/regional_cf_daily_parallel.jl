using MixedLayerModel
using FileIO
using NCDatasets
using Plots
using Statistics
using StatsBase
using Random
# include("mlm_solve_funcs.jl")

using Distributed
addprocs(8; exeflags = "--project=experiments/")

include("one_day_cf_regional.jl")


path = "experiments/figures/20221219_CFdaily_parallel/";
mkpath(path);

# load boundary conditions from file
file = "experiments/data/regional_daily_good_BCs_JJA_NEP_subonly.nc";
ds = Dataset(file, "r");
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
skip1 = 200 #20, 5
skip2 = 200 #10, 5
lon = ds["lon"][1:skip1:end]
lat = ds["lat"][1:skip2:end]
N, M = length(lon), length(lat)

# import xarray as xr
# ds = xr.open_dataset("data/box_BCs_daily_JJA_NEP_subonly.nc")
# x = ds.where((ds.D500 > 0) & (ds.D500 < 10e-6) 
#                 & (ds.EIS > 0) & (ds.EIS < 20) 
#                 & (ds.RH500 < 1) 
#                 & (ds.WS > 0) & (ds.WS < 15)
#                 & (ds.sst == ds.sst))
# x.to_netcdf("data/regional_daily_good_BCs_JJA_NEP_subonly.nc")

Ndays = 10 # about 55% success rate
Random.seed!(1234)
days_indices = sample(1:length(time), Ndays, replace=false, ordered=true)
CFsave = fill(NaN, (N, M, Ndays))
println(size(CFsave))

# for (i, dayi) in enumerate(days_indices)
#     println(i, ": ", dayi)
#     CFsave[:,:,i] = one_day(dayi, ds, skip1, skip2, par)
# end

out = pmap(x -> RegionalCF.one_day(x, ds, skip1, skip2, par), days_indices)
print(out)

# CFsave[:,:,:] = hcat(out...)

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
