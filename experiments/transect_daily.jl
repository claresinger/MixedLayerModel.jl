using MixedLayerModel
using MixedLayerModel: Cp, g
using FileIO
using NCDatasets
using Plots
using Statistics
using StatsBase
using Random

path = "experiments/figures/20221110_dailytransect/"
mkpath(path)

include("mlm_solve_funcs.jl")

# set up MLM params
par = climatology();
par.rtype = varRad();
par.stype = fixSST();
par.ftype = varFlux();
par.fttype = fixEIS();
par.etype = enBal();

par.decoup_slope = 8;
par.α_vent = 1.22e-3;
par.Cd = 3e-4;

# load boundary conditions from file
file = "experiments/data/transect_BCs_all_JJA_NEP.nc";
ds = Dataset(file, "r");
time = ds["time"]

Ndays = 100
Random.seed!(1234)
days_indices = sample(1:length(time), Ndays, replace=false, ordered=true)

skipi = 1
lon = ds["lon"][1:skipi:end]
Nlon = length(lon)
println(Nlon, " out of ", length(ds["lon"]), " longitudes")

uf_save = zeros(Ndays, Nlon, 5)
cf_monthly_mean = zeros(Nlon)

for (i,loni) in enumerate(lon)
    local j = i*skipi
    par.SST0 = ds["sst_mean"][j];
    par.V = ds["WS_mean"][j];
    par.D = max(ds["D500_mean"][j], 1e-6);
    par.RHft = ds["RH500_mean"][j];
    par.EIS0 = ds["EIS_mean"][j];
    dt, tmax = 24*5, 100; # hours, days
    u0, sol = run_mlm(par, init=1, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
    local uf = sol.u[end];
    cf_monthly_mean[i] = uf[5];
end
println(cf_monthly_mean)

for (id, day) in enumerate(days_indices)
    println(id, " / ", Ndays)
    for (i,loni) in enumerate(lon)
        local j = i*skipi
        par.SST0 = ds["sst"][j, day];
        par.V = ds["WS"][j, day];
        par.D = max(ds["D500"][j, day], 1e-6);
        par.RHft = ds["RH500"][j, day];
        par.EIS0 = ds["EIS"][j, day];
        dt, tmax = 24*5, 100; # hours, days
        u0, sol = run_mlm(par, init=1, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
        local uf = sol.u[end];
        # println(uf);
        uf_save[id, i, :] = uf;
    end
end

# create output netcdf file
isfile(path*"transect_output.nc") ? rm(path*"transect_output.nc") : "no file"
ds_save = Dataset(path*"transect_output.nc","c")

# Define the dimension "lon" of size Nindex
defDim(ds_save, "lon", Nlon)
defVar(ds_save, "lon", lon, ("lon",))

# Define the dimension "time" of size Ndays
defDim(ds_save, "time", Ndays)
defVar(ds_save, "time", ds["time"][days_indices], ("time",))

# Define the variables
v = defVar(ds_save,"zi",uf_save[:,:,1],("time","lon"))
v.attrib["units"] = "m"
v.attrib["long_name"] = "cloud top height"

v = defVar(ds_save,"sM",uf_save[:,:,2],("time","lon"))
v.attrib["units"] = "J/kg"
v.attrib["long_name"] = "liquid static energy"

v = defVar(ds_save,"qtM",uf_save[:,:,3],("time","lon"))
v.attrib["units"] = "kg/kg"
v.attrib["long_name"] = "total water specific humidity"

v = defVar(ds_save,"sst",uf_save[:,:,4],("time","lon"))
v.attrib["units"] = "K"
v.attrib["long_name"] = "sea surface temperature"

v = defVar(ds_save,"cf",uf_save[:,:,5],("time","lon"))
v.attrib["units"] = "-"
v.attrib["long_name"] = "cloud fraction"

v = defVar(ds_save,"cf_monthly_mean",cf_monthly_mean,("lon",))
v.attrib["units"] = "-"
v.attrib["long_name"] = "cloud fraction"

v = defVar(ds_save,"obs_cf_mean",ds["allsc_mean"][1:skipi:end],("lon",))
v.attrib["units"] = "-"
v.attrib["long_name"] = "observed cloud fraction mean (CASCCAD)"

v = defVar(ds_save,"obs_cf_std",ds["allsc_std"][1:skipi:end],("lon",))
v.attrib["units"] = "-"
v.attrib["long_name"] = "observed cloud fraction std (CASCCAD)"

close(ds)
close(ds_save)
