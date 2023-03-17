using MixedLayerModel
using MixedLayerModel: Cp, g
using FileIO
using NCDatasets
using Plots
using Statistics
using StatsBase
using Random

path = "experiments/figures/20230215_dailytransect_subonly_100days_skip1_1var/"
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
par.α_vent = 1.69e-3;
par.Cd = 6e-4; #7.9e-4;

# load boundary conditions from file
file = "experiments/data/transect_BCs_all_JJA_NEP_subonly.nc";
ds = Dataset(file, "r");
time = ds["time"]

Ndays = 100
Random.seed!(1234)
days_indices = sample(1:length(time), Ndays, replace=false, ordered=true)

skipi = 1
lon = ds["lon"][1:skipi:end]
Nlon = length(lon)
println(Nlon, " out of ", length(ds["lon"]), " longitudes")

uf_save = fill(NaN, (Ndays, Nlon, 5));
bc_vars = ["sst", "WS", "D500", "RH500", "EIS"]
Nvar = length(bc_vars)
cf_1var = fill(NaN, (Ndays, Nlon, Nvar));

# Threads.@threads for id in 1:Ndays
    # GC.safepoint(); # garbage collecting for threads
for id in 1:Ndays
    day = days_indices[id]
    println(id, " / ", Ndays, ": ", day)
    @time begin
        for (i,loni) in enumerate(lon)
            local j = (i-1)*skipi+1
            par.SST0 = ds["sst"][j, day];
            par.V = ds["WS"][j, day];
            par.D = ds["D500"][j, day];
            par.RHft = ds["RH500"][j, day];
            par.EIS0 = ds["EIS"][j, day];
            dt, tmax = 5*24*3600, 100*24*3600; # seconds
            try
                local _, sol = run_mlm(par, init=1, dt=dt, tspan=(0.0,tmax), quiet=true);
                uf_save[id, i, :] = sol.u[end];
            catch
                println("fail! ", day, ", ", loni, ", all")
            end

            for (v,var) in enumerate(bc_vars)
                par.SST0 = var=="sst" ? ds["sst"][j,day] : mean(ds["sst"][1,:]);
                par.V = var=="WS" ? ds["WS"][j,day] : mean(ds["WS"][1,:]);
                par.D = var=="D500" ? ds["D500"][j,day] : mean(ds["D500"][1,:]);
                par.RHft = var=="RH500" ? ds["RH500"][j,day] : mean(ds["RH500"][1,:]);
                par.EIS0 = var=="EIS" ? ds["EIS"][j,day] : mean(ds["EIS"][1,:]);
                try
                    local _, sol = run_mlm(par, init=1, dt=dt, tspan=(0.0,tmax), quiet=true);
                    cf_1var[id, i, v] = sol.u[end][5];
                catch
                    println("fail! ", day, ", ", loni, ", ", var)
                end
            end
        end
        println(uf_save[id, :, 5])
    end
end

cf_mean = fill(NaN, (Nlon))
for (i,loni) in enumerate(lon)
    local j = (i-1)*skipi+1
    par.SST0 = mean(ds["sst"][j, :]);
    par.V = mean(ds["WS"][j, :]);
    par.D = mean(ds["D500"][j, :]);
    par.RHft = mean(ds["RH500"][j, :]);
    par.EIS0 = mean(ds["EIS"][j, :]);
    dt, tmax = 5*24*3600, 100*24*3600; # seconds
    local _, sol = run_mlm(par, init=1, dt=dt, tspan=(0.0,tmax), quiet=true);
    cf_mean[i] = sol.u[end][5];
end

# create output netcdf file
filename = "transect_output_all.nc"
isfile(path*filename) ? rm(path*filename) : "no file"
ds_save = Dataset(path*filename,"c")

# Define the dimension "lon" of size Nindex
defDim(ds_save, "lon", Nlon)
defVar(ds_save, "lon", lon, ("lon",))

# Define the dimension "time" of size Ndays
defDim(ds_save, "time", Ndays)
defVar(ds_save, "time", ds["time"][days_indices], ("time",))

# Define the dimension "var" of size Nvar
defDim(ds_save, "var", Nvar)
defVar(ds_save, "var", bc_vars, ("var",))

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

v = defVar(ds_save,"cf_mean",cf_mean,("lon",))
v.attrib["units"] = "-"
v.attrib["long_name"] = "cloud fraction"

v = defVar(ds_save,"cf_1var",cf_1var,("time","lon","var"))
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
