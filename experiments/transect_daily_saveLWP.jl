using NCDatasets
using MixedLayerModel: calc_LCL, incloud_LWP
using MixedLayerModel: Cp, g

path = "experiments/figures/20230215_dailytransect_subonly_100days_skip1_1var/"
srcfile = "transect_output_all.nc"
newfile = "transect_output_all_LWP.nc"
cp(path*srcfile, path*newfile; force=true)  # overwrite if needed

ds = Dataset(path*newfile, "a")  # "a" = append/update

# Read dimension sizes
Ndays = length(ds["time"][:])
Nlon = length(ds["lon"][:])
Nvar = length(ds["var"][:])

# Preallocate uf_save (Float64 or whatever your original type was)
u = zeros(Float64, Ndays, Nlon, Nvar)

# Assign each variable slice back into uf_save[:, :, i]
var_names = ["zi", "sM", "qtM", "sst", "cf"]
for (i, name) in enumerate(var_names)
    u[:,:,i] = ds[name]
end

zi, sM, qtM, SST, CF = u;
zb = zeros(Float64, Ndays, Nlon)
LWP = zeros(Float64, Ndays, Nlon)
for d in 1:Ndays
    for l in 1:Nlon
        try
            zb[d,l] = calc_LCL(u[d,l,:])
            LWP[d,l] = incloud_LWP(u[d,l,:], zb[d,l])
        catch
            zb[d,l] = NaN
            LWP[d,l] = NaN
        end
    end
end

v = defVar(ds,"zb",zb,("time","lon"))
v.attrib["units"] = "m"
v.attrib["long_name"] = "cloud base altitude"

v = defVar(ds,"icLWP",LWP,("time","lon"))
v.attrib["units"] = "kg/m2"
v.attrib["long_name"] = "in-cloud liquid water path"

close(ds)