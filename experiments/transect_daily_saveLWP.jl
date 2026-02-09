using NCDatasets
using MixedLayerModel

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
dR = zeros(Float64, Ndays, Nlon)
LHF = zeros(Float64, Ndays, Nlon)
De = zeros(Float64, Ndays, Nlon)
for d in 1:Ndays
    for l in 1:Nlon
        try
            zb[d,l] = calc_LCL(u[d,l,:])
            LWP[d,l] = incloud_LWP(u[d,l,:], zb[d,l])
            dR[d,l] = calc_cloudtop_RAD(u[d,l,:], par, LWP[d,l], par.rtype);
            LHF[d,l] = calc_LHF(u[d,l,:], par)
            De[d,l] = calc_decoupling(u[d,l,:], par, zb[d,l], LWP[d,l])
        catch
            zb[d,l] = NaN
            LWP[d,l] = NaN
            dR[d,l] = NaN
            LHF[d,l] = NaN
            De[d,l] = NaN
        end
    end
end

v = defVar(ds,"zb",zb,("time","lon"))
v.attrib["units"] = "m"
v.attrib["long_name"] = "cloud base altitude"

v = defVar(ds,"icLWP",LWP,("time","lon"))
v.attrib["units"] = "kg/m2"
v.attrib["long_name"] = "in-cloud liquid water path"

v = defVar(ds,"dR",dR,("time","lon"))
v.attrib["units"] = "W/m2"
v.attrib["long_name"] = "cloud-top radiative cooling"

v = defVar(ds,"LHF",LHF,("time","lon"))
v.attrib["units"] = "W/m2"
v.attrib["long_name"] = "surface latent heat flux"

v = defVar(ds,"De",De,("time","lon"))
v.attrib["units"] = "-"
v.attrib["long_name"] = "decoupling parameter"

close(ds)