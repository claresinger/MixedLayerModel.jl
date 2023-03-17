using MixedLayerModel
using FileIO
using NCDatasets
using Plots
using Statistics
include("mlm_solve_funcs.jl")

path = "experiments/figures/20221211_perturb_from_mean/";
mkpath(path);

function allsky_lwp(u0, par)
    dt, tmaximum = 5, 100; # days
    _, sol = run_mlm_from_init(u0, par, dt=3600.0*24.0*dt, tspan=(0.0,3600.0*24.0*tmaximum), quiet=true);
    uf = sol.u[end];
    cf_ss = uf[5];
    zb_ss = calc_LCL(uf);
    ic_lwp_ss = incloud_LWP(uf, zb_ss);
    lwp = ic_lwp_ss*cf_ss;
    swcre = calc_SWCRE(uf);
    return lwp, swcre
end

# load boundary conditions from file
file = "experiments/data/box_BCs_dailystats_JJA_NEP_subonly.nc";
ds = Dataset(file, "r");

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
skip1 = 5
skip2 = 5
lon = ds["lon"][1:skip1:end]
lat = ds["lat"][1:skip2:end]
N, M = length(lon), length(lat)
println(N, " x ", M, " = ", N*M)

vars = ["sst","WS","D500","RH500","EIS","CO2"]
CF0 = fill(NaN, (N,M))
dLWPdX = fill(NaN, (N,M,length(vars)))
dCREdX = fill(NaN, (N,M,length(vars)))

for (i1,loni) in enumerate(lon)
    for (i2, lati) in enumerate(lat)
        local j1 = 1 + (i1-1)*skip1
        local j2 = 1 + (i2-1)*skip2
        # baseline
        if typeof(ds["sst_mean"][j1,j2]) == Missing
            println(loni, " ", lati)
            println(ds["sst_mean"][j1,j2])
        else
            par.SST0 = ds["sst_mean"][j1,j2] # (K)
            par.V = ds["WS_mean"][j1,j2] # m/s
            par.D = ds["D500_mean"][j1,j2] # (1/s)
            par.RHft = ds["RH500_mean"][j1,j2] # (-)
            par.EIS0 = ds["EIS_mean"][j1,j2] # (K)
            par.CO2 = 400 # (ppm)
            _, sol = run_mlm(par, dt=3600.0*24.0*5, tspan=(0.0,3600.0*24.0*100), quiet=true);
            u0 = sol.u[end];
            CF0[i1,i2] = u0[5];

            if CF0[i1,i2] > 0.5
                # loop for all vars and do + and - perturbations
                for (k,var) in enumerate(vars)
                    # decrease
                    if var == "CO2"
                        meanminus = 400 - 50;
                    else
                        meanminus = ds[var*"_mean"][j1,j2] - ds[var*"_std"][j1,j2] / 2;
                    end
                    par.SST0 = (var == "sst") ? meanminus : ds["sst_mean"][j1,j2] # (K)
                    par.V = (var == "WS") ? meanminus : ds["WS_mean"][j1,j2] # m/s
                    par.D = (var == "D500") ? meanminus : ds["D500_mean"][j1,j2] # (1/s)
                    par.RHft = (var == "RH500") ? meanminus : ds["RH500_mean"][j1,j2] # (-)
                    par.EIS0 = (var == "EIS") ? meanminus : ds["EIS_mean"][j1,j2] # (K)
                    par.CO2 = (var == "CO2") ? meanminus : 400 # (ppm)
                    lwpm, swcrem = allsky_lwp(u0, par);
                    
                    # increase
                    if var == "CO2"
                        meanplus = 400 + 50;
                    else
                        meanplus = ds[var*"_mean"][j1,j2] + ds[var*"_std"][j1,j2] / 2;
                    end
                    par.SST0 = (var == "sst") ? meanplus : ds["sst_mean"][j1,j2] # (K)
                    par.V = (var == "WS") ? meanplus : ds["WS_mean"][j1,j2] # m/s
                    par.D = (var == "D500") ? meanplus : ds["D500_mean"][j1,j2] # (1/s)
                    par.RHft = (var == "RH500") ? meanplus : ds["RH500_mean"][j1,j2] # (-)
                    par.EIS0 = (var == "EIS") ? meanplus : ds["EIS_mean"][j1,j2] # (K)
                    par.CO2 = (var == "CO2") ? meanplus : 400 # (ppm)
                    lwpp, swcrep = allsky_lwp(u0, par);

                    # save
                    dLWPdX[i1,i2,k] = lwpp - lwpm;
                    dCREdX[i1,i2,k] = swcrep - swcrem;
                end
            end
        end
    end
end
close(ds) # close data file

# create output netcdf file
ds = Dataset(path*"partial_derivatives.nc","c")

# Define the dimension "lon" and "lat" with the size N and M resp.
defDim(ds,"lon",N)
defDim(ds,"lat",M)
defDim(ds,"var",length(vars))
defVar(ds,"lon",lon,("lon",))
defVar(ds,"lat",lat,("lat",))
defVar(ds,"var",vars,("var",))

# # Define the variables
v = defVar(ds,"CF0",CF0,("lon","lat"))
v.attrib["long_name"] = "baseline cloud fraction"
v = defVar(ds,"dLWPdX",dLWPdX,("lon","lat","var"))
v.attrib["long_name"] = "ΔLWP / σX (or 100ppmv CO2)"
v = defVar(ds,"dCREdX",dCREdX,("lon","lat","var"))
v.attrib["long_name"] = "ΔCRE / σX (or 100ppmv CO2)"
close(ds)

# make preliminary plot
cmap = :bluesreds
clims = (-0.5,0.5)

p1 = contourf(lon, lat, dLWPdX[:,:,1]', c=cmap, clims=clims, 
    title="ΔLWP / σSST", xlabel="lon", ylabel="lat",
    left_margin=15Plots.mm, right_margin=10Plots.mm,
    bottom_margin=15Plots.mm, top_margin=10Plots.mm)
p2 = contourf(lon, lat, dLWPdX[:,:,2]', c=cmap, clims=clims, 
    title="ΔLWP / σWS", xlabel="lon", ylabel="lat",
    left_margin=15Plots.mm, right_margin=10Plots.mm,
    bottom_margin=15Plots.mm, top_margin=10Plots.mm)
p3 = contourf(lon, lat, dLWPdX[:,:,3]', c=cmap, clims=clims, 
    title="ΔLWP / σD500", xlabel="lon", ylabel="lat",
    left_margin=15Plots.mm, right_margin=10Plots.mm,
    bottom_margin=15Plots.mm, top_margin=10Plots.mm)
p4 = contourf(lon, lat, dLWPdX[:,:,4]', c=cmap, clims=clims, 
    title="ΔLWP / σRH500", xlabel="lon", ylabel="lat",
    left_margin=15Plots.mm, right_margin=10Plots.mm,
    bottom_margin=15Plots.mm, top_margin=10Plots.mm)
p5 = contourf(lon, lat, dLWPdX[:,:,5]', c=cmap, clims=clims, 
    title="ΔLWP / σEIS", xlabel="lon", ylabel="lat",
    left_margin=15Plots.mm, right_margin=10Plots.mm,
    bottom_margin=15Plots.mm, top_margin=10Plots.mm)
p6 = contourf(lon, lat, dLWPdX[:,:,6]', c=cmap, clims=clims, 
    title="ΔLWP / 100 ppm CO2", xlabel="lon", ylabel="lat",
    left_margin=15Plots.mm, right_margin=10Plots.mm,
    bottom_margin=15Plots.mm, top_margin=10Plots.mm)

plot(p1, p2, p3, p4, p5, p6, layout=(2,3), size=(1600,800), dpi=200)
savefig(path*"partial_derivatives_JJA_NEP.png")
