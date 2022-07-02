push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using MixedLayerModel
using FileIO
using NCDatasets
using Plots
using Statistics

include("mlm_solve_funcs.jl")

function allsky_lwp(par)
    dt, tmaximum = 24, 50;
    u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmaximum));
    uf = sol.u[end];
    cf_ss = uf[5];
    zb_ss = calc_LCL(uf);
    ic_lwp_ss = incloud_LWP(uf, zb_ss);
    lwp = ic_lwp_ss*cf_ss;
    return lwp
end

# load boundary conditions from file
file = "experiments/data/box_BCs_JJA_NEP.nc";
ds = Dataset(file, "r");
# println(ds)

par = climatology();
par.rtype = varRad();
par.stype = fixSST();
par.ftype = varFlux();
par.fttype = fixEIS();
par.etype = enBal();

skip1 = 1
skip2 = 1
lon = ds["lon"][1:skip1:end]
lat = ds["lat"][1:skip2:end]
N, M = length(lon), length(lat)
println(N, " ", M)
ddSST = zeros(N, M)
ddV = zeros(N, M)
ddD = zeros(N, M)
ddRH = zeros(N, M)
ddEIS = zeros(N, M)
ddCO2 = zeros(N, M)

for (i1,loni) in enumerate(lon)
    for (i2, lati) in enumerate(lat)
        local j1 = i1*skip1
        local j2 = i2*skip2
        if i1 > 0
            println(i1, ", ", i2)
            println(j1, ", ", j2)

            par.SST0 = ds["sst"][j1,j2];
            par.V = ds["WS"][j1,j2];
            # par.LHF = ds["LHF"][j1,j2];
            # par.SHF = ds["SHF"][j1,j2];
            par.D = ds["D500"][j1,j2];
            par.RHft = ds["RH500"][j1,j2];
            par.EIS = 10.0;

            # base state
            lwp0 = allsky_lwp(par);

            # sst
            stdSST = ds["sst_std"][j1,j2];
            par.SST0 = ds["sst"][j1,j2] - stdSST;
            lwpm = allsky_lwp(par);
            par.SST0 = ds["sst"][j1,j2] + stdSST;
            lwpp = allsky_lwp(par);
            ddSST[i1,i2] = mean(diff([lwpm, lwp0, lwpp]));

            # V
            stdV = ds["WS_std"][j1,j2];
            par.V = ds["WS"][j1,j2] - stdV;
            lwpm = allsky_lwp(par);
            par.V = ds["WS"][j1,j2] + stdV;
            lwpp = allsky_lwp(par);
            ddV[i1,i2] = mean(diff([lwpm, lwp0, lwpp]));

            # D
            stdD = ds["D500_std"][j1,j2] / 10;
            par.D = ds["D500"][j1,j2] - stdD;
            lwpm = allsky_lwp(par);
            par.D = ds["D500"][j1,j2] + stdD;
            lwpp = allsky_lwp(par);
            ddD[i1,i2] = mean(diff([lwpm, lwp0, lwpp]));

            # RH
            stdRH = ds["RH500_std"][j1,j2];
            par.RHft = ds["RH500"][j1,j2] - stdRH;
            lwpm = allsky_lwp(par);
            par.RHft = ds["RH500"][j1,j2] + stdRH;
            lwpp = allsky_lwp(par);
            ddRH[i1,i2] = mean(diff([lwpm, lwp0, lwpp]));

            # EIS
            EIS0 = 10.0;
            stdEIS = 1.0;
            par.EIS = EIS0 - stdEIS;
            lwpm = allsky_lwp(par);
            par.EIS = EIS0 + stdEIS;
            lwpp = allsky_lwp(par);
            ddEIS[i1,i2] = mean(diff([lwpm, lwp0, lwpp]));

            # CO2
            stdCO2 = 100.0;
            par.CO2 = par0.CO2 - stdCO2;
            lwpm = allsky_lwp(par);
            par.CO2 = par0.CO2 + stdCO2;
            lwpp = allsky_lwp(par);
            ddCO2[i1,i2] = mean(diff([lwpm, lwp0, lwpp]));
        end
    end
end

println(lon)
println(lat)
println(ddSST)
println(ddV)
println(ddD)
println(ddRH)
println(ddEIS)
println(ddCO2)

cmap = :bluesreds
clims = (-0.2,0.2)

p1 = contourf(lon, lat, ddSST', c=cmap, clims=clims, 
    title="ΔLWP / σSST", xlabel="lon", ylabel="lat",
    left_margin=15Plots.mm, right_margin=10Plots.mm,
    bottom_margin=15Plots.mm, top_margin=10Plots.mm)
p2 = contourf(lon, lat, ddV', c=cmap, clims=clims, 
    title="ΔLWP / σV", xlabel="lon", ylabel="lat",
    left_margin=15Plots.mm, right_margin=10Plots.mm,
    bottom_margin=15Plots.mm, top_margin=10Plots.mm)
p3 = contourf(lon, lat, ddD', c=cmap, clims=clims, 
    title="ΔLWP / (σD500/10)", xlabel="lon", ylabel="lat",
    left_margin=15Plots.mm, right_margin=10Plots.mm,
    bottom_margin=15Plots.mm, top_margin=10Plots.mm)
p4 = contourf(lon, lat, ddRH', c=cmap, clims=clims, 
    title="ΔLWP / σRH500", xlabel="lon", ylabel="lat",
    left_margin=15Plots.mm, right_margin=10Plots.mm,
    bottom_margin=15Plots.mm, top_margin=10Plots.mm)
p5 = contourf(lon, lat, ddEIS', c=cmap, clims=clims, 
    title="ΔLWP / σEIS", xlabel="lon", ylabel="lat",
    left_margin=15Plots.mm, right_margin=10Plots.mm,
    bottom_margin=15Plots.mm, top_margin=10Plots.mm)
p6 = contourf(lon, lat, ddCO2', c=cmap, clims=clims, 
    title="ΔLWP / 100 ppm CO2", xlabel="lon", ylabel="lat",
    left_margin=15Plots.mm, right_margin=10Plots.mm,
    bottom_margin=15Plots.mm, top_margin=10Plots.mm)

plot(p1, p2, p3, p4, p5, p6, layout=(2,3), size=(1600,800), dpi=200)
savefig("experiments/figures/box_linear_perturb_JJA_NEP.png")

close(ds)



# create output netcdf file
mkpath("experiments/tmp/")
ds = Dataset("experiments/tmp/partial_derivatives.nc","c")

# Define the dimension "lon" and "lat" with the size N and M resp.
defDim(ds,"lon",N)
defDim(ds,"lat",M)
defVar(ds,"lon",lon,("lon",))
defVar(ds,"lat",lat,("lat",))

# Define the variables
v = defVar(ds,"dLWPdSST",ddSST,("lon","lat"))
v.attrib["long_name"] = "ΔLWP / σSST"

v = defVar(ds,"dLWPdV",ddV,("lon","lat"))
v.attrib["long_name"] = "ΔLWP / σV"

v = defVar(ds,"dLWPdD",ddD,("lon","lat"))
v.attrib["long_name"] = "ΔLWP / (σD500/10)"

v = defVar(ds,"dLWPdRH",ddRH,("lon","lat"))
v.attrib["long_name"] = "ΔLWP / σRH500"

v = defVar(ds,"dLWPdEIS",ddEIS,("lon","lat"))
v.attrib["long_name"] = "ΔLWP / σEIS"

v = defVar(ds,"dLWPdCO2",ddCO2,("lon","lat"))
v.attrib["long_name"] = "ΔLWP / 100 ppm CO2"

print(ds)
close(ds)