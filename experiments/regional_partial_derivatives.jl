push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using MixedLayerModel
using MixedLayerModel: Cp, g
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

function slope(xi,yi)
    return sum((xi.-mean(xi)).*(yi.-mean(yi))) / sum((xi.-mean(xi)).^2)
end

# load boundary conditions from file
file = "experiments/data/box_BCs_JJA_NEP.nc";
ds = Dataset(file, "r");
# println(ds)

par0 = climatology();

par = climatology();
par.rtype = varRad();
par.stype = fixSST();
par.ftype = varFlux();
par.fttype = fixedFT();
par.etype = enBal();

skip1 = 4
skip2 = 2
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
    
            # par.OHU = -5;
            # par.LHF = ds["LHF"][j1,j2];
            # par.SHF = 1.0;
            # par.SHF = ds["SHF"][j1,j2];
            
            par.D = ds["D500"][j1,j2];
            par.RHft = ds["RH500"][j1,j2];

            EIS = 10.0; # (K)
            par.sft0 = Cp*(par.SST0 + EIS); # 10 (K) jump
            par.Gamma_s = (Cp*-5e-3) + g; # (K/m)

            # base state
            lwp0 = allsky_lwp(par);

            # sst
            stdSST = 1.0;
            par.SST0 = ds["sst"][i1,i2] - stdSST;
            lwpm = allsky_lwp(par);
            par.SST0 = ds["sst"][i1,i2] + stdSST;
            lwpp = allsky_lwp(par);
            ddSST[i1,i2] = mean(diff([lwpm, lwp0, lwpp]));
            # ddSST[i1,i2] = slope(ds["sst"][i1,i2] .* [0.99, 1, 1.01], [lwpm, lwp0, lwpp])

            # V
            stdV = 0.1;
            par.V = par0.V - stdV;
            lwpm = allsky_lwp(par);
            par.V = par0.V + stdV;
            lwpp = allsky_lwp(par);
            ddV[i1,i2] = mean(diff([lwpm, lwp0, lwpp]));
            # ddV[i1,i2] = slope(par0.V .* [0.99, 1, 1.01], [lwpm, lwp0, lwpp])

            # D
            stdD = 0.1e-6;
            par.D = ds["D500"][i1,i2] - stdD;
            lwpm = allsky_lwp(par);
            par.D = ds["D500"][i1,i2] + stdD;
            lwpp = allsky_lwp(par);
            ddD[i1,i2] = mean(diff([lwpm, lwp0, lwpp]));
            # ddD[i1,i2] = slope(ds["D500"][i1,i2] .* [0.99, 1, 1.01], [lwpm, lwp0, lwpp])

            # RH
            stdRH = 0.01
            par.RHft = ds["RH500"][i1,i2] - stdRH;
            lwpm = allsky_lwp(par);
            par.RHft = ds["RH500"][i1,i2] + stdRH;
            lwpp = allsky_lwp(par);
            ddRH[i1,i2] = mean(diff([lwpm, lwp0, lwpp]));
            # ddRH[i1,i2] = slope(ds["RH500"][i1,i2] .* [0.99, 1, 1.01], [lwpm, lwp0, lwpp])

            # EIS
            EIS0 = 10.0;
            stdEIS = 1.0;
            par.sft0 = Cp*(par.SST0 + EIS0 - stdEIS);
            lwpm = allsky_lwp(par);
            par.sft0 = Cp*(par.SST0 + EIS0 + stdEIS);
            lwpp = allsky_lwp(par);
            ddEIS[i1,i2] = mean(diff([lwpm, lwp0, lwpp]));
            # ddEIS[i1,i2] = slope(EIS0 .* [0.99, 1, 1.01], [lwpm, lwp0, lwpp])

            # CO2
            stdCO2 = 5.0;
            par.CO2 = par0.CO2 - stdCO2;
            lwpm = allsky_lwp(par);
            par.CO2 = par0.CO2 + stdCO2;
            lwpp = allsky_lwp(par);
            ddCO2[i1,i2] = mean(diff([lwpm, lwp0, lwpp]));
            # ddCO2[i1,i2] = slope(par0.CO2 .* [0.99, 1, 1.01], [lwpm, lwp0, lwpp])
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

x = ddSST
clims = (-maximum(abs.(x)), maximum(abs.(x)))
println(clims)
p1 = contourf(lon, lat, ddSST, c=cmap, clims=clims, 
    title="ΔLWP / σSST", xlabel="lon", ylabel="lat",
    left_margin=15Plots.mm, right_margin=10Plots.mm,
    bottom_margin=15Plots.mm, top_margin=10Plots.mm)
x = ddV
clims = (-maximum(abs.(x)), maximum(abs.(x)))
println(clims)
p2 = contourf(lon, lat, ddV, c=cmap, clims=clims, 
    title="ΔLWP / σV", xlabel="lon", ylabel="lat",
    left_margin=15Plots.mm, right_margin=10Plots.mm,
    bottom_margin=15Plots.mm, top_margin=10Plots.mm)
x = ddD
clims = (-maximum(abs.(x)), maximum(abs.(x)))
println(clims)
p3 = contourf(lon, lat, ddD, c=cmap, clims=clims, 
    title="ΔLWP / σD", xlabel="lon", ylabel="lat",
    left_margin=15Plots.mm, right_margin=10Plots.mm,
    bottom_margin=15Plots.mm, top_margin=10Plots.mm)
x = ddRH
clims = (-maximum(abs.(x)), maximum(abs.(x)))
println(clims)
p4 = contourf(lon, lat, ddRH, c=cmap, clims=clims, 
    title="ΔLWP / σRH", xlabel="lon", ylabel="lat",
    left_margin=15Plots.mm, right_margin=10Plots.mm,
    bottom_margin=15Plots.mm, top_margin=10Plots.mm)
x = ddEIS
clims = (-maximum(abs.(x)), maximum(abs.(x)))
println(clims)
p5 = contourf(lon, lat, ddEIS, c=cmap, clims=clims, 
    title="ΔLWP / σEIS", xlabel="lon", ylabel="lat",
    left_margin=15Plots.mm, right_margin=10Plots.mm,
    bottom_margin=15Plots.mm, top_margin=10Plots.mm)
x = ddCO2
clims = (-maximum(abs.(x)), maximum(abs.(x)))
println(clims)
p6 = contourf(lon, lat, ddCO2, c=cmap, clims=clims, 
    title="ΔLWP / σCO2", xlabel="lon", ylabel="lat",
    left_margin=15Plots.mm, right_margin=10Plots.mm,
    bottom_margin=15Plots.mm, top_margin=10Plots.mm)

plot(p1, p2, p3, p4, p5, p6, layout=(2,3), size=(1600,800), dpi=200)
savefig("experiments/figures/box_linear_perturb_JJA_NEP.png")

close(ds)
