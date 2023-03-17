using FileIO
using NCDatasets
using MixedLayerModel
using Plots
using Statistics

Cd = 4

path = "experiments/figures/20221115_dailytransect_subonly/"
file = path*"transect_output_Cd"*string(Cd)*".nc";
ds = Dataset(file, "r");
uf_save = cat(ds["zi"], ds["sM"], ds["qtM"], ds["sst"], ds["cf"], dims=3)
println(size(uf_save))
Ndays, Nlons, _ = size(uf_save)

lon = ds["lon"]

# set up MLM params
par = climatology();
par.rtype = varRad();
par.stype = fixSST();
par.ftype = varFlux();
par.fttype = fixEIS();
par.etype = enBal();
par.decoup_slope = 8;
par.α_vent = 1.22e-3;
par.Cd = Cd * 1e-4;

LHF = zeros((Ndays,Nlons))
for i in 1:Ndays
    for j in 1:Nlons
        uf = uf_save[i,j,:]
        LHF[i,j] = calc_LHF(uf, par)
    end
end
plot(dpi=500)
println(mean(LHF, dims=1))
plot!(lon, mean(LHF, dims=1)', lw=3, color=:black, label=false)
for i in 1:Ndays
    plot!(lon, LHF[i,:], label=false)
end
xlabel!("Longitude")
xlims!(-150,-120)
ylabel!("LHF [W/m2]")
savefig(path*"LHF_Cd"*string(Cd)*".png")