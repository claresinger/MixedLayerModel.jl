using FileIO
using NCDatasets
using Plots
using LaTeXStrings
using MixedLayerModel

# plot LES
ds = Dataset("experiments/LES_steadystate_all_upsteps.nc");
max = 4
co2 = Float64.(ds["CO2"][1:max]);
zi = Float64.(ds["zi"][1:max]);
zb = Float64.(ds["zb"][1:max]);
we = zi.*6e-6*1e3;
lwp = Float64.(ds["lwp"][1:max]*1e3);
sst = Float64.(ds["sst"][1:max]);
lhf = Float64.(ds["lhf"][1:max]);

p1 = scatter(co2, zi, marker=:x, markersize=5, label="LES", ylabel="Inversion height, zi [m]")
p2 = scatter(co2, zb, marker=:x, markersize=5, label="", ylabel="Cloud base, zb [m]")
p3 = scatter(co2, we, marker=:x, markersize=5, label="", ylabel="Entrainment velocity, we [mm/s]")
p4 = scatter(co2, lwp, marker=:x, markersize=5, label="", xlabel="CO2 [ppmv]", ylabel="LWP [g/m2]")
p5 = scatter(co2, sst, marker=:x, markersize=5, label="", xlabel="CO2 [ppmv]", ylabel="SST [K]")
p6 = scatter(co2, lhf, marker=:x, markersize=5, label="", xlabel="CO2 [ppmv]", ylabel="LHF [W/m2]")

# plot MLM
# exp_path = "new_alpha/"
# co2 = [400, 500, 600, 700, 800];
exp_path = "new_alpha_enBal/"
co2 = [400, 500, 600, 700, 800, 900, 1000];

zi, zb, we = zeros(length(co2)), zeros(length(co2)), zeros(length(co2));
lwp, sst, lhf = zeros(length(co2)), zeros(length(co2)), zeros(length(co2));
for (i, co2i) in enumerate(co2)
    if co2i == 400.0
        file = "experiments/output/"*exp_path*"co2_400.jld2"
    else
        file = "experiments/output/"*exp_path*"co2_upstep_"*string(co2i)*".jld2"
    end
    dat = load(file);
    zii, hM, qM, ssti = dat["uf"];
    zbi = dat["zb"];
    lwpi = calc_LWP(zii, hM, qM) * 1e3;
    zi[i], zb[i], lwp[i], sst[i] = zii, zbi, lwpi, ssti;
    lhf[i], we[i] = dat["LHF"], dat["we"]*1e3;
end
scatter!(p1, co2, zi, marker=:o, markersize=5, label="MLM with slab ocean")
scatter!(p2, co2, zb, marker=:o, markersize=5, label="")
scatter!(p3, co2, we, marker=:o, markersize=5, label="")
scatter!(p4, co2, lwp, marker=:o, markersize=5, label="")
scatter!(p5, co2, sst, marker=:o, markersize=5, label="")
scatter!(p6, co2, lhf, marker=:o, markersize=5, label="")

# save plot
p = plot(p1,p2,p3,p4,p5,p6, layout=(2,3), link=:x, size=(1000,500), dpi=300);
mkpath("experiments/figures/"*exp_path)
savefig(p, "experiments/figures/"*exp_path*"steady-state.png")