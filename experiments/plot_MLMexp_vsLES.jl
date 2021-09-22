using FileIO
using NCDatasets
using Plots
using LaTeXStrings
using MixedLayerModel

ENV["GKSwstype"]="nul"

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
dR = Float64.(ds["deltaR"][1:max]);
S = (lhf./dR).*((zi.-zb)./zi);

p1 = scatter(co2, zi, marker=:x, markersize=5, label="LES", ylabel="Inversion height, zi [m]")
p2 = scatter(co2, zb, marker=:x, markersize=5, label="", ylabel="Cloud base, zb [m]")
p3 = scatter(co2, zi-zb, marker=:x, markersize=5, label="", ylabel="Cloud depth, zc [m]")

p4 = scatter(co2, we, marker=:x, markersize=5, label="", ylabel="Entrainment, we [mm/s]")
p5 = scatter(co2, lwp, marker=:x, markersize=5, label="", xlabel="CO2 [ppmv]", ylabel="LWP [g/m2]")
p6 = scatter(co2, sst, marker=:x, markersize=5, label="", xlabel="CO2 [ppmv]", ylabel="SST [K]")

p7 = scatter(co2, lhf, marker=:x, markersize=5, label="", xlabel="CO2 [ppmv]", ylabel="LHF [W/m2]")
p8 = scatter(co2, dR, marker=:x, markersize=5, label="", xlabel="CO2 [ppmv]", ylabel="dR [W/m2]")
p9 = scatter(co2, S, marker=:x, markersize=5, label="", xlabel="CO2 [ppmv]", ylabel="Stability param., S")

# plot MLM
# exp_path = "new_alpha/"
# co2 = [400, 500, 600, 700, 800];
exp_path = "new_alpha_enBal_invco2/"
co2 = [400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400];

zi, zb, we = zeros(length(co2)), zeros(length(co2)), zeros(length(co2));
lwp, sst, lhf = zeros(length(co2)), zeros(length(co2)), zeros(length(co2));
dR = zeros(length(co2));
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
    lhf[i], we[i], dR[i] = dat["LHF"], dat["we"]*1e3, dat["ΔR"];
end
S = (lhf./dR).*((zi.-zb)./zi);

scatter!(p1, co2, zi, marker=:circle, markersize=5, markerstrokewidth=0, label="MLM with slab ocean")
scatter!(p2, co2, zb, marker=:circle, markersize=5, markerstrokewidth=0, label="")
scatter!(p3, co2, zi-zb, marker=:circle, markersize=5, markerstrokewidth=0, label="")

scatter!(p4, co2, we, marker=:circle, markersize=5, markerstrokewidth=0, label="")
scatter!(p5, co2, lwp, marker=:circle, markersize=5, markerstrokewidth=0, label="")
scatter!(p6, co2, sst, marker=:circle, markersize=5, markerstrokewidth=0, label="")

scatter!(p7, co2, lhf, marker=:circle, markersize=5, markerstrokewidth=0, label="")
scatter!(p8, co2, dR, marker=:circle, markersize=5, markerstrokewidth=0, label="")
scatter!(p9, co2, S, marker=:circle, markersize=5, markerstrokewidth=0, label="")

# save plot
p = plot(p1,p2,p3,p4,p5,p6,p7,p8,p9, layout=(3,3), 
    link=:x, size=(1000,700), dpi=300,
    left_margin=10Plots.mm, bottom_margin=5Plots.mm, top_margin=5Plots.mm);
mkpath("experiments/figures/"*exp_path)
savefig(p, "experiments/figures/"*exp_path*"steady-state.png")
