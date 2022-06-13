push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using FileIO
using NCDatasets
using Plots
using LaTeXStrings
using MixedLayerModel
using MixedLayerModel: Rd, Rv, L0, T0, Cp, δ, ϵ, μ

ENV["GKSwstype"]="nul"

# plot LES
ds = Dataset("experiments/LES_steadystate_all_upsteps.nc");
max = 6;
co2 = Float64.(ds["CO2"][1:max]);
zi = Float64.(ds["zi"][1:max]);
zb = Float64.(ds["zb"][1:max]);
ent = zi.*6e-6*1e3;
lwp = Float64.(ds["lwp"][1:max]*1e3);
sst = Float64.(ds["sst"][1:max]);
lhf = Float64.(ds["lhf"][1:max]);
dR = Float64.(ds["deltaR"][1:max]);
dR = [74.1,69.6,66.5,61.9,4.7,2.5];
S = (lhf./dR).*((zi.-zb)./zi);
hplus_les = Float64.(ds["h_plus"][1:max]);
hM_les = Float64.(ds["hM"][1:max]);
qplus_les = Float64.(ds["qt_plus"][1:max]);
qM_les = Float64.(ds["qtM"][1:max]);
Δs_vli = ((hplus_les-hM_les) - μ*L0*(qplus_les-qM_les))*1e-3;

p1 = scatter(co2, zi, marker=:x, markersize=5, label="LES", ylabel="")
scatter!(co2, zb, marker=:x, color="green", markersize=5, label="", ylabel="Cloud top/base [m]")
p2 = scatter(co2, zi-zb, marker=:x, markersize=5, label="", ylabel="Cloud depth, zc [m]")
p3 = scatter(co2, lwp, marker=:x, markersize=5, label="", ylabel="LWP [g/m2]")

p4 = scatter(co2, ent, marker=:x, markersize=5, label="", ylabel="Entrainment, we [mm/s]")
p5 = scatter(co2, Δs_vli, marker=:x, markersize=5, label="", ylabel="Inv. strength [kJ/kg]")
p6 = scatter(co2, sst, marker=:x, markersize=5, label="", ylabel="SST [K]")

p7 = scatter(co2, lhf, marker=:x, markersize=5, label="", xlabel="CO2 [ppmv]", ylabel="LHF [W/m2]")
p8 = scatter(co2, dR, marker=:x, markersize=5, label="", xlabel="CO2 [ppmv]", ylabel="dR [W/m2]")
p9 = scatter(co2, S, marker=:x, markersize=5, label="", 
                yscale=:log10, yticks=[0.1,0.5,1,5,10], yformatter=:plain,
                xlabel="CO2 [ppmv]", ylabel="Stability param., S")

# p1 = scatter(co2, zi, marker=:x, color="white", markersize=5, label="", ylabel="")
# scatter!(co2, zb, marker=:x, color="white", markersize=5, label="", ylabel="Cloud top/base [m]")
# p2 = scatter(co2, zi-zb, marker=:x, color="white", markersize=5, label="", ylabel="Cloud depth, zc [m]")
# p3 = scatter(co2, lwp, marker=:x, color="white", markersize=5, label="", ylabel="LWP [g/m2]")

# p4 = scatter(co2, we, marker=:x, color="white", markersize=5, label="", ylabel="Entrainment, we [mm/s]")
# p5 = scatter(co2, Δs_vli, marker=:x, color="white", markersize=5, label="", ylabel="Inv. strength [kJ/kg]")
# p6 = scatter(co2, sst, marker=:x, color="white", markersize=5, label="", ylabel="SST [K]")

# p7 = scatter(co2, lhf, marker=:x, color="white", markersize=5, label="", xlabel="CO2 [ppmv]", ylabel="LHF [W/m2]")
# p8 = scatter(co2, dR, marker=:x, color="white", markersize=5, label="", xlabel="CO2 [ppmv]", ylabel="dR [W/m2]")
# p9 = scatter(co2, S, marker=:x, color="white", markersize=5, label="", xlabel="CO2 [ppmv]", ylabel="Stability param., S")

# plot MLM
# exp_path = "new_alpha/"
# co2 = [400, 500, 600, 700, 800];
# exp_path = "new_alpha_enBal_invco2/"
# co2 = [400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400];
# exp_path = "enBal_restart/"
# co2 = [400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000];

# exp_path = "dampSST/"
# co2 = [400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600];
# exp_path = "twocol/"
# co2 = [400, 600, 800, 1000, 1200, 1600, 1800, 2000, 2200, 2400, 2500, 2600];
exp_path = "modCF/"
co2 = [400, 600, 800, 1000, 1200, 1300, 1400, 1600];

zi, zb, ent = zeros(length(co2)), zeros(length(co2)), zeros(length(co2));
cf, lwp, sst, lhf = zeros(length(co2)), zeros(length(co2)), zeros(length(co2)), zeros(length(co2));
dR, Δs_vli = zeros(length(co2)), zeros(length(co2));
for (i, co2i) in enumerate(co2)
    if co2i == 400.0
        file = "experiments/output/"*exp_path*"co2_400.jld2"
    else
        file = "experiments/output/"*exp_path*"co2_upstep_"*string(co2i)*".jld2"
    end
    dat = load(file);
    uf = dat["uf"];
    p = dat["p"];
    zii, sM, qM, ssti, cfi = uf;
    zbi = dat["zb"];
    zi[i], zb[i], sst[i] = zii, zbi, ssti;
    lhf[i], ent[i], dR[i] = dat["LHF"], dat["we"]*1e3, dat["ΔR"];
    LWPi = incloud_LWP(uf, zbi);
    Δs_vli[i] = sv_jump(uf, p, LWPi)*1e-3;

    cf[i] = cfi;
    lwp[i] = incloud_LWP(uf, zb[i]) * 1e3;
end
S = (lhf./dR).*((zi.-zb)./zi);

ms = 4
scatter!(p1, co2, zi, marker=:circle, color=2, markersize=ms, markerstrokewidth=0, label="MLM with slab ocean")
scatter!(p1, co2, zb, marker=:circle, color="red", markersize=ms, markerstrokewidth=0, label="")
scatter!(p2, co2, zi-zb, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p3, co2, lwp .* cf, marker=:circle, markersize=ms, markerstrokewidth=0, label="")

scatter!(p4, co2, ent, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p5, co2, Δs_vli, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p6, co2, sst, marker=:circle, markersize=ms, markerstrokewidth=0, label="")

scatter!(p7, co2, lhf, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p8, co2, dR, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p9, co2, S, marker=:circle, markersize=ms, markerstrokewidth=0, label="")

# save plot
p = plot(p1,p2,p3,p4,p5,p6,p7,p8,p9, layout=(3,3), 
    link=:x, size=(1000,700), dpi=300,
    legend=:bottomleft,
    left_margin=10Plots.mm, bottom_margin=5Plots.mm, top_margin=5Plots.mm);
mkpath("experiments/figures/"*exp_path)
savefig(p, "experiments/figures/"*exp_path*"steady-state.png")