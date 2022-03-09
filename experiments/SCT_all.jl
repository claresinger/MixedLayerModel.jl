push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using FileIO
using Plots
using MixedLayerModel

# # run transition from coast to West
# ARGS = ["1"]; include("SCT_NEP.jl");
# ARGS = ["2"]; include("SCT_NEP.jl");
# ARGS = ["3"]; include("SCT_NEP.jl");
# ARGS = ["4"]; include("SCT_NEP.jl");

# plot values along transect
exp_path = "SCT_NEP/"
site = [1, 2, 3, 4];
N = length(site);

zi, zb, ent = zeros(N), zeros(N), zeros(N);
cf, lwp, sst, lhf = zeros(N), zeros(N), zeros(N), zeros(N);
dR, Δs_vli = zeros(N), zeros(N);

for (i, sitei) in enumerate(site)
    file = "experiments/output/"*exp_path*"sol_"*string(sitei)*".jld2"
    dat = load(file);
    uf = dat["uf"];
    par = dat["p"];
    zii, hM, qM, ssti, cfi = uf;
    zbi = dat["zb"];
    zi[i], zb[i], sst[i], cf[i] = zii, zbi, ssti, cfi;
    lhf[i], ent[i], dR[i] = dat["LHF"], dat["we"]*1e3, dat["ΔR"];
    lwp[i] = incloud_LWP(uf, zb[i]) * 1e3;
    Δs_vli[i] = Δs(uf, par, lwp[i])*1e-3;
end
I = (lhf./dR).*((zi.-zb)./zi);

ms = 8
c = "crimson"
p1 = scatter(site, sst, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, 
                label="", ylabel="SST [K]")
p2 = scatter(site, lhf, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, 
                label="", ylabel="LHF [W/m\$^2\$]")
p3 = scatter(site, I, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, 
                label="", ylabel="I [-]")
p4 = scatter(site, cf*100, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, 
                label="", ylabel="CF [%]")
p5 = scatter(site, zi, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, 
                label="", ylabel="Inversion height [m]")
p6 = scatter(site, zi-zb, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, 
                label="", ylabel="Cloud depth [m]")


### observations
sst_obs = [290.2, 293.5, 297.2, 299.5];
cf_obs = [21.60, 62.26, 31.77, 13.18];

c = "dodgerblue"
scatter!(p1, site, sst_obs, color=c, marker=:x, markersize=ms, markerstrokewidth=0, label="")
#scatter!(p2, site, lhf_obs, color=c, marker=:x, markersize=ms, markerstrokewidth=0, label="")
#scatter!(p3, site, I_obs, color=c, marker=:x, markersize=ms, markerstrokewidth=0, label="")
scatter!(p4, site, cf_obs, color=c, marker=:x, markersize=ms, markerstrokewidth=0, label="")
#scatter!(p5, site, zi_obs, color=c, marker=:x, markersize=ms, markerstrokewidth=0, label="")
#scatter!(p6, site, zi_obs - zb_obs, color=c, marker=:x, markersize=ms, markerstrokewidth=0, label="")

p = plot(p1,p2,p3,p4,p5,p6, layout=(2,3), 
    link=:x, size=(1000,650), dpi=300,
    legend=:topright, legendfontsize=12, legendfont=font(12),
    left_margin=10Plots.mm, bottom_margin=5Plots.mm, top_margin=5Plots.mm);
mkpath("experiments/figures/"*exp_path)
savefig(p, "experiments/figures/"*exp_path*"NEP_transect.png")