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
S = (lhf./dR).*((zi.-zb)./zi);

ms = 8
c = "crimson"
p1 = scatter(site, dR, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, 
                label="", ylabel="ΔR [W/m\$^2\$]")
p2 = scatter(site, S, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, 
                label="", ylabel="Stability, \$S\$")
p3 = scatter(site, cf*100, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, 
                label="", ylabel="CF [%]")
p4 = scatter(site, sst, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, 
                label="", ylabel="SST [K]")
plot!(p1, site, dR, color=c, linewidth=2, label="")
plot!(p2, site, S, color=c, linewidth=2, label="")
plot!(p3, site, cf*100, color=c, linewidth=2, label="")
plot!(p4, site, sst, color=c, linewidth=2, label="")

p = plot(p1,p2,p3,p4, layout=(2,2), 
    link=:x, size=(1000,650), dpi=300,
    legend=:topright, legendfontsize=12, legendfont=font(12),
    left_margin=10Plots.mm, bottom_margin=5Plots.mm, top_margin=5Plots.mm);
mkpath("experiments/figures/"*exp_path)
savefig(p, "experiments/figures/"*exp_path*"NEP_transect.png")