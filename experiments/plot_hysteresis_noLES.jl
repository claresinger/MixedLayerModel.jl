using FileIO
using Plots
using MixedLayerModel
using MixedLayerModel: Rd, Rv, L0, T0, Cp, δ, ϵ, μ

ENV["GKSwstype"]="nul"
Plots.scalefontsizes(1.8)

##############
#plot modified results

exp_path = ARGS[1];
co2u = ARGS[2];
co2u = parse.(Int64, split(replace(replace(co2u, "]"=>""), "["=>""),","))
co2d = ARGS[3];
co2d = parse.(Int64, split(replace(replace(co2d, "]"=>""), "["=>""),","))
print(co2u)
println(typeof(co2u))

N = length(co2u);
zi, zb, ent = zeros(N), zeros(N), zeros(N);
cf, lwp, sst, lhf = zeros(N), zeros(N), zeros(N), zeros(N);
dR = zeros(N);

for (i, co2i) in enumerate(co2u)
    file = "experiments/output/"*exp_path*"co2_upstep_"*string(co2i)*".jld2"
    dat = load(file);
    uf = dat["uf"];
    par = dat["p"];
    zii, sM, qM, ssti, cfi = uf;
    zbi = dat["zb"];
    zi[i], zb[i], sst[i], cf[i] = zii, zbi, ssti, cfi;
    lhf[i], ent[i], dR[i] = dat["LHF"], dat["we"]*1e3, dat["ΔR"];
    lwp[i] = incloud_LWP(uf, zb[i]) * 1e3;
end
S = (lhf./dR).*((zi.-zb)./zi);

ms = 8
c = "crimson"
p_dR = scatter(co2u, dR, color=c, marker=:circle, markersize=ms, markerstrokewidth=0,
                label="", ylabel="ΔR [W/m²]",
                xticks=xtks, xlim=xrange, ylim=[0,100])
p_decoup = scatter(co2u, S, color=c, marker=:circle, markersize=ms, markerstrokewidth=0,
                label="", xlabel="CO₂ [ppmv]", ylabel="Decoupling, \$\\mathcal{D}\$",
                yscale=:log10, yticks=([0.1,1,10], ["0.1","1","10"]), ylim=[0.1,10],
                xticks=xtks, xlim=xrange)
p_sst = scatter(co2u, sst, color=c, marker=:circle, markersize=ms, markerstrokewidth=0,
                label="", ylabel="SST [K]",
                xticks=xtks, xlim=xrange, ylim=SSTrange)
p_lhf = scatter(co2u, lhf, color=c, marker=:circle, markersize=ms, markerstrokewidth=0,
                label="", xlabel="CO₂ [ppmv]", ylabel="LHF [W/m²]",
                xticks=xtks, xlim=xrange, ylim=LHFrange)
p_zi = scatter(co2u, zi, color=c, marker=:circle, markersize=ms, markerstrokewidth=0,
                label="", ylabel="zᵢ [m]",
                xticks=xtks, xlim=xrange, ylim=[0, 1000])
p_cf = scatter(co2u, cf*100, color=c, marker=:circle, markersize=ms, markerstrokewidth=0,
                label="", xlabel="CO₂ [ppmv]", ylabel="CF [%]",
                xticks=xtks, xlim=xrange, ylim=[0,100])

plot!(p_dR, co2u, dR, color=c, linewidth=2, label="")
plot!(p_decoup, co2u, S, color=c, linewidth=2, label="")
plot!(p_sst, co2u, sst, color=c, linewidth=2, label="")
plot!(p_lhf, co2u, lhf, color=c, linewidth=2, label="")
plot!(p_zi, co2u, zi, color=c, linewidth=2, label="")
plot!(p_cf, co2u, cf*100, color=c, linewidth=2, label="")

N = length(co2d);
zi, zb, ent = zeros(N), zeros(N), zeros(N);
cf, lwp, sst, lhf = zeros(N), zeros(N), zeros(N), zeros(N);
dR = zeros(N);

for (i, co2i) in enumerate(co2d)
    file = "experiments/output/"*exp_path*"co2_downstep_"*string(co2i)*".jld2"
    dat = load(file);
    uf = dat["uf"];
    par = dat["p"];
    zii, sM, qM, ssti, cfi = uf;
    zbi = dat["zb"];
    zi[i], zb[i], sst[i], cf[i] = zii, zbi, ssti, cfi;
    lhf[i], ent[i], dR[i] = dat["LHF"], dat["we"]*1e3, dat["ΔR"];
    lwp[i] = incloud_LWP(uf, zb[i]) * 1e3;
end
S = (lhf./dR).*((zi.-zb)./zi);

c = "royalblue"
scatter!(p_dR, co2d, dR, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_decoup, co2d, S, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_sst, co2d, sst, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_lhf, co2d, lhf, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_zi, co2d, zi, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_cf, co2d, cf*100, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
plot!(p_dR, co2d, dR, color=c, linewidth=2, label="")
plot!(p_decoup, co2d, S, color=c, linewidth=2, label="")
plot!(p_sst, co2d, sst, color=c, linewidth=2, label="")
plot!(p_lhf, co2d, lhf, color=c, linewidth=2, label="")
plot!(p_zi, co2d, zi, color=c, linewidth=2, label="")
plot!(p_cf, co2d, cf*100, color=c, linewidth=2, label="")


##############
# plot orig results
orig_co2u = parse.(Int64, split(replace(replace(orig_co2u, "]"=>""), "["=>""),","))
orig_co2d = parse.(Int64, split(replace(replace(orig_co2d, "]"=>""), "["=>""),","))

N = length(orig_co2u);
zi, zb, ent = zeros(N), zeros(N), zeros(N);
cf, lwp, sst, lhf = zeros(N), zeros(N), zeros(N), zeros(N);
dR = zeros(N);

for (i, co2i) in enumerate(orig_co2u)
    file = "experiments/output/"*orig_path*"co2_upstep_"*string(co2i)*".jld2"
    dat = load(file);
    uf = dat["uf"];
    par = dat["p"];
    zii, sM, qM, ssti, cfi = uf;
    zbi = dat["zb"];
    zi[i], zb[i], sst[i], cf[i] = zii, zbi, ssti, cfi;
    lhf[i], ent[i], dR[i] = dat["LHF"], dat["we"]*1e3, dat["ΔR"];
    lwp[i] = incloud_LWP(uf, zb[i]) * 1e3;
end
S = (lhf./dR).*((zi.-zb)./zi);

c = "crimson"
ms = 3
scatter!(p_dR, orig_co2u, dR, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_decoup, orig_co2u, S, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_sst, orig_co2u, sst, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_lhf, orig_co2u, lhf, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_zi, orig_co2u, zi, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_cf, orig_co2u, cf*100, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
plot!(p_dR, orig_co2u, dR, color=c, linewidth=1, alpha=0.5, label="")
plot!(p_decoup, orig_co2u, S, color=c, linewidth=1, alpha=0.5, label="")
plot!(p_sst, orig_co2u, sst, color=c, linewidth=1, alpha=0.5, label="")
plot!(p_lhf, orig_co2u, lhf, color=c, linewidth=1, alpha=0.5, label="")
plot!(p_zi, orig_co2u, zi, color=c, linewidth=1, alpha=0.5, label="")
plot!(p_cf, orig_co2u, cf*100, color=c, linewidth=1, alpha=0.5, label="")


N = length(orig_co2d);
zi, zb, ent = zeros(N), zeros(N), zeros(N);
cf, lwp, sst, lhf = zeros(N), zeros(N), zeros(N), zeros(N);
dR = zeros(N);

for (i, co2i) in enumerate(orig_co2d)
    file = "experiments/output/"*orig_path*"co2_downstep_"*string(co2i)*".jld2"
    dat = load(file);
    uf = dat["uf"];
    par = dat["p"];
    zii, sM, qM, ssti, cfi = uf;
    zbi = dat["zb"];
    zi[i], zb[i], sst[i], cf[i] = zii, zbi, ssti, cfi;
    lhf[i], ent[i], dR[i] = dat["LHF"], dat["we"]*1e3, dat["ΔR"];
    lwp[i] = incloud_LWP(uf, zb[i]) * 1e3;
end
S = (lhf./dR).*((zi.-zb)./zi);

c = "royalblue"
ms = 3
scatter!(p_dR, orig_co2d, dR, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_decoup, orig_co2d, S, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_sst, orig_co2d, sst, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_lhf, orig_co2d, lhf, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_zi, orig_co2d, zi, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_cf, orig_co2d, cf*100, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
plot!(p_dR, orig_co2d, dR, color=c, linewidth=1, alpha=0.5, label="")
plot!(p_decoup, orig_co2d, S, color=c, linewidth=1, alpha=0.5, label="")
plot!(p_sst, orig_co2d, sst, color=c, linewidth=1, alpha=0.5, label="")
plot!(p_lhf, orig_co2d, lhf, color=c, linewidth=1, alpha=0.5, label="")
plot!(p_zi, orig_co2d, zi, color=c, linewidth=1, alpha=0.5, label="")
plot!(p_cf, orig_co2d, cf*100, color=c, linewidth=1, alpha=0.5, label="")

##############

# # legend
# plot!(p_dR, [-1, -1], [-1, -1], color="black", alpha=0.5, label="Orig")
# scatter!(p_dR, [-1], [-1], color="black", marker=:circle, markersize=1, markerstrokewidth=0, label="WV off")

# save plot
p = plot(p_dR,p_sst,p_zi,p_decoup,p_lhf,p_cf, layout=(2,3), 
    link=:x, size=(1300,650), dpi=300,
    legend=:topright, legendfontsize=12, legendfont=font(12),
    left_margin=10Plots.mm, bottom_margin=7Plots.mm, top_margin=5Plots.mm);
mkpath("experiments/figures/"*exp_path)
savefig(p, "experiments/figures/"*exp_path*"hystersis_plot.png")

# save another plot
p = plot(p_dR,p_sst,p_lhf,p_cf, layout=(2,2), 
    link=:x, size=(900,500), dpi=300,
    legend=:topright, legendfontsize=10, legendfont=font(10),
    left_margin=10Plots.mm, bottom_margin=7Plots.mm, top_margin=5Plots.mm);
mkpath("experiments/figures/"*exp_path)
savefig(p, "experiments/figures/"*exp_path*"hystersis_plot_"*title*".png")

# reset fontsizes
Plots.scalefontsizes(1/1.8)