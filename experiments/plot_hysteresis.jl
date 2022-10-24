push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using FileIO
using Plots
using MixedLayerModel
using MixedLayerModel: Rd, Rv, L0, T0, Cp, δ, ϵ, μ

ENV["GKSwstype"]="nul"
Plots.scalefontsizes(1.8)

# plot LES
co2 = [200,300,400,800,1000,1200,1300,1400,1600];
sst = [287.7,289.1,290.0,292.2,293.2,294.3,304.5,305.8,308.0];
cf = [1.00,1.00,1.00,1.00,1.00,0.99,0.25,0.22,0.19];
zi = [1442,1349,1266,1078,1011,972,834,781,703];
zb = [996,939,895,717,650,621,429,413,392];
lhf = [97.8,103.9,107.1,112.8,115.3,120.5,208.7,213.5,220.9];
dR = [80.2,76.0,74.1,69.6,66.5,61.9,4.7,3.9,2.5];
S = (lhf./dR).*((zi.-zb)./zi);

ms = 10
c = "crimson"
p_dR = scatter(co2, dR, color=c, marker=:x, markersize=ms, label="", 
                ylabel="ΔR [W/m²]",
                xticks=([400, 800, 1200, 1600]), xlim=[100,2000], ylim=[0,100])
p_decoup = scatter(co2, S, color=c, marker=:x, markersize=ms, label="", 
                xlabel="CO₂ [ppmv]", ylabel="Decoupling, \$\\mathcal{D}\$",
                yscale=:log10, yticks=([0.1,1,10], ["0.1","1","10"]), ylim=[0.1,40],
                xticks=([400, 800, 1200, 1600]), xlim=[100,2000])
p_sst = scatter(co2, sst, color=c, marker=:x, markersize=ms, label="", 
                ylabel="SST [K]",
                xticks=([400, 800, 1200, 1600]), xlim=[100,2000], ylim=[286, 312])
p_lhf = scatter(co2, lhf, color=c, marker=:x, markersize=ms, label="",
                xlabel="CO₂ [ppmv]", ylabel="LHF [W/m²]",
                xticks=([400, 800, 1200, 1600]), xlim=[100,2000], ylim=[50, 250])
p_zi = scatter(co2, zi, color=c, marker=:x, markersize=ms, label="",
                ylabel="zᵢ [m]",
                xticks=([400, 800, 1200, 1600]), xlim=[100,2000], ylim=[500, 1500])
p_cf = scatter(co2, cf*100, color=c, marker=:x, markersize=ms, label="", 
                xlabel="CO₂ [ppmv]", ylabel="CF [%]",
                xticks=([400, 800, 1200, 1600]), xlim=[100,2000], ylim=[0,110])

plot!(p_dR, co2, dR, linewidth=2, linestyle=:dot, color=c, label="")
plot!(p_decoup, co2, S, linewidth=2, linestyle=:dot, color=c, label="")
plot!(p_sst, co2, sst, linewidth=2, linestyle=:dot, color=c, label="")
plot!(p_lhf, co2, lhf, linewidth=2, linestyle=:dot, color=c, label="")
plot!(p_zi, co2, zi, linewidth=2, linestyle=:dot, color=c, label="")
plot!(p_cf, co2, cf*100, linewidth=2, linestyle=:dot, color=c, label="")

# plot LES downsteps
co2 = [200,300,400,800,1000,1200,1300,1400];
sst = [287.6,296.8,297.9,300.9,302.0,303.7,304.7,306.2];
cf = [1.00,0.46,0.42,0.32,0.30,0.26,0.24,0.22];
zi = [1416,1340,1230,1013,949,866,826,766];
zb = [978,589,544,481,462,437,425,410];
lhf = [96.1,179.1,183.7,195.2,199.7,206.0,209.4,214.4];
dR = [82.1,9.6,9.7,7.2,6.6,5.1,4.5,3.6];
S = (lhf./dR).*((zi.-zb)./zi);

c = "royalblue"
scatter!(p_dR, co2, dR, color=c, marker=:x, markersize=ms, label="")
scatter!(p_decoup, co2, S, color=c, marker=:x, markersize=ms, label="")
scatter!(p_sst, co2, sst, color=c, marker=:x, markersize=ms, label="")
scatter!(p_lhf, co2, lhf, color=c, marker=:x, markersize=ms, label="")
scatter!(p_zi, co2, zi, color=c, marker=:x, markersize=ms, label="")
scatter!(p_cf, co2, cf*100, color=c, marker=:x, markersize=ms, label="")
plot!(p_dR, co2, dR, linewidth=2, linestyle=:dot, color=c, label="")
plot!(p_decoup, co2, S, linewidth=2, linestyle=:dot, color=c, label="")
plot!(p_sst, co2, sst, linewidth=2, linestyle=:dot, color=c, label="")
plot!(p_lhf, co2, lhf, linewidth=2, linestyle=:dot, color=c, label="")
plot!(p_zi, co2, zi, linewidth=2, linestyle=:dot, color=c, label="")
plot!(p_cf, co2, cf*100, linewidth=2, linestyle=:dot, color=c, label="")

###################

exp_path = ARGS[1];
co2u = ARGS[2];
co2u = parse.(Int64, split(replace(replace(co2u, "]"=>""), "["=>""),","))
co2d = ARGS[3];
co2d = parse.(Int64, split(replace(replace(co2d, "]"=>""), "["=>""),","))
print(co2u)
println(typeof(co2u))
# co2d = [1800, 1600, 1400, 1200, 1000, 800, 600, 400, 300]#, 200];

N = length(co2u);
zi, zb, ent = zeros(N), zeros(N), zeros(N);
cf, lwp, sst, lhf = zeros(N), zeros(N), zeros(N), zeros(N);
dR = zeros(N);

for (i, co2i) in enumerate(co2u)
    if co2i == 400
        file = "experiments/output/"*exp_path*"co2_400.jld2"
    else
        file = "experiments/output/"*exp_path*"co2_upstep_"*string(co2i)*".jld2"
    end
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
scatter!(p_dR, co2u, dR, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_decoup, co2u, S, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_sst, co2u, sst, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_lhf, co2u, lhf, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_zi, co2u, zi, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p_cf, co2u, cf*100, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
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

# legend
scatter!(p_dR, [-1], [-1], color="black", marker=:x, markersize=1, markerstrokewidth=0.1, label="LES")
scatter!(p_dR, [-1], [-1], color="black", marker=:circle, markersize=1, markerstrokewidth=0, label="Bulk model")

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
    legend=:topright, legendfontsize=12, legendfont=font(12),
    left_margin=10Plots.mm, bottom_margin=7Plots.mm, top_margin=5Plots.mm);
mkpath("experiments/figures/"*exp_path)
savefig(p, "experiments/figures/"*exp_path*"hystersis_plot_min.png")

# reset fontsizes
Plots.scalefontsizes(1/1.8)