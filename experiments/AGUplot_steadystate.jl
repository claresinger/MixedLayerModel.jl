push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using FileIO
using Plots
using MixedLayerModel
using MixedLayerModel: Rd, Rv, L0, T0, Cp, δ, ϵ, μ

ENV["GKSwstype"]="nul"
Plots.scalefontsizes(2)

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
p1 = scatter(co2, dR, color=c, marker=:x, markersize=ms, label="", ylabel="ΔR [W/m\$^2\$]",
                xticks=([400, 800, 1200, 1600]), xlim=[100,1700], ylim=[0,100])
p2 = scatter(co2, S, color=c, marker=:x, markersize=ms, label="", ylabel="Stability, \$S\$",
                yscale=:log10, yticks=([0.2,0.5,2,5,20], ["0.2","0.5","2","5","20"]),
                ylim=[0.2,40], xticks=([400, 800, 1200, 1600]), xlim=[100,1700])
p3 = scatter(co2, cf*100, color=c, marker=:x, markersize=ms, label="", xlabel="CO\$_2\$ [ppmv]", ylabel="CF [%]",
                xticks=([400, 800, 1200, 1600]), xlim=[100,1700], ylim=[0,110])
p4 = scatter(co2, sst, color=c, marker=:x, markersize=ms, label="", xlabel="CO\$_2\$ [ppmv]", ylabel="SST [K]",
                xticks=([400, 800, 1200, 1600]), xlim=[100,1700], ylim=[285, 312])
plot!(p1, co2, dR, linewidth=2, linestyle=:dot, color=c, label="")
plot!(p2, co2, S, linewidth=2, linestyle=:dot, color=c, label="")
plot!(p3, co2, cf*100, linewidth=2, linestyle=:dot, color=c, label="")
plot!(p4, co2, sst, linewidth=2, linestyle=:dot, color=c, label="")

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
scatter!(p1, co2, dR, color=c, marker=:x, markersize=ms, label="")
scatter!(p2, co2, S, color=c, marker=:x, markersize=ms, label="")
scatter!(p3, co2, cf*100, color=c, marker=:x, markersize=ms, label="")
scatter!(p4, co2, sst, color=c, marker=:x, markersize=ms, label="")
plot!(p1, co2, dR, linewidth=2, linestyle=:dot, color=c, label="")
plot!(p2, co2, S, linewidth=2, linestyle=:dot, color=c, label="")
plot!(p3, co2, cf*100, linewidth=2, linestyle=:dot, color=c, label="")
plot!(p4, co2, sst, linewidth=2, linestyle=:dot, color=c, label="")

###################

# exp_path = "Tt_01/"
# co2 = [200, 400, 600, 800, 1000, 1100, 1200];
# exp_path = "Tt02_m5/"
# co2 = [400, 600, 800, 1000, 1100, 1200, 1300];
# exp_path = "Tt02_m8/"
# co2 = [400, 600, 700, 800, 900, 1000];
exp_path = "Tt02_m6/"
co2 = [400, 600, 800, 900, 1000, 1100];
N = length(co2);

zi, zb, ent = zeros(N), zeros(N), zeros(N);
cf, lwp, sst, lhf = zeros(N), zeros(N), zeros(N), zeros(N);
dR = zeros(N);

for (i, co2i) in enumerate(co2)
    if co2i == 400
        file = "experiments/output/"*exp_path*"co2_400.jld2"
    else
        file = "experiments/output/"*exp_path*"co2_upstep_"*string(co2i)*".jld2"
    end
    dat = load(file);
    uf = dat["uf"];
    par = dat["p"];
    zii, hM, qM, ssti, cfi = uf;
    zbi = dat["zb"];
    zi[i], zb[i], sst[i], cf[i] = zii, zbi, ssti, cfi;
    lhf[i], ent[i], dR[i] = dat["LHF"], dat["we"]*1e3, dat["ΔR"];
    lwp[i] = incloud_LWP(uf, zb[i]) * 1e3;
end
S = (lhf./dR).*((zi.-zb)./zi);

ms = 8
c = "crimson"
scatter!(p1, co2, dR, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p2, co2, S, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p3, co2, cf*100, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p4, co2, sst, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
plot!(p1, co2, dR, color=c, linewidth=2, label="")
plot!(p2, co2, S, color=c, linewidth=2, label="")
plot!(p3, co2, cf*100, color=c, linewidth=2, label="")
plot!(p4, co2, sst, color=c, linewidth=2, label="")

# co2 = [1200, 1000, 800, 700, 600, 400];
# co2 = [900, 800, 600, 400, 200, 100];
co2 = [1000, 800, 600, 400, 300, 200];
N = length(co2);

zi, zb, ent = zeros(N), zeros(N), zeros(N);
cf, lwp, sst, lhf = zeros(N), zeros(N), zeros(N), zeros(N);
dR = zeros(N);

for (i, co2i) in enumerate(co2)
    file = "experiments/output/"*exp_path*"co2_downstep_"*string(co2i)*".jld2"
    dat = load(file);
    uf = dat["uf"];
    par = dat["p"];
    zii, hM, qM, ssti, cfi = uf;
    zbi = dat["zb"];
    zi[i], zb[i], sst[i], cf[i] = zii, zbi, ssti, cfi;
    lhf[i], ent[i], dR[i] = dat["LHF"], dat["we"]*1e3, dat["ΔR"];
    lwp[i] = incloud_LWP(uf, zb[i]) * 1e3;
end
S = (lhf./dR).*((zi.-zb)./zi);

c = "royalblue"
scatter!(p1, co2, dR, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p2, co2, S, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p3, co2, cf*100, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
scatter!(p4, co2, sst, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
plot!(p1, co2, dR, color=c, linewidth=2, label="")
plot!(p2, co2, S, color=c, linewidth=2, label="")
plot!(p3, co2, cf*100, color=c, linewidth=2, label="")
plot!(p4, co2, sst, color=c, linewidth=2, label="")

# legend
scatter!(p1, [-1], [-1], color="black", marker=:x, markersize=1, markerstrokewidth=0, label="LES")
scatter!(p1, [-1], [-1], color="black", marker=:circle, markersize=1, markerstrokewidth=0, label="MLM")

# save plot
p = plot(p1,p2,p3,p4, layout=(2,2), 
    link=:x, size=(1000,650), dpi=300,
    legend=:topright, legendfontsize=12, legendfont=font(12),
    left_margin=10Plots.mm, bottom_margin=5Plots.mm, top_margin=5Plots.mm);
mkpath("experiments/figures/"*exp_path)
savefig(p, "experiments/figures/"*exp_path*"AGU-steady-state.png")
Plots.scalefontsizes(1/2)