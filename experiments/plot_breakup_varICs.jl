using FileIO
using Plots
using Statistics
using MixedLayerModel
using MixedLayerModel: Rd, Rv, L0, T0, Cp, δ, ϵ, μ


exp_path = ARGS[1];
co2u = ARGS[2];
expN = ARGS[3];

ENV["GKSwstype"]="nul"
Plots.scalefontsizes(1.2)

c = "crimson"
xtks = ([400, 800, 1200, 1600, 2000])
xrange = [300,2100]
SSTrange = [286, 312]
LHFrange = [50, 250]
p_dR = scatter([], [], label="", ylabel="ΔR [W/m²]",
    annotation = (xrange[1]*1.1, 100*1.1, text("a)", fontsize=12)),
    xticks=xtks, xlim=xrange, ylim=[0,100])
p_decoup = scatter([], [], label="", xlabel="CO₂ [ppmv]", ylabel="Decoupling, \$\\mathcal{D}\$",
    yscale=:log10, yticks=([0.1,1,10], ["0.1","1","10"]), ylim=[0.1,10],
    xticks=xtks, xlim=xrange)
p_sst = scatter([], [], label="", ylabel="SST [K]",
    annotation = (xrange[1]*1.1, (SSTrange[2] - SSTrange[1])*1.1 + SSTrange[1], text("b)", fontsize=12)),
    xticks=xtks, xlim=xrange, ylim=SSTrange)
p_lhf = scatter([], [], label="", xlabel="CO₂ [ppmv]", ylabel="LHF [W/m²]",
    annotation = (xrange[1]*1.1, (LHFrange[2] - LHFrange[1])*1.1 + LHFrange[1], text("c)", fontsize=12)),
    xticks=xtks, xlim=xrange, ylim=LHFrange)
p_zi = scatter([], [], label="", ylabel="zᵢ [m]",
    xticks=xtks, xlim=xrange, ylim=[0, 1000])
p_cf = scatter([], [], label="", xlabel="CO₂ [ppmv]", ylabel="CF [%]",
    annotation = (xrange[1]*1.1, 110*1.1, text("d)", fontsize=12)),
    xticks=xtks, xlim=xrange, ylim=[0,110])


N = length(co2u);
all_dR = zeros(N, expN)
all_S = zeros(N, expN)
all_sst = zeros(N, expN)
all_lhf = zeros(N, expN)
all_zi = zeros(N, expN)
all_cf = zeros(N, expN)

# loop over each initial condition variation
for expi in 1:expN
    zi, zb, ent = zeros(N), zeros(N), zeros(N);
    cf, lwp, sst, lhf = zeros(N), zeros(N), zeros(N), zeros(N);
    dR = zeros(N);

    for (i, co2i) in enumerate(co2u)
        file = "experiments/output/"*exp_path*string(expi)*"co2_upstep_"*string(co2i)*".jld2"
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

    # Store results for this experiment
    all_dR[:, expi] = dR
    all_S[:, expi] = S
    all_sst[:, expi] = sst
    all_lhf[:, expi] = lhf
    all_zi[:, expi] = zi
    all_cf[:, expi] = cf

    # ms = 2
    # plot!(p_dR, co2u, dR, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
    # plot!(p_decoup, co2u, S, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
    # plot!(p_sst, co2u, sst, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
    # plot!(p_lhf, co2u, lhf, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
    # plot!(p_zi, co2u, zi, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
    # plot!(p_cf, co2u, cf*100, color=c, marker=:circle, markersize=ms, markerstrokewidth=0, label="")

    lw = 0.5
    plot!(p_dR, co2u, dR, color=c, linewidth=lw, label="")
    plot!(p_decoup, co2u, S, color=c, linewidth=lw, label="")
    plot!(p_sst, co2u, sst, color=c, linewidth=lw, label="")
    plot!(p_lhf, co2u, lhf, color=c, linewidth=lw, label="")
    plot!(p_zi, co2u, zi, color=c, linewidth=lw, label="")
    plot!(p_cf, co2u, cf*100, color=c, linewidth=lw, label="")
end

# Compute and plot means
dR_mean = mean(all_dR, dims=2)[:]
S_mean = mean(all_S, dims=2)[:]
sst_mean = mean(all_sst, dims=2)[:]
lhf_mean = mean(all_lhf, dims=2)[:]
zi_mean = mean(all_zi, dims=2)[:]
cf_mean = mean(all_cf, dims=2)[:]

# Plot the averages with a different style (thicker, darker line)
lw_mean = 3
c_mean = "darkred"
ms = 3
plot!(p_dR, co2u, dR_mean, color=c_mean, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
plot!(p_decoup, co2u, S_mean, color=c_mean, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
plot!(p_sst, co2u, sst_mean, color=c_mean, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
plot!(p_lhf, co2u, lhf_mean, color=c_mean, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
plot!(p_zi, co2u, zi_mean, color=c_mean, marker=:circle, markersize=ms, markerstrokewidth=0, label="")
plot!(p_cf, co2u, cf_mean*100, color=c_mean, marker=:circle, markersize=ms, markerstrokewidth=0, label="")

plot!(p_dR, co2u, dR_mean, color=c_mean, linewidth=lw_mean, label="")
plot!(p_decoup, co2u, S_mean, color=c_mean, linewidth=lw_mean, label="")
plot!(p_sst, co2u, sst_mean, color=c_mean, linewidth=lw_mean, label="")
plot!(p_lhf, co2u, lhf_mean, color=c_mean, linewidth=lw_mean, label="")
plot!(p_zi, co2u, zi_mean, color=c_mean, linewidth=lw_mean, label="")
plot!(p_cf, co2u, cf_mean*100, color=c_mean, linewidth=lw_mean, label="")


# save plot
p_final = plot(
    p_dR, p_sst, p_lhf, p_cf,
    layout = (2,2),
    link = :x,
    size = (900,500),
    dpi = 300,
    legend = :topright,
    legendfontsize = 10,
    legendfont = Plots.font(10),
    left_margin   = 7 * Plots.mm,
    right_margin  = 2 * Plots.mm,
    bottom_margin = 6 * Plots.mm,
    top_margin    = 6 * Plots.mm,
);
mkpath("experiments/figures/"*exp_path)
savefig(p_final, "experiments/figures/"*exp_path*"breakup_plot.png")
savefig(p_final, "experiments/figures/"*exp_path*"breakup_plot.pdf")

# reset fontsizes
Plots.scalefontsizes(1/1.2)
