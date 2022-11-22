PERFECT=false

# Import modules
using Distributed
using EnsembleKalmanProcesses
using Statistics
using Plots
using StatsPlots
using JLD2

include("forward_model.jl") # calls whole set of CO2 experiments in MLM

CO2updn_list = [200,300,400,800,1000,1200,1300,1400,1600,1400,1300,1200,1000,800,400,300,200];
SSTupdn_list = [287.7,289.1,290.0,292.2,293.2,294.3,304.5,305.8,308.0,306.2,304.7,303.7,302.0,300.9,297.9,296.8,287.6];
LHFupdn_list = [97.8,103.9,107.1,112.8,115.3,120.5,208.7,213.5,220.9,214.4,209.4,206.0,199.7,195.2,183.7,179.1,96.1];
nd = length(CO2updn_list);

# load ekiobj data
homedir = pwd()
N_ens = 5 # number of ensemble members
N_iter = 2 # number of EKI iterations
NNstring = "Nens" *string(N_ens) * "_Niter" * string(N_iter)
save_directory = homedir * "/experiments/ekp/20221121_LES_noise5pct_" * NNstring * "/"
@load save_directory * "ekiobj.jld2" ekiobj
@load save_directory * "priors.jld2" priors
@load save_directory * "prior_posterior.jld2" SSTi LHFi SSTf LHFf
@load save_directory * "truth.jld2" truth
if PERFECT
    @load save_directory * "params_true.jld2" params_true
end
####################

# plot error vs. iterations
plot(size=(600,300), layout=(1,1), dpi=400)
plot!(
    1:N_iter,
    get_error(ekiobj),
    xaxis="Iterations",
    yaxis="Error",
    yscale=:log10,
    marker=:o,
    ms=5,
    label=""
)
savefig(save_directory * "error_iterations.png")
####################

# plot prior and posterior on calibrated variables
plot(size=(800,400), layout=(1,2), dpi=300, left_margin=5Plots.mm, bottom_margin = 5Plots.mm)

yt = hcat(truth.samples...);
samples = vcat(GModel.unnormalize_data(yt[1:nd,:], "SST"), GModel.unnormalize_data(yt[nd+1:end,:], "LHF"));

plot!(CO2updn_list, mean(samples,dims=2)[1:nd], 
    yerr = std(samples,dims=2)[1:nd], markerstrokecolor=:black,
    linestyle=:solid, lw=3, color=:black, label="Truth", subplot=1)
plot!(CO2updn_list, mean(samples,dims=2)[nd+1:end], 
    yerr = std(samples,dims=2)[nd+1:end], markerstrokecolor=:black,
    linestyle=:solid, lw=3, color=:black, label=false, subplot=2)

gi = get_g(ekiobj, 1)
σi = zeros(nd,2)
for i in 1:nd
    SST = GModel.unnormalize_data(gi[i,:], "SST")
    LHF = GModel.unnormalize_data(gi[i,:], "LHF")
    σi[i,1] = std(collect(skipmissing(replace(SST, NaN=>missing))))
    σi[i,2] = std(collect(skipmissing(replace(LHF, NaN=>missing))))
end

gf = get_g_final(ekiobj)
σf = zeros(nd,2)
for i in 1:nd
    SST = GModel.unnormalize_data(gf[i,:], "SST")
    LHF = GModel.unnormalize_data(gf[i,:], "LHF")
    σf[i,1] = std(collect(skipmissing(replace(SST, NaN=>missing))))
    σf[i,2] = std(collect(skipmissing(replace(LHF, NaN=>missing))))
end

plot!(CO2updn_list, SSTi, yerr=σi[:,1], subplot=1, linestyle=:solid, 
    color=:pink, markerstrokecolor=:pink,
    label="Prior", xaxis="CO₂ [ppmv]", yaxis="SST [K]", legend=:topleft)
plot!(CO2updn_list, LHFi, yerr=σi[:,2], subplot=2, linestyle=:solid, 
    color=:pink, markerstrokecolor=:pink,
    label=false, xaxis="CO₂ [ppmv]", yaxis="LHF [W/m²]")

plot!(CO2updn_list, SSTf, yerr=σf[:,1], subplot=1, label="Posterior",
    linestyle=:solid, color=:green, markerstrokecolor=:green)
plot!(CO2updn_list, LHFf, yerr=σf[:,2], subplot=2, label=false,
    linestyle=:solid, color=:green, markerstrokecolor=:green)
savefig(save_directory * "prior_posterior.png")
####################

# plot SST
plot(size=(600,400), layout=(1,1), dpi=200, palette = palette(:heat, N_iter+1))
for i in 1:N_iter
    g_i = get_g(ekiobj, i)
    for j in 1:N_ens
        g_ij = GModel.unnormalize_data(g_i[1:nd,j], "SST");
        scatter!(
            CO2updn_list, 
            g_ij, 
            marker=:o, 
            ms=5, 
            xaxis="CO2 [ppmv]",
            yaxis="SST [K]",
            color=i+1, 
            label= j==1 ? "Iter. "*string(i) : false,
        )            
        plot!(CO2updn_list, g_ij, linestyle=:dot, color=i+1, label=false)
    end
end
if PERFECT
    plot!(CO2updn_list, GModel.unnormalize_data(truth.mean[1:nd], "SST"), 
        linestyle=:solid, lw=3, color=:black, label="Truth", legend=:topleft)
else
    plot!(CO2updn_list, SSTupdn_list, 
        linestyle=:solid, lw=3, color=:black, label="LES", legend=:topleft)
end
savefig(save_directory * "SST_hysteresis_loop.png")
plot!(ylims=[280,320])
savefig(save_directory * "SST_hysteresis_loop_ylim.png")
####################

# plot LHF
plot(size=(600,400), layout=(1,1), dpi=200, palette = palette(:heat, N_iter+1))
for i in 1:N_iter
    g_i = get_g(ekiobj, i)
    for j in 1:N_ens
        g_ij = GModel.unnormalize_data(g_i[nd+1:end,j], "LHF")
        scatter!(
            CO2updn_list, 
            g_ij, 
            marker=:o, 
            ms=5, 
            xaxis="CO2 [ppmv]",
            yaxis="LHF [W/m2]",
            color=i+1, 
            label= j==1 ? "Iter. "*string(i) : false,
        )            
        plot!(CO2updn_list, g_ij, linestyle=:dot, color=i+1, label=false)
    end
end
if PERFECT
    plot!(CO2updn_list, GModel.unnormalize_data(truth.mean[nd+1:end], "LHF"), 
        linestyle=:solid, lw=3, color=:black, label="LES", legend=:topleft)
else
    plot!(CO2updn_list, LHFupdn_list,
        linestyle=:solid, lw=3, color=:black, label="LES", legend=:topleft)
end
savefig(save_directory * "LHF_hysteresis_loop.png")
plot!(ylims=[50,250])
savefig(save_directory * "LHF_hysteresis_loop_ylim.png")
####################

# plot SST and LHF for each iteration
plot(size=(800,200*N_iter), layout=(N_iter,2), dpi=200, left_margin=10Plots.mm)
for i in 1:N_iter
    g_i = get_g(ekiobj, i)
    for j in 1:N_ens
        g_ij = GModel.unnormalize_data(g_i[1:nd,j], "SST")
        scatter!(
            CO2updn_list, 
            g_ij, 
            marker=:o, 
            ms=5, 
            xaxis="CO2 [ppmv]",
            yaxis="SST [K]",
            subplot=2*i-1,
            label=false
        )  
        plot!(CO2updn_list, g_ij, linestyle=:dot, subplot=2*i-1, label=false)
        
        g_ij = GModel.unnormalize_data(g_i[nd+1:end,j], "LHF")
        scatter!(
            CO2updn_list, 
            g_ij, 
            marker=:o, 
            ms=5, 
            xaxis="CO2 [ppmv]",
            yaxis="LHF [W/m2]",
            subplot=2*i,
            label=false
        )            
        plot!(CO2updn_list, g_ij, linestyle=:dot, subplot=2*i, label=false)
    end
end
if PERFECT
    for i in 1:N_iter
        plot!(CO2updn_list, GModel.unnormalize_data(truth.mean[1:nd], "SST"), subplot=2*i-1,
            linestyle=:solid, lw=3, color=:black, label="truth", legend=:topleft)
        plot!(CO2updn_list, GModel.unnormalize_data(truth.mean[nd+1:end], "LHF"), subplot=2*i,
            linestyle=:solid, lw=3, color=:black, label=false)
    end
else
    for i in 1:N_iter
        plot!(CO2updn_list, SSTupdn_list, subplot=2*i-1,
            linestyle=:solid, lw=3, color=:black, label="LES", legend=:topleft)
        plot!(CO2updn_list, LHFupdn_list, subplot=2*i,
            linestyle=:solid, lw=3, color=:black, label=false)
    end
end
savefig(save_directory * "hysteresis_loop_iterations.png")
# plot!(ylims=[50,250], subplot=)
# savefig(save_directory * "hysteresis_loop_iterations_ylim.png")
####################

# plot u parameter convergence
# u_init = get_u_prior(ekiobj)
ϕ_all = get_ϕ(priors, ekiobj)
ϕ_all = reshape(vcat(ϕ_all...),(6,N_ens,N_iter+1));
param_name = ["\$C_d \\times 10^4\$ [-]", "\$\\alpha_{\\mathrm{vent}} \\times 10^3\$ [m s⁻¹]", 
    "\$a_T\$ [K]", "\$b_T\$ [K]", "\$c_T\$ [K]", "\$b_{\\mathrm{SW}}\$ [W m⁻²]"]
scale = [10^4, 10^3, 1, 1, 1, 1]
for i in 0:N_iter
    plot(size=(1200, 400), layout=(1,3), dpi=200, left_margin=10Plots.mm, bottom_margin = 10Plots.mm)
    for pl in 1:3
        plot!(
            ϕ_all[pl*2, :, i+1] * scale[pl*2],
            ϕ_all[pl*2-1, :, i+1] * scale[pl*2-1],
            seriestype = :scatter, 
            xlims = extrema(ϕ_all[pl*2, :, :]) .* scale[pl*2],
            ylims = extrema(ϕ_all[pl*2-1, :, :]) .* scale[pl*2-1], 
            label = false,
            subplot=pl,
            xlabel = param_name[pl*2],
            ylabel = param_name[pl*2-1],
            title = pl==2 ? "EKI iteration = " * string(i) : "",
        )
        if PERFECT
            plot!(
                [params_true[pl*2-1]],
                xaxis = "u"*string(pl*2),
                yaxis = "u"*string(pl*2-1),
                seriestype = "hline",
                linestyle = :dash,
                linecolor = :red,
                label = false,
                subplot=pl,
            )
            plot!(
                [params_true[pl*2]], 
                seriestype = "vline", 
                linestyle = :dash, 
                linecolor = :red, 
                label = false,
                subplot=pl,
            )
        end
    end
    figpath = joinpath(save_directory, "parameters_EKP_it_$(i).png")
    savefig(figpath)
end
####################

# plot u parameter convergence
plot(size=(1200, 400), layout=(1,3), dpi=200, left_margin=10Plots.mm, bottom_margin = 10Plots.mm)
for (i,it) in enumerate([0, N_iter])
    for pl in 1:3
        plot!(
            ϕ_all[pl*2, :, it+1] * scale[pl*2],
            ϕ_all[pl*2-1, :, it+1] * scale[pl*2-1],
            seriestype = :scatter, 
            xlims = extrema(ϕ_all[pl*2, :, :]) .* [0.95,1.05] .* scale[pl*2],
            ylims = extrema(ϕ_all[pl*2-1, :, :]) .* [0.95,1.05] .* scale[pl*2-1], 
            label = pl == 1 ? ["Initial","Final"][i] : false,
            color = [:black, :red][i],
            subplot=pl,
            xlabel = param_name[pl*2],
            ylabel = param_name[pl*2-1],
        )
    end
    figpath = joinpath(save_directory, "parameters_EKP_firstlast.png")
    savefig(figpath)
end
####################