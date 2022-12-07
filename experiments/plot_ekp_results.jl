PERFECT=false

# Import modules
using EnsembleKalmanProcesses
using Statistics
using Plots
using StatsPlots
using JLD2

include("normalize.jl")

# point to path
N_ens = 90 # number of ensemble members
N_iter = 15 # number of EKI iterations
N_params = 3; # number of calibrated parameters
NNstring = "Nens" *string(N_ens) * "_Niter" * string(N_iter)
save_directory = "experiments/ekp/20221205_LES_10pct_jumps_constraints_3params_" * NNstring * "/"

# load ekiobj data
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
err = get_error(ekiobj);
lims = (floor(minimum(err)/100)*100, ceil(maximum(err)/100)*100);
plot!(
    1:N_iter,
    err,
    xaxis="Iterations",
    yaxis="Error",
    yscale=:log10,
    yticks=(lims[1]:200:lims[2], Int.(lims[1]:200:lims[2])),
    ylims=lims,
    marker=:o,
    color=:black,
    ms=5,
    label=""
)
savefig(save_directory * "error_iterations.png")
####################

# plot prior and posterior on calibrated variables
plot(size=(800,300), layout=(1,2), dpi=300, left_margin=5Plots.mm, bottom_margin = 5Plots.mm)

yt = hcat(truth.samples...);
samples = vcat(unnormalize_data(yt[1:nd,:], "SST"), unnormalize_data(yt[nd+1:end,:], "LHF"));

plot!(CO2updn_list, mean(samples,dims=2)[1:nd], 
    yerr = std(samples,dims=2)[1:nd], markerstrokecolor=:black,
    linestyle=:solid, lw=2, color=:black, label= PERFECT ? "Truth" : "LES", 
    subplot=1, title="a)", titleloc = :left, titlefont = font(10))
plot!(CO2updn_list, mean(samples,dims=2)[nd+1:end], 
    yerr = std(samples,dims=2)[nd+1:end], markerstrokecolor=:black,
    linestyle=:solid, lw=2, color=:black, label=false, 
    subplot=2, title="b)", titleloc = :left, titlefont = font(10))

gi = get_g(ekiobj, 1)
σi = zeros(nd,2)
for i in 1:nd
    SST = unnormalize_data(gi[i,:], "SST")
    LHF = unnormalize_data(gi[i+nd,:], "LHF")
    σi[i,1] = std(collect(skipmissing(replace(SST, NaN=>missing))))
    σi[i,2] = std(collect(skipmissing(replace(LHF, NaN=>missing))))
end

gf = get_g_final(ekiobj)
σf = zeros(nd,2)
for i in 1:nd
    SST = unnormalize_data(gf[i,:], "SST")
    LHF = unnormalize_data(gf[i+nd,:], "LHF")
    σf[i,1] = std(collect(skipmissing(replace(SST, NaN=>missing))))
    σf[i,2] = std(collect(skipmissing(replace(LHF, NaN=>missing))))
end

plot!(CO2updn_list, SSTi, yerr=σi[:,1], subplot=1, linestyle=:solid, 
    color=:grey, markerstrokecolor=:grey,
    label="Prior", xaxis="CO₂ [ppmv]", yaxis="SST [K]", legend=:topleft)
plot!(CO2updn_list, LHFi, yerr=σi[:,2], subplot=2, linestyle=:solid, 
    color=:grey, markerstrokecolor=:grey,
    label=false, xaxis="CO₂ [ppmv]", yaxis="LHF [W/m²]")

plot!(CO2updn_list, SSTf, yerr=σf[:,1], subplot=1, label="Optimal",
    linestyle=:solid, color=:red3, markerstrokecolor=:red3, lw=2)
plot!(CO2updn_list, LHFf, yerr=σf[:,2], subplot=2, label=false,
    linestyle=:solid, color=:red3, markerstrokecolor=:red3, lw=2)
savefig(save_directory * "prior_posterior.png")
####################

# # plot SST
# plot(size=(600,400), layout=(1,1), dpi=200, palette = palette(:heat, N_iter+1))
# for i in 1:N_iter
#     g_i = get_g(ekiobj, i)
#     for j in 1:N_ens
#         g_ij = unnormalize_data(g_i[1:nd,j], "SST");
#         scatter!(
#             CO2updn_list, 
#             g_ij, 
#             marker=:o, 
#             ms=5, 
#             xaxis="CO2 [ppmv]",
#             yaxis="SST [K]",
#             color=i+1, 
#             label= j==1 ? "Iter. "*string(i) : false,
#         )            
#         plot!(CO2updn_list, g_ij, linestyle=:dot, color=i+1, label=false)
#     end
# end
# if PERFECT
#     plot!(CO2updn_list, unnormalize_data(truth.mean[1:nd], "SST"), 
#         linestyle=:solid, lw=3, color=:black, label="Truth", legend=:topleft)
# else
#     plot!(CO2updn_list, SSTupdn_list, 
#         linestyle=:solid, lw=3, color=:black, label="LES", legend=:topleft)
# end
# savefig(save_directory * "SST_hysteresis_loop.png")
# plot!(ylims=[280,320])
# savefig(save_directory * "SST_hysteresis_loop_ylim.png")
# ####################

# # plot LHF
# plot(size=(600,400), layout=(1,1), dpi=200, palette = palette(:heat, N_iter+1))
# for i in 1:N_iter
#     g_i = get_g(ekiobj, i)
#     for j in 1:N_ens
#         g_ij = unnormalize_data(g_i[nd+1:end,j], "LHF")
#         scatter!(
#             CO2updn_list, 
#             g_ij, 
#             marker=:o, 
#             ms=5, 
#             xaxis="CO2 [ppmv]",
#             yaxis="LHF [W/m2]",
#             color=i+1, 
#             label= j==1 ? "Iter. "*string(i) : false,
#         )            
#         plot!(CO2updn_list, g_ij, linestyle=:dot, color=i+1, label=false)
#     end
# end
# if PERFECT
#     plot!(CO2updn_list, unnormalize_data(truth.mean[nd+1:end], "LHF"), 
#         linestyle=:solid, lw=3, color=:black, label="LES", legend=:topleft)
# else
#     plot!(CO2updn_list, LHFupdn_list,
#         linestyle=:solid, lw=3, color=:black, label="LES", legend=:topleft)
# end
# savefig(save_directory * "LHF_hysteresis_loop.png")
# plot!(ylims=[50,250])
# savefig(save_directory * "LHF_hysteresis_loop_ylim.png")
# ####################

# # plot SST and LHF for each iteration
# plot(size=(800,200*N_iter), layout=(N_iter,2), dpi=200, left_margin=10Plots.mm)
# for i in 1:N_iter
#     g_i = get_g(ekiobj, i)
#     for j in 1:N_ens
#         g_ij = unnormalize_data(g_i[1:nd,j], "SST")
#         scatter!(
#             CO2updn_list,
#             g_ij,
#             marker=:o,
#             ms=5,
#             xaxis="CO2 [ppmv]",
#             yaxis="SST [K]",
#             subplot=2*i-1,
#             label=false
#         )  
#         plot!(CO2updn_list, g_ij, linestyle=:dot, subplot=2*i-1, label=false)
        
#         g_ij = unnormalize_data(g_i[nd+1:end,j], "LHF")
#         scatter!(
#             CO2updn_list,
#             g_ij,
#             marker=:o,
#             ms=5,
#             xaxis="CO2 [ppmv]",
#             yaxis="LHF [W/m2]",
#             subplot=2*i,
#             label=false
#         )            
#         plot!(CO2updn_list, g_ij, linestyle=:dot, subplot=2*i, label=false)
#     end
# end
# if PERFECT
#     for i in 1:N_iter
#         plot!(CO2updn_list, unnormalize_data(truth.mean[1:nd], "SST"), subplot=2*i-1,
#             linestyle=:solid, lw=3, color=:black, label="truth", legend=:topleft)
#         plot!(CO2updn_list, unnormalize_data(truth.mean[nd+1:end], "LHF"), subplot=2*i,
#             linestyle=:solid, lw=3, color=:black, label=false)
#     end
# else
#     for i in 1:N_iter
#         plot!(CO2updn_list, SSTupdn_list, subplot=2*i-1,
#             linestyle=:solid, lw=3, color=:black, label="LES", legend=:topleft)
#         plot!(CO2updn_list, LHFupdn_list, subplot=2*i,
#             linestyle=:solid, lw=3, color=:black, label=false)
#     end
# end
# savefig(save_directory * "hysteresis_loop_iterations.png")

# # plot!(ylims=[50,250], subplot=)
# # savefig(save_directory * "hysteresis_loop_iterations_ylim.png")
# ####################

# plot SST and LHF for each iteration
plot(size=(800,300), layout=(1,2), dpi=200, left_margin=7Plots.mm, bottom_margin=5Plots.mm)
for (ind,i) in enumerate([1,N_iter])
    g_i = get_g(ekiobj, i)
    for j in 1:N_ens
        g_ij = unnormalize_data(g_i[1:nd,j], "SST")
        scatter!(
            CO2updn_list,
            g_ij,
            marker=:o,
            ms=5,
            color=[:grey, :red3][ind],
            xaxis="CO2 [ppmv]",
            yaxis="SST [K]",
            subplot=1,
            label=j == 1 ? ["Prior", "Posterior"][ind] : false
        )  
        plot!(CO2updn_list, g_ij, linestyle=:dot, color=[:grey, :red3][ind], subplot=1, label=false)
        
        g_ij = unnormalize_data(g_i[nd+1:end,j], "LHF")
        scatter!(
            CO2updn_list,
            g_ij,
            marker=:o,
            ms=5,
            color=[:grey, :red3][ind],
            xaxis="CO2 [ppmv]",
            yaxis="LHF [W/m2]",
            subplot=2,
            label=false
        )            
        plot!(CO2updn_list, g_ij, linestyle=:dot, color=[:grey, :red3][ind], subplot=2, label=false)
    end
end
if PERFECT
    plot!(CO2updn_list, unnormalize_data(truth.mean[1:nd], "SST"), subplot=1,
        linestyle=:solid, lw=2, color=:black, label="truth", legend=:topleft)
    plot!(CO2updn_list, unnormalize_data(truth.mean[nd+1:end], "LHF"), subplot=2,
        linestyle=:solid, lw=2, color=:black, label=false)
else
    plot!(CO2updn_list, SSTupdn_list, subplot=1,
        linestyle=:solid, lw=2, color=:black, label="LES", legend=:topleft)
    plot!(CO2updn_list, LHFupdn_list, subplot=2,
        linestyle=:solid, lw=2, color=:black, label=false)
end
savefig(save_directory * "hysteresis_loop_firstlast.png")
####################

# plot ϕ parameter convergence
ϕ_list = get_ϕ(priors, ekiobj)
ϕ_all = zeros(N_params,N_ens,N_iter+1);
for (i,ϕi) in enumerate(ϕ_list)
    ϕ_all[:,:,i] = ϕi;
end
param_name = ["\$C_d \\times 10^4\$ [-]", "\$\\alpha_{\\mathrm{vent}} \\times 10^3\$ [m s⁻¹]", 
    "\$a_T\$ [K]", "\$b_T\$ [K]", "\$c_T\$ [K]", "\$b_{\\mathrm{SW}}\$ [W m⁻²]"]
scale = [10^4, 10^3, 1, 1, 1, 1]
for i in 0:N_iter
    plot(size=(1200, 400), layout=(1,3), dpi=200, left_margin=10Plots.mm, bottom_margin = 10Plots.mm)
    for pl in 1:1
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

# plot ϕ parameter convergence
plot(size=(400, 400), layout=(1,1), dpi=200, left_margin=5Plots.mm, bottom_margin = 5Plots.mm)
for (i,it) in enumerate([0, N_iter])
    for pl in 1:1
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
    savefig(save_directory * "parameters_EKP_firstlast.png")
end

####################

println(get_ϕ_mean_final(priors, ekiobj))

# plot ϕ parameter convergence
param_name = ["\$C_d \\times 10^4\$ [-]", "\$\\alpha_{\\mathrm{vent}} \\times 10^3\$ [m s⁻¹]", 
    "\$b_{\\mathrm{SW}}\$ [W m⁻²]"]
scale = [10^4, 10^3, 1]
plot(size=(900, 900), layout=(N_params,N_params), dpi=200,
    left_margin=5Plots.mm, bottom_margin = 5Plots.mm)
for (i,it) in enumerate([0, N_iter])
    for par1 in 1:N_params
        for par2 in 1:N_params
            if par1 == par2
                h1 = histogram!(
                    ϕ_all[par1, :, it+1] * scale[par1], 
                    bins = 10,
                    label = par1*par2 == 1 ? ["Initial","Final"][i] : false,
                    linecolor = [:grey, :red3][i],
                    color = [:grey, :red3][i],
                    subplot = par1 + (par2-1)*N_params,
                    xlabel = par2 == 3 ? param_name[par1] : "",
                    xlims = extrema(ϕ_all[par1, :, :]) .* [0.95,1.05] .* scale[par1],
                    yaxis = false,
                    legend = :topleft,
                )
            elseif par1 < par2
                plot!(
                    ϕ_all[par1, :, it+1] * scale[par1],
                    ϕ_all[par2, :, it+1] * scale[par2],
                    seriestype = :scatter,
                    xlims = extrema(ϕ_all[par1, :, :]) .* [0.95,1.05] .* scale[par1],
                    ylims = extrema(ϕ_all[par2, :, :]) .* [0.95,1.05] .* scale[par2], 
                    label = par1*par2 == 1 ? ["Initial","Final"][i] : false,
                    color = [:grey, :red3][i],
                    markerstrokewidth = 0,
                    subplot = par1 + (par2-1)*N_params,
                    xlabel = par2 == 3 ? param_name[par1] : "",
                    ylabel = par1 == 1 ? param_name[par2] : "",
                )
            elseif par1 > par2
                plot!(legend=false,grid=false,yaxis=false,xaxis=false,foreground_color_subplot=:white,subplot=par1+(par2-1)*N_params)
            end
        end
    end
    savefig(save_directory * "parameters_EKP_corner.png")
end

####################