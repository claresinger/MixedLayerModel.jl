include("forward_model.jl") # calls whole set of CO2 experiments in MLM

# Import modules
using Distributions  # probability distributions and associated functions
using LinearAlgebra
using StatsPlots
using Plots
using Random
using JLD2

using EnsembleKalmanProcesses
using EnsembleKalmanProcesses.Observations
using EnsembleKalmanProcesses.DataContainers
using EnsembleKalmanProcesses.ParameterDistributions

const EKP = EnsembleKalmanProcesses

rng_seed = 4137
Random.seed!(rng_seed)

# LES data, truth
co2_up = [200,300,400,800,1000,1200,1300,1400,1600];
sst_up = [287.7,289.1,290.0,292.2,293.2,294.3,304.5,305.8,308.0];
co2_dn = [200,300,400,800,1000,1200,1300,1400];
sst_dn = [287.6,296.8,297.9,300.9,302.0,303.7,304.7,306.2];

###
###  Define the (true) parameters
###
# Define the parameters that we want to learn
decoup_slope_true = 10;
α_vent_true = 1.08e-3;
EIS0_true = 10;
ECS_true = 3;
Eexport_true = 15;
SW_b_true = 150;

params_true = [decoup_slope_true, α_vent_true, EIS0_true, ECS_true, Eexport_true, SW_b_true]
param_names = ["decoup_slope","α_vent", "EIS0", "ECS", "Eexport", "SW_b"]

n_param = length(param_names)
params_true = reshape(params_true, (n_param, 1))
println("ϕ_true: ", params_true)

###
###  Define the parameter priors
###
prior_decoup = constrained_gaussian("prior_decoup", 10, 1, 0.0, Inf)
prior_α = constrained_gaussian("prior_α", 1.08e-3, 0.01e-3, 0.0, Inf)
prior_EIS = constrained_gaussian("prior_EIS", 10, 1, 0.0, Inf)
prior_ECS = constrained_gaussian("prior_ECS", 3, 0.5, 0.0, Inf)
prior_Eexport = constrained_gaussian("prior_Eexport", 12, 2, 0.0, Inf)
prior_SW = constrained_gaussian("prior_SW", 100, 30, 0.0, Inf)
priors = combine_distributions([prior_decoup, prior_α, prior_EIS, prior_ECS, prior_Eexport, prior_SW])

###
###  Define the data from which we want to learn the parameters
###
# data_names = ["y0", "y1"]
data_names = ["SST"]

###
###  Generate (artificial) truth samples
###  Note: The observables y are related to the parameters θ by:
###        y = G(θ) + η
###

# Input: params: [N_params, N_ens]
# Output: gt: [N_data, N_ens]
# Dropdims of the output since the forward model is only being run with N_ens=1 
# corresponding to the truth construction
gt = dropdims(run_forward(params_true), dims = 2)

n_samples = 100
yt = zeros(length(gt), n_samples)
noise_level = 1e-2
Γy = noise_level * convert(Array, Diagonal(gt))
μ = zeros(length(gt))
# Add noise
for i in 1:n_samples
    yt[:, i] = gt .+ rand(MvNormal(μ, Γy))
end

# Construct observation object
truth = Observations.Observation(yt, Γy, data_names)
truth_sample = truth.mean

###
###  Calibrate: Ensemble Kalman Inversion
###

N_ens = 50 # number of ensemble members
N_iter = 10 # number of EKI iterations
# initial parameters: N_params x N_ens
initial_params = construct_initial_ensemble(priors, N_ens; rng_seed = rng_seed)
ϕ_init_mean = transform_unconstrained_to_constrained(priors, mean(initial_params, dims=2));
println("ϕ_init_mean: ", ϕ_init_mean)

ekiobj = EKP.EnsembleKalmanProcess(initial_params, truth_sample, truth.obs_noise_cov, Inversion())

# EKI iterations
println("EKP inversion error:")
err = zeros(N_iter)
for i in 1:N_iter
    #ϕ_i = get_ϕ_final(priors, ekiobj)
    ϕ_i = transform_unconstrained_to_constrained(priors, get_u_final(ekiobj))
    g_ens = run_ensembles(ϕ_i, N_ens)
    EKP.update_ensemble!(ekiobj, g_ens)
    err[i] = get_error(ekiobj)[end] #mean((params_true - mean(params_i,dims=2)).^2)
    println("Iteration: " * string(i) * ", Error: " * string(err[i]))
end

# u_final = mean(get_u_final(ekiobj), dims=2)
# println("u_final: ", u_final)

# ϕ_final = get_ϕ_mean_final(priors, ekiobj)
ϕ_final = transform_unconstrained_to_constrained(priors, mean(get_u_final(ekiobj), dims=2))
println("ϕ_final: ", ϕ_final)


########
# Output and save data/figures
########

# Output figure save directory
homedir = pwd()
save_directory = homedir * "/experiments/ekp/20220927_perfectmodel/"
if ~isdir(save_directory)
    mkpath(save_directory)
end

u_stored = get_u(ekiobj, return_array = false)
g_stored = get_g(ekiobj, return_array = false)
println(g_stored)
@save save_directory * "parameter_storage.jld2" u_stored
@save save_directory * "data_storage.jld2" g_stored
@save save_directory * "ekiobj.jld2" ekiobj

# @load save_directory * "ekiobj.jld2" ekiobj

# plot SST
plot(size=(600,400), layout=(1,1), dpi=200, palette = palette(:heat, N_iter+1))
CO2updn_list = [200,300,400,800,1000,1200,1300,1400,1600,1400,1300,1200,1000,800,400,300,200];
for i in 1:N_iter
    g_i = get_g(ekiobj, i)
    for j in 1:N_ens
        g_ij = g_i[:,j]
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
plot!(CO2updn_list, truth_sample, linestyle=:solid, color=:black, label="Truth")
savefig(save_directory * "SST_hysteresis_loop.png")

# plot u parameter convergence
u_init = get_u_prior(ekiobj)
for i in 1:N_iter
    plot(size=(1200, 400), layout=(1,3), dpi=200, left_margin=10Plots.mm, bottom_margin = 10Plots.mm)
    u_i = get_u(ekiobj, i)
    for pl in 1:3
        plot!(
            u_i[pl*2-1, :], 
            u_i[pl*2, :], 
            seriestype = :scatter, 
            xlims = extrema(u_init[pl*2-1, :]), 
            ylims = extrema(u_init[pl*2, :]),
            label = false,
            subplot=pl,
        )
        plot!(
            [params_true[pl*2-1]],
            xaxis = "u"*string(pl*2-1),
            yaxis = "u"*string(pl*2),
            seriestype = "vline",
            linestyle = :dash,
            linecolor = :red,
            label = false,
            title = pl==2 ? "EKI iteration = " * string(i) : "",
            subplot=pl,
        )
        plot!(
            [params_true[pl*2]], 
            seriestype = "hline", 
            linestyle = :dash, 
            linecolor = :red, 
            label = false,
            subplot=pl,
        )
    end
    figpath = joinpath(save_directory, "posterior_EKP_it_$(i).png")
    savefig(figpath)
    sleep(0.5)
end