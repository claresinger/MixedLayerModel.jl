using Distributed
addprocs(1; exeflags = "--project=experiments/")

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

###
# PERFECT MODEL SETTING
###

###
###  Define the (true) parameters
###
# Define the parameters that we want to learn
# decoup_slope_true = 8;
Cd_true = 0.8e-3;
α_vent_true = 1.08e-3;
EIS0_true = 10;
ECS_true = 3;
Eexport_true = 15;
SW_b_true = 150;

# params_true = [decoup_slope_true, α_vent_true, EIS0_true, ECS_true, Eexport_true, SW_b_true]
params_true = [Cd_true, α_vent_true, EIS0_true, ECS_true, Eexport_true, SW_b_true]

n_param = length(params_true)
params_true = reshape(params_true, (n_param, 1))
println("ϕ_true: ", params_true)

###
###  Define the parameter priors
###
# prior_decoup = constrained_gaussian("prior_decoup", 10, 2, 0.0, Inf)
# prior_α = constrained_gaussian("prior_α", 1.08e-3, 0.02e-3, 0.0, Inf)
# prior_EIS = constrained_gaussian("prior_EIS", 10, 2, 0.0, Inf)
# prior_ECS = constrained_gaussian("prior_ECS", 3, 1, 0.0, Inf)
# prior_Eexport = constrained_gaussian("prior_Eexport", 12, 4, 0.0, Inf)
# prior_SW = constrained_gaussian("prior_SW", 120, 50, 0.0, Inf)
# prior_decoup = constrained_gaussian("prior_decoup", 10, 5, 0.0, Inf)
# prior_α = constrained_gaussian("prior_α", 1e-3, 1e-3, 0.0, Inf)
# prior_EIS = constrained_gaussian("prior_EIS", 10, 5, 0.0, Inf)
# prior_ECS = constrained_gaussian("prior_ECS", 3, 2, 0.0, Inf)
# prior_Eexport = constrained_gaussian("prior_Eexport", 12, 5, 0.0, Inf)
# prior_SW = constrained_gaussian("prior_SW", 150, 100, 0.0, Inf)
# priors = combine_distributions([prior_decoup, prior_α, prior_EIS, prior_ECS, prior_Eexport, prior_SW])

prior_Cd = constrained_gaussian("prior_Cd", 0.8e-3, 0.1e-3, 0.0, Inf)
prior_α = constrained_gaussian("prior_α", 1.0e-3, 0.1e-3, 0.0, Inf)
prior_EIS = constrained_gaussian("prior_EIS", 10, 1, 0.0, Inf)
prior_ECS = constrained_gaussian("prior_ECS", 3, 1, 0.0, Inf)
prior_Eexport = constrained_gaussian("prior_Eexport", 15, 1, 0.0, Inf)
prior_SW = constrained_gaussian("prior_SW", 120, 10, 0.0, Inf)
priors = combine_distributions([prior_Cd, prior_α, prior_EIS, prior_ECS, prior_Eexport, prior_SW])

# Input: params: [N_params, N_ens]
# Output: gt: [N_data, N_ens]
# Dropdims of the output since the forward model is only being run with N_ens=1 
# corresponding to the truth construction
# gt = dropdims(GModel.run_forward(params_true), dims = 2)

# n_samples = 100
# yt = zeros(length(gt), n_samples)
# noise_level = 1e-2
# Γy = noise_level * convert(Array, Diagonal(abs.(gt)))
# μy = zeros(length(gt))
# # Add noise
# for i in 1:n_samples
#     yt[:, i] = gt .+ rand(MvNormal(μy, Γy))
# end


# Construct truth
# gt = dropdims(GModel.run_forward(params_true), dims = 2);
Gt = dropdims(GModel.run_forward(params_true), dims = 2);
nd = Int(length(Gt) / 2);
gt = vcat(GModel.unnormalize_data(Gt[1:nd], "SST"), GModel.unnormalize_data(Gt[nd+1:end], "LHF"));

# Construct normalized observations yt and error Γy
n_samples = 100;
noise_level = 1e-2;
s = noise_level * randn((length(gt), n_samples)) .+ 1
samples = gt .* s
yt = vcat(GModel.normalize_data(samples[1:nd,:], "SST"), GModel.normalize_data(samples[nd+1:end,:], "LHF"));

# Construct observation object
data_names = ["SST","LHF"]
truth = Observations.Observation(yt, data_names)
truth_sample = truth.mean

###
###  Calibrate: Ensemble Kalman Inversion
###

N_ens = 5 # number of ensemble members
N_iter = 3 # number of EKI iterations
initial_params = construct_initial_ensemble(priors, N_ens; rng_seed = rng_seed) # initial parameters: N_params x N_ens
ϕ_init_mean = transform_unconstrained_to_constrained(priors, mean(initial_params, dims=2));
println("ϕ_init_mean: ", ϕ_init_mean)

ekiobj = EKP.EnsembleKalmanProcess(
    initial_params, 
    truth.mean, 
    truth.obs_noise_cov, 
    Inversion(),
    failure_handler_method = SampleSuccGauss(),
 )

# EKI iterations
@time begin
    println("EKP inversion error:")
    err = zeros(N_iter)
    for i in 1:N_iter
        ϕ_i = get_ϕ_final(priors, ekiobj)
        @time g_ens = GModel.run_ensembles(ϕ_i, N_ens)
        EKP.update_ensemble!(ekiobj, g_ens)
        err[i] = get_error(ekiobj)[end]
        println("Iteration: " * string(i) * ", Error: " * string(err[i]))
        println(get_ϕ_mean_final(priors, ekiobj))
    end
end

ϕ_final = get_ϕ_mean_final(priors, ekiobj)
println("ϕ_final: ", ϕ_final)


# ########
# # Output and save data/figures
# ########

# Output figure save directory
homedir = pwd()
NNstring = "Nens" *string(N_ens) * "_Niter" * string(N_iter)
save_directory = homedir * "/experiments/ekp/20221101_perf_normalizeY_" * NNstring * "/"
if ~isdir(save_directory)
    mkpath(save_directory)
end

u_stored = get_u(ekiobj, return_array = false)
g_stored = get_g(ekiobj, return_array = false)

CO2updn_list = [200,300,400,800,1000,1200,1300,1400,1600,1400,1300,1200,1000,800,400,300,200];
nd = length(CO2updn_list);
ϕi = get_ϕ_mean(priors, ekiobj, 1)
gi = GModel.run_ensembles(ϕi, 1)
SSTi = GModel.unnormalize_data(gi[1:nd,1], "SST")
LHFi = GModel.unnormalize_data(gi[nd+1:end,1], "LHF")
ϕf = get_ϕ_mean_final(priors, ekiobj)
gf = GModel.run_ensembles(ϕf, 1)
SSTf = GModel.unnormalize_data(gf[1:nd,1], "SST")
LHFf = GModel.unnormalize_data(gf[nd+1:end,1], "LHF")

@save save_directory * "parameter_storage.jld2" u_stored
@save save_directory * "data_storage.jld2" g_stored
@save save_directory * "ekiobj.jld2" ekiobj
@save save_directory * "prior_posterior.jld2" SSTi LHFi SSTf LHFf
@save save_directory * "truth_sample.jld2" truth_sample
@save save_directory * "params_true.jld2" params_true
