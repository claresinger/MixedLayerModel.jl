using Distributed
addprocs(1; exeflags = "--project=experiments/")

include("forward_model.jl") # calls whole set of CO2 experiments in MLM

# Import modules
using Distributions  # probability distributions and associated functions
using Statistics
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
# REAL LES LEARNING
###
# LES data, truth
CO2updn_list = [200,300,400,800,1000,1200,1300,1400,1600,1400,1300,1200,1000,800,400,300,200];
SSTupdn_list = [287.7,289.1,290.0,292.2,293.2,294.3,304.5,305.8,308.0,306.2,304.7,303.7,302.0,300.9,297.9,296.8,287.6];
LHFupdn_list = [97.8,103.9,107.1,112.8,115.3,120.5,208.7,213.5,220.9,214.4,209.4,206.0,199.7,195.2,183.7,179.1,96.1];
ziupdn_list = [1442,1349,1266,1078,1011,972,834,781,703,766,826,866,949,1013,1230,1340,1416];
nd = length(CO2updn_list);

# gt = vcat(GModel.normalize_data(SSTupdn_list, "SST"), GModel.normalize_data(LHFupdn_list, "LHF"));

# n_samples = 100
# yt = zeros(length(gt), n_samples)
# noise_level = 1e-2 .* ones(length(gt))
# # jumps = [6,7,16,17];
# # noise_level[jumps] .-= 0.8e-2;
# Γy = noise_level .* convert(Array, Diagonal(abs.(gt)))
# μy = zeros(length(gt))
# # Add noise
# for i in 1:n_samples
#     yt[:, i] = gt .+ rand(MvNormal(μy, Γy))
# end

# Construct normalized observations yt and error Γy
gt = vcat(SSTupdn_list, LHFupdn_list);
n_samples = 100;
samples = zeros(length(gt), n_samples);
noise = 1e-2 .* ones(length(gt));
for i in 1:n_samples
    samples[:, i] = gt .+ rand(MvNormal(zeros(length(gt)), noise .* gt));
end
yt = vcat(GModel.normalize_data(samples[1:nd,:], "SST"), GModel.normalize_data(samples[nd+1:end,:], "LHF"));

# Construct observation object
data_names = ["SST","LHF"]
truth = Observations.Observation(yt, data_names)

###
###  Define the parameter priors
###
# prior_decoup = constrained_gaussian("prior_decoup", 8, 2, 0.0, Inf)
prior_Cd = constrained_gaussian("prior_Cd", 0.8e-3, 0.1e-3, 0.0, Inf)
prior_α = constrained_gaussian("prior_α", 1.08e-3, 0.1e-3, 0.0, Inf)
prior_EIS = constrained_gaussian("prior_EIS", 10, 0.1, 0.0, Inf)
prior_ECS = constrained_gaussian("prior_ECS", 3, 0.1, 0.0, Inf)
prior_Eexport = constrained_gaussian("prior_Eexport", 15, 0.1, 0.0, Inf)
prior_SW = constrained_gaussian("prior_SW", 150, 10, 0.0, Inf)
priors = combine_distributions([prior_Cd, prior_α, prior_EIS, prior_ECS, prior_Eexport, prior_SW])

###
###  Calibrate: Ensemble Kalman Inversion
###

N_ens = 20 # number of ensemble members
N_iter = 5 # number of EKI iterations
# initial parameters: N_params x N_ens
initial_params = construct_initial_ensemble(priors, N_ens; rng_seed = rng_seed)
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
    end
end

ϕ_final = get_ϕ_mean_final(priors, ekiobj)
println("ϕ_final: ", ϕ_final)


########
# Output and save data/figures
########

# Output figure save directory
homedir = pwd()
NNstring = "Nens" *string(N_ens) * "_Niter" * string(N_iter)
save_directory = homedir * "/experiments/ekp/20221101_LES_narrowprior_uniformerror_" * NNstring * "/"
if ~isdir(save_directory)
    mkpath(save_directory)
end

u_stored = get_u(ekiobj, return_array = false)
g_stored = get_g(ekiobj, return_array = false)

ϕi = get_ϕ_mean(priors, ekiobj, 1)
gi = GModel.run_ensembles(ϕi, 1)
SSTi = GModel.unnormalize_data(gi[1:nd,1], "SST")
LHFi = GModel.unnormalize_data(gi[1:nd,1], "LHF")
ϕf = get_ϕ_mean_final(priors, ekiobj)
gf = GModel.run_ensembles(ϕf, 1)
SSTf = GModel.unnormalize_data(gf[1:nd,1], "SST")
LHFf = GModel.unnormalize_data(gf[1:nd,1], "LHF")

@save save_directory * "parameter_storage.jld2" u_stored
@save save_directory * "data_storage.jld2" g_stored
@save save_directory * "ekiobj.jld2" ekiobj
@save save_directory * "prior_posterior.jld2" SSTi LHFi SSTf LHFf
