using Distributed
addprocs(8; exeflags = "--project=experiments/")

include("forward_model.jl") # calls whole set of CO2 experiments in MLM
include("normalize.jl")

# Import modules
using Distributions  # probability distributions and associated functions
using Statistics
using LinearAlgebra
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
# construct truth (normalized)
gt = vcat(normalize_data(SSTupdn_list, "SST"), normalize_data(LHFupdn_list, "LHF"));

# construct observations
n_samples = 100
yt = zeros(length(gt), n_samples)
noise_level = 0.1 # 10% on diagonals
# Γy = noise_level .* convert(Array, Diagonal(abs.(gt)))
Γy = noise_level .* mean(abs.(gt)) .* convert(Array, Diagonal(ones(nd*2)));
ijump = [9, 17, 26, 34]
for i in ijump
    Γy[i,i] /= 20;
end
# scatter(vcat(CO2updn_list, CO2updn_list), diag(Γy))
# savefig("out.png")

μy = zeros(length(gt))
# Add noise
for i in 1:n_samples
    yt[:, i] = gt .+ rand(MvNormal(μy, Γy))
end

# Construct observation object
data_names = ["SST","LHF"]
truth = Observations.Observation(yt, data_names)

###
###  Define the parameter priors
###
# prior_Cd = constrained_gaussian("prior_Cd", 8e-4, 1e-4, 0, 5e-3)
# prior_α = constrained_gaussian("prior_α", 1.22e-3, 0.05e-3, 0, 5e-3)
# prior_EIS = constrained_gaussian("prior_EIS", 10, 1, 0.0, Inf)
# prior_ECS = constrained_gaussian("prior_ECS", 1.5, 0.5, 0.0, Inf)
# prior_Eexport = constrained_gaussian("prior_Eexport", 10, 1, 0.0, Inf)
# prior_SW = constrained_gaussian("prior_SW", 150, 10, 0.0, Inf)
# priors = combine_distributions([prior_Cd, prior_α, prior_EIS, prior_ECS, prior_Eexport, prior_SW])
prior_Cd = ParameterDistribution(Parameterized(Normal(8e-4, 0.5e-4)), bounded(1e-4, 1e-2), "prior_Cd")
prior_fluxα = ParameterDistribution(Parameterized(Normal(0.45, 0.05)), bounded(0,1), "prior_α")
prior_SW = ParameterDistribution(Parameterized(Normal(140, 10)), bounded(50,250), "prior_SW")
priors = combine_distributions([prior_Cd, prior_fluxα, prior_SW])

###
###  Calibrate: Ensemble Kalman Inversion
###
N_ens = 50 # number of ensemble members
N_iter = 3 # number of EKI iterations
println(N_ens, " ", N_iter)
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
println()
err = zeros(N_iter)
for i in 1:N_iter
    ϕ_i = get_ϕ_final(priors, ekiobj)
    @time g_ens = GModel.run_ensembles(ϕ_i, N_ens)
    EKP.update_ensemble!(ekiobj, g_ens)
    err[i] = get_error(ekiobj)[end]
    println("Iteration: " * string(i) * ", Error: " * string(err[i]))
    println(get_ϕ_mean_final(priors, ekiobj))
end

ϕ_final = get_ϕ_mean_final(priors, ekiobj)
println("ϕ_final: ", ϕ_final)


########
# Output and save data/figures
########

# Output figure save directory
homedir = pwd()
NNstring = "Nens" *string(N_ens) * "_Niter" * string(N_iter)
save_directory = homedir * "/experiments/ekp/20231024_LES_10pct_jumps_constraints_3params_" * NNstring * "/"
if ~isdir(save_directory)
    mkpath(save_directory)
end

u_stored = get_u(ekiobj, return_array = false)
ϕ_stored = get_ϕ(priors, ekiobj)
g_stored = get_g(ekiobj, return_array = false)

ϕi = get_ϕ_mean(priors, ekiobj, 1)
gi = GModel.run_ensembles(ϕi, 1)
SSTi = unnormalize_data(gi[1:nd,1], "SST")
LHFi = unnormalize_data(gi[nd+1:end,1], "LHF")
ϕf = get_ϕ_mean_final(priors, ekiobj)
gf = GModel.run_ensembles(ϕf, 1)
SSTf = unnormalize_data(gf[1:nd,1], "SST")
LHFf = unnormalize_data(gf[nd+1:end,1], "LHF")

@save save_directory * "raw_parameter_storage.jld2" u_stored
@save save_directory * "real_parameter_storage.jld2" ϕ_stored
@save save_directory * "data_storage.jld2" g_stored
@save save_directory * "ekiobj.jld2" ekiobj
@save save_directory * "prior_posterior.jld2" SSTi LHFi SSTf LHFf
@save save_directory * "truth.jld2" truth
@save save_directory * "priors.jld2" priors
