using Test

push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using MixedLayerModel
using MixedLayerModel:psurf, Cp, L0

@testset "Thermo Tests" begin
    include("test_thermodynamics.jl")
end

@testset "Config Tests" begin
    include("test_configs.jl")
end