using Test

@testset "Thermo Tests" begin
    include("test_thermodynamics.jl")
end

@testset "Config Tests" begin
    include("test_configs.jl")
end