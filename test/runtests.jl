using Test

@testset "Thermodynamics" begin
    include("test_thermodynamics.jl")
end

@testset "SteadyState" begin
    include("test_steadystate.jl")
end
