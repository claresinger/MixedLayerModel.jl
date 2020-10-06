using MixedLayerModel
include("../src/Definitions.jl")

@test pres(0.0, 290.0) == psurf
@test pres(0.0, 300.0) == psurf
@test pres(200.0, 300.0) < psurf

@test q_v(200.0, 288.0, 12e-3) == q_sat(200.0, 288.0)
@test q_v(200.0, 288.0, 8e-3) == 8e-3

@test q_l(200.0, 288.0, 12e-3) > 0
@test q_l(200.0, 288.0, 8e-3) == 0.0

z = collect(0:10:500);
qtM = 0.8 * q_sat(0.0, 290.0)
hM = Cp * (288.0) + L0 * qtM
@test all(theta(z, hM, qtM) .> temp(z, hM, qtM))