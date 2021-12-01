push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using MixedLayerModel
using MixedLayerModel:psurf, Cp, L0

@test pres(0.0, 290.0) == psurf
@test pres(0.0, 300.0) == psurf
@test pres(200.0, 300.0) < psurf

@test q_v(200.0, 288.0, 12e-3) == q_sat(200.0, 288.0)
@test q_v(200.0, 288.0, 8e-3) == 8e-3

@test q_l(200.0, 288.0, 12e-3) > 0
@test q_l(200.0, 288.0, 8e-3) == 0.0

# test potential temperature increases
# and temperature decreases with altitude
z = collect(0:10:500);
T = 290.0;
qtM = 0.8 * q_sat(0.0, T);
hM = Cp * T + L0 * qtM;
@test all(diff(theta.(z, hM, qtM)) .> 0)
@test all(diff(temp.(z, hM, qtM)) .< 0)

# test LCL calculation in all limits
zi = 900.0
u = [zi, hM, qtM, 290.0, 1.0];
zb = calc_LCL(u);
@test 0 < zb < zi

zi = 100.0;
u = [zi, hM, qtM, 290.0, 1.0];
zb = calc_LCL(u);
@test zb == zi

u = [900.0, hM, qtM*3, 290.0, 1.0];
zb = calc_LCL(u);
@test zb == 0.0