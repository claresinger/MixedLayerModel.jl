# test pressure <= psurf
@test pres(0.0, 290.0) == psurf
@test pres(0.0, 300.0) == psurf
@test pres(200.0, 300.0) < psurf

# test liquid water specific humidity >= 0
@test q_l(200.0, 288.0, 12e-3) > 0
@test q_l(200.0, 288.0, 8e-3) == 0.0

# test temperature decreases with altitude
z = collect(0:10:500);
T = 290.0;
qtM = 0.8 * q_sat(0.0, T);
sM = Cp * T;
@test all(diff(temp.(z, sM, qtM)) .< 0)

# test LCL calculation in all limits
zi = 900.0
u = [zi, sM, qtM, 290.0, 1.0];
zb = calc_LCL(u);
@test 0 < zb < zi

zi = 100.0;
u = [zi, sM, qtM, 290.0, 1.0];
zb = calc_LCL(u);
@test zb == zi

u = [900.0, sM, qtM*3, 290.0, 1.0];
zb = calc_LCL(u);
@test zb == 0.0

# test moist adiabat
p = upCO2();
zft = 5000;
Tsurf = 300;
T, z = moist_adiabat(Tsurf, zft, p);
@test T == sort!(T)
Tft = temp_ft(Tsurf, zft, p);
@test Tft < Tsurf

zft = 400;
T, z = moist_adiabat(Tsurf, zft, p);
@test T == sort!(T)
Tft = temp_ft(Tsurf, zft, p);
@test Tft < Tsurf

p.RHtrop0 = 1.01;
zft = 5000;
T, z = moist_adiabat(Tsurf, zft, p);
@test T == sort!(T)
Tft = temp_ft(Tsurf, zft, p);
@test Tft < Tsurf