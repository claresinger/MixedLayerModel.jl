# test albedo and emissivity is in bounds
LWP = 0:0.1:100
@test all(cloud_albedo.(LWP) .>= 0)
@test all(cloud_albedo.(LWP) .<= 1)

@test all(cloud_emissivity.(LWP) .>= 0)
@test all(cloud_emissivity.(LWP) .<= 1)

# test downwelling emission temp is reasonable
par = upCO2();
qM = 10e-3;
u = [1000.0, Cp*290, qM, 290, 1.0];
zb = calc_LCL(u);
LWP = incloud_LWP(u, zb);
T = ΔTa(u, p, LWP, p.wvtype);
@test -200 < T <= 0

# and gives same answer when qft = p.qft_rad
qft = qM + qjump(u, p, LWP, p.fttype);
par.wvtype = wvRADOFF();
par.qft_rad = qft;
@test T == ΔTa(u, p, LWP, p.wvtype);