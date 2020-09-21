include("MLMparams.jl")
include("MLMrun.jl")

###
# run to equilibrium with elevated CO2 of 600 ppm
###

par = interact_surf_params();
par.CO2 = 600.0;
par.dSST = 0.0;
u0, uf = run_with_output(par);

OHU = calc_surf_RAD(uf,par) - calc_SHF(uf,par) - calc_LHF(uf,par);
println("OHU: ",OHU," (W/m^2)");

par = interact_surf_params();
par.CO2 = 600.0;
par.dSST = 1.0;
par.OHU = OHU;
println(atmos_emissivity(par))
u0, uf = run_with_output(par);