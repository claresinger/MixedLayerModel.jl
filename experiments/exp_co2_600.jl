include("../src/Radiation.jl")
include("../src/MLMparams.jl")
include("../src/MLMrun.jl")

###
# run to equilibrium with elevated CO2 of 600 ppm
###

par = basic_params();
par.rtype = varRad();
par.CO2 = 600.0;

u0, uf = run_with_output(par);

OHU = calc_surf_RAD(uf,par) - calc_SHF(uf,par) - calc_LHF(uf,par);
println("OHU: ",OHU," (W/m^2)");

par = basic_params();
par.rtype = varRad();
par.CO2 = 600.0;
par.dSST = 1.0;
par.OHU = OHU;

u0, uf = run_with_output(par);