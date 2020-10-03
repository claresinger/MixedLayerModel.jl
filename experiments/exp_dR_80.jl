include("../src/Radiation.jl")
include("../src/MLMparams.jl")
include("../src/MLMrun.jl")

###
# run to equilibrium in base state with ΔR = 80 W/m2
###

par = basic_params();

u0, uf = run_with_output(par);

OHU = calc_surf_RAD(uf,par) - calc_SHF(uf,par) - calc_LHF(uf,par);
println("OHU: ",OHU," (W/m^2)");

par = basic_params();
par.dSST = 1.0;
par.OHU = OHU;

u0, uf = run_with_output(par);