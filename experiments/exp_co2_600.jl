using MixedLayerModel

###
# run to equilibrium with elevated CO2 of 600 ppm
###

# define OHU from 400 ppm simulation
par = basic_params();
par.rtype = varRad();
par.CO2 = 400.0;
u0, uf = run_with_output(par);
OHU = calc_surf_RAD(uf,par) - calc_SHF(uf,par) - calc_LHF(uf,par);
println("OHU: ",OHU," (W/m^2)");

# now increase CO2 and check
par = basic_params();
par.rtype = varRad();
par.CO2 = 600.0;
par.dSST = 1.0;
par.OHU = OHU;
u0, uf = run_with_output(par);