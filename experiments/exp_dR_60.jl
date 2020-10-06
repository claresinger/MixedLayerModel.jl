using MixedLayerModel

###
# run to equilibrium in base state with ΔR = 80 W/m2
###

# define OHU from baseline 80 W/m2 simulation
par = basic_params();
u0, uf = run_with_output(par);
OHU = calc_surf_RAD(uf,par) - calc_SHF(uf,par) - calc_LHF(uf,par);
println("OHU: ",OHU," (W/m^2)");

# then decrease cloud-top longwave cooling and check
par = basic_params();
par.ΔR = 60.0;
par.dSST = 1.0;
par.OHU = OHU;
u0, uf = run_with_output(par);