using MixedLayerModel

###
# run to equilibrium in base state with 400 ppm CO2
###
filename = "exp_co2_400.txt"
open(string("experiments/output/",filename), "w") do io 
    write(io, "\n")
end;

# define OHU from 400 ppm simulation
par = basic_params();
par.rtype = varRad();
u0, uf = run_with_output(par, filename);
OHU = calc_surf_RAD(uf,par) - calc_SHF(uf,par) - calc_LHF(uf,par);
open(string("experiments/output/",filename), "a") do io
    write(io, string("OHU: ",OHU," (W/m^2)\n\n"))
end;

# now run with evolving SST and check that it is consistent
par = basic_params();
par.rtype = varRad();
par.dSST = 1.0;
par.OHU = OHU;
u0, uf = run_with_output(par, filename);