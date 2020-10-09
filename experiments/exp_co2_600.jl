# using MixedLayerModel

###
# run to equilibrium with elevated CO2 of 600 ppm
###
filename = "exp_co2_600.txt"
open(string("experiments/output/",filename), "w") do io 
    write(io, "\n")
end;
# define OHU from 400 ppm simulation
par = basic_params();
par.rtype = varRad();
par.CO2 = 400.0;
u0, uf = run(par, filename);

# set OHU, increase CO2, let SST evolve and check cloud changes
par.OHU = calc_OHU(uf, par, par.stype);
par.CO2 = 600.0;
par.stype = varSST();
u0, uf = run_from_init(uf, par, filename);