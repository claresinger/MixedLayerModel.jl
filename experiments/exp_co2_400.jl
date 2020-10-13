using MixedLayerModel

# create or overwrite output file
filename = "exp_co2_400.txt"
open(string("experiments/output/",filename), "w") do io 
    write(io, "")
end;

# define OHU from 400 ppm simulation
par = basic_params();
par.rtype = varRad();
u0, uf = run_mlm(par, filename);

# set OHU, increase CO2, let SST evolve and check cloud changes
par.OHU = calc_OHU(uf, par, par.stype);
par.stype = varSST();
u0, uf = run_mlm_from_init(uf, par, filename);