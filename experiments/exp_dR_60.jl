using MixedLayerModel

# create or overwrite output file
filename = "exp_dR_60.txt"
open(string("experiments/output/",filename), "w") do io 
    write(io, "")
end;

# run to equilibrium in base state with ΔR = 80 W/m2
par = basic_params();
u0, uf = run_mlm(par, filename);

# set OHU, decrease cloud-top longwave cooling, let SST evolve and check cloud changes
par.OHU = calc_OHU(uf, par, par.stype);
par.ΔR = 60.0;
par.stype = varSST();
u0, uf = run_mlm_from_init(uf, par, filename);