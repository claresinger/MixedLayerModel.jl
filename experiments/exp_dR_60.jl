# using MixedLayerModel

###
# run to equilibrium in base state with ΔR = 60 W/m2
###
filename = "exp_dR_60.txt"
open(string("experiments/output/",filename), "w") do io 
    write(io, "\n")
end;

# define OHU from baseline 80 W/m2 simulation
par = basic_params();
u0, uf = run_with_output(par, filename);
OHU = calc_surf_RAD(uf,par) - calc_SHF(uf,par) - calc_LHF(uf,par);
open(string("experiments/output/",filename), "a") do io
    write(io, string("OHU: ",OHU," (W/m^2)\n\n"))
end;

# then decrease cloud-top longwave cooling and check
par = basic_params();
par.ΔR = 60.0;
par.dSST = 1.0;
par.OHU = OHU;
u0, uf = run_with_output(par, filename);