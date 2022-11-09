exp_path = "20221109_fixSST_co2loop/";

# create parameters
using MixedLayerModel
par = upCO2();
par.etype = enBal();
par.fttype = co2EIS();
par.rtype = varRad();
par.stype = fixSST();
dt, tmax = 24.0*10.0, 100.0;

# adjust tunable parameters
par.Cd = 8e-4;
par.α_vent = 1e-3;
par.EIS0 = 8.0;
par.ECS = 3.0;
par.Eexport = 15.0;
par.SW_b = 160;

# baseline
ARGS = [exp_path]; include("co2_400.jl")

# upsteps
ARGS = ["200.0", exp_path]; include("co2_upstep.jl")
ARGS = ["400.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["600.0", exp_path]; include("co2_upstep.jl")
ARGS = ["800.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["1000.0", exp_path]; include("co2_upstep.jl")
ARGS = ["1200.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["1400.0", exp_path]; include("co2_upstep.jl")
ARGS = ["1600.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["1800.0", exp_path]; include("co2_upstep.jl")
ARGS = ["2000.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["2200.0", exp_path]; include("co2_upstep.jl")
ARGS = ["2400.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["2600.0", exp_path]; include("co2_upstep.jl")
ARGS = ["2800.0", exp_path]; include("co2_upstep.jl")
ARGS = ["3200.0", exp_path]; include("co2_upstep.jl")
ARGS = ["3600.0", exp_path]; include("co2_upstep.jl")
ARGS = ["4000.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["4400.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["4800.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["5200.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["5600.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["6000.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["6400.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["6800.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["7200.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["7600.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["8000.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["8400.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["8800.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["9200.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["9600.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["10000.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["20000.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["30000.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["40000.0", exp_path]; include("co2_upstep.jl")


# downsteps
ARGS = ["4000.0", exp_path]; include("co2_downstep.jl")
ARGS = ["3600.0", exp_path]; include("co2_downstep.jl")
ARGS = ["3200.0", exp_path]; include("co2_downstep.jl")
ARGS = ["2800.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["2600.0", exp_path]; include("co2_downstep.jl")
ARGS = ["2400.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["2200.0", exp_path]; include("co2_downstep.jl")
ARGS = ["2000.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["1800.0", exp_path]; include("co2_downstep.jl")
ARGS = ["1600.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["1400.0", exp_path]; include("co2_downstep.jl")
ARGS = ["1200.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["1000.0", exp_path]; include("co2_downstep.jl")
ARGS = ["800.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["600.0", exp_path]; include("co2_downstep.jl")
ARGS = ["400.0", exp_path]; include("co2_downstep.jl")
ARGS = ["200.0", exp_path]; include("co2_downstep.jl")

# plot
co2u = "[200, 400, 800, 1200, 1600, 2000, 2400, 2800, 3200, 3600, 4000]";#, 4400, 4800, 5200, 6000, 7200, 8000, 8800, 9200, 10000]";
co2d = "[4000, 3600, 3200, 2800, 2400, 2000, 1600, 1200, 800, 400, 200]";
ARGS = [exp_path, co2u, co2d];
xtks = ([0, 1000, 2000, 3000, 4000])
xrange = [0,4000]
SSTrange = [286, 294]
LHFrange = [0, 150]

orig_co2u = "[200, 300, 400, 600, 800, 1000, 1200, 1400, 1600, 1800]";
orig_co2d = "[1800, 1600, 1400, 1200, 1000, 800, 600, 400, 300, 200, 100]";
orig_path = "20221102_prior/";
include("plot_hysteresis_noLES.jl")
