exp_path = "20221028_fixSST_fixEIS_co2loop/";

# create parameters
using MixedLayerModel
par = upCO2();
par.etype = enBal();
par.fttype = fixEIS();
par.rtype = varRad();
par.stype = fixSST();
dt, tmax = 48.0, 100.0;

par.Cd = 1e-3;

# adjust tunable parameters
par.SW_b = 150;
par.decoup_slope = 8;
par.α_vent = 0.95e-3;
par.EIS0 = 10.0;
par.ECS = 3.0;
par.Eexport = 15.0;

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
co2u = "[200, 400, 800, 1200, 1600, 2000, 2400, 2800, 3200, 3600, 4000]";
co2d = "[4000, 3600, 3200, 2800, 2400, 2000, 1600, 1200, 800, 400, 200]";
ARGS = [exp_path, co2u, co2d];
xtks = ([0, 1000, 2000, 3000, 4000])
xrange = [0,4500]
SSTrange = [286, 294]
LHFrange = [0, 150]
include("plot_hysteresis_noLES.jl")
