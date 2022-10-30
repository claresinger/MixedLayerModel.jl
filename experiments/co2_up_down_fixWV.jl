exp_path = "20221028_fixWVrad_co2loop/";

using MixedLayerModel

# create parameters
par = upCO2();
par.etype = enBal();
par.fttype = co2EIS();
par.rtype = varRad();
par.stype = fixSST();
par.wvtype = wvRADOFF();
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

par.stype = varSST();
# upsteps
ARGS = ["200.0", exp_path]; include("co2_upstep.jl")
ARGS = ["400.0", exp_path]; include("co2_upstep.jl")
ARGS = ["600.0", exp_path]; include("co2_upstep.jl")
ARGS = ["800.0", exp_path]; include("co2_upstep.jl")
ARGS = ["1000.0", exp_path]; include("co2_upstep.jl")
ARGS = ["1200.0", exp_path]; include("co2_upstep.jl")
ARGS = ["1400.0", exp_path]; include("co2_upstep.jl")
ARGS = ["1600.0", exp_path]; include("co2_upstep.jl")
ARGS = ["1800.0", exp_path]; include("co2_upstep.jl")
ARGS = ["2000.0", exp_path]; include("co2_upstep.jl")

# downsteps
ARGS = ["2000.0", exp_path]; include("co2_downstep.jl")
ARGS = ["1800.0", exp_path]; include("co2_downstep.jl")
ARGS = ["1600.0", exp_path]; include("co2_downstep.jl")
ARGS = ["1400.0", exp_path]; include("co2_downstep.jl")
ARGS = ["1200.0", exp_path]; include("co2_downstep.jl")
ARGS = ["1000.0", exp_path]; include("co2_downstep.jl")
ARGS = ["800.0", exp_path]; include("co2_downstep.jl")
ARGS = ["600.0", exp_path]; include("co2_downstep.jl")
ARGS = ["400.0", exp_path]; include("co2_downstep.jl")
ARGS = ["200.0", exp_path]; include("co2_downstep.jl")

# plot
co2u = "[200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]";
co2d = "[2000, 1800, 1600, 1400, 1200, 1000, 800, 600, 400, 200]";
ARGS = [exp_path, co2u, co2d];
xtks = ([0, 500, 1000, 1500, 2000])
xrange = [0,2200]
SSTrange = [286, 310]
LHFrange = [50, 220]
include("plot_hysteresis_noLES.jl")
