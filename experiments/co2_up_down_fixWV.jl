exp_path = "20221110_fixWVrad_co2loop_Cd8/";

using MixedLayerModel

# create parameters
par = upCO2();
par.etype = enBal();
par.fttype = co2EIS();
par.rtype = varRad();
par.stype = fixSST();
par.wvtype = wvRADOFF();
dt, tmax = 10.0, 100.0; # days

par.τCF = 2.0 # 2 days
par.Hw = 0.2; # 10 days
# adjust tunable parameters
par.α_vent = 1.22e-3;
par.Cd = 8e-4;
par.EIS0 = 8.0;
par.ECS = 1.5;
par.Eexport = 10.0;
par.SW_b = 150;

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
ARGS = ["1600.0", exp_path]; include("co2_upstep.jl")
ARGS = ["2000.0", exp_path]; include("co2_upstep.jl")
ARGS = ["2400.0", exp_path]; include("co2_upstep.jl")
ARGS = ["2800.0", exp_path]; include("co2_upstep.jl")
ARGS = ["3200.0", exp_path]; include("co2_upstep.jl")
ARGS = ["3600.0", exp_path]; include("co2_upstep.jl")

# downsteps
ARGS = ["3600.0", exp_path]; include("co2_downstep.jl")
ARGS = ["3200.0", exp_path]; include("co2_downstep.jl")
ARGS = ["2800.0", exp_path]; include("co2_downstep.jl")
ARGS = ["2400.0", exp_path]; include("co2_downstep.jl")
ARGS = ["2000.0", exp_path]; include("co2_downstep.jl")
ARGS = ["1600.0", exp_path]; include("co2_downstep.jl")
ARGS = ["1200.0", exp_path]; include("co2_downstep.jl")
ARGS = ["1000.0", exp_path]; include("co2_downstep.jl")
ARGS = ["800.0", exp_path]; include("co2_downstep.jl")
ARGS = ["600.0", exp_path]; include("co2_downstep.jl")
ARGS = ["400.0", exp_path]; include("co2_downstep.jl")
ARGS = ["200.0", exp_path]; include("co2_downstep.jl")

# plot
co2u = "[200, 400, 600, 800, 1000, 1200, 1600, 2000, 2400, 2800, 3200, 3600]";
co2d = "[3600, 3200, 2800, 2400, 2000, 1600, 1200, 1000, 800, 600, 400, 200]";
ARGS = [exp_path, co2u, co2d];
xtks = ([0, 1000, 2000, 3000])
xrange = [0,3800]
SSTrange = [285, 310]
LHFrange = [50, 250]

orig_co2u = "[200, 300, 400, 600, 800, 1000, 1200, 1400, 1600, 1800]";
orig_co2d = "[1800, 1600, 1400, 1200, 1000, 800, 600, 400, 300, 200, 100]";
orig_path = "20221109_Hw_α1.22/";
include("plot_hysteresis_noLES.jl")
