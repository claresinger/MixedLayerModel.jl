exp_path = "20221110_fixSST_fixEIS_co2loop_Cd8/";

# create parameters
using MixedLayerModel
par = upCO2();
par.etype = enBal();
par.fttype = fixEIS(); # no tropical warming from fixEIS()
par.rtype = varRad();
par.stype = fixSST(); # no subtropical warming from fixSST()
dt, tmax = 24.0, 100.0; # days

par.τCF = 2.0 # 2 days
par.Hw = 0.2; # 10 days
# adjust tunable parameters
par.α_vent = 1.22e-3;
par.Cd = 8e-4;
par.EIS0 = 8.0;
par.SW_b = 150;

# # baseline
# ARGS = [exp_path]; include("co2_400.jl")

# # upsteps
# ARGS = ["200.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["400.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["800.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["1200.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["2000.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["2800.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["3600.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["4400.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["5200.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["6000.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["6800.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["7600.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["8000.0", exp_path]; include("co2_upstep.jl")

# # downsteps
# ARGS = ["8000.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["7600.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["6800.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["6000.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["5200.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["4400.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["3600.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["2800.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["2000.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["1200.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["800.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["400.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["200.0", exp_path]; include("co2_downstep.jl")

# plot
co2u = "[200, 400, 800, 1200, 2000, 2800, 3600, 4400, 5200, 6000, 6800, 7600, 8000]";
co2d = "[8000, 7600, 6800, 6000, 5200, 4400, 3600, 2800, 2000, 1200, 800, 400, 200]";
ARGS = [exp_path, co2u, co2d];
xtks = ([0, 2000, 4000, 6000, 8000])
xrange = [0,8200]
SSTrange = [285, 310]
LHFrange = [50, 250]

orig_co2u = "[200, 300, 400, 600, 800, 1000, 1200, 1400, 1600, 1800]";
orig_co2d = "[1800, 1600, 1400, 1200, 1000, 800, 600, 400, 300, 200, 100]";
orig_path = "20221109_Hw_α1.22/";
include("plot_hysteresis_noLES.jl")
