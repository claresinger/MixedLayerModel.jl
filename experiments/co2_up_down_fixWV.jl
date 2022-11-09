exp_path = "20221108_fixWVrad_co2loop/";

using MixedLayerModel

# create parameters
par = upCO2();
par.etype = enBal();
par.fttype = co2EIS();
par.rtype = varRad();
par.stype = fixSST();
par.wvtype = wvRADOFF();
dt, tmax = 24.0*10, 100.0;

# adjust tunable parameters
par.Cd = 8e-4;
par.α_vent = 1e-3;
par.EIS0 = 8.0;
par.ECS = 3.0;
par.Eexport = 15.0;
par.SW_b = 160;

# # baseline
# ARGS = [exp_path]; include("co2_400.jl")

# par.stype = varSST();
# # upsteps
# ARGS = ["200.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["400.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["600.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["800.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["1000.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["1200.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["1400.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["1600.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["1800.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["2000.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["2200.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["2400.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["2600.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["2800.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["3000.0", exp_path]; include("co2_upstep.jl")

# # downsteps
# ARGS = ["3000.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["2800.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["2600.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["2400.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["2200.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["2000.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["1800.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["1600.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["1400.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["1200.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["1000.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["800.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["600.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["400.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["200.0", exp_path]; include("co2_downstep.jl")

# plot
co2u = "[200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000]";
co2d = "[3000, 2800, 2600, 2400, 2200, 2000, 1800, 1600, 1400, 1200, 1000, 800, 600, 400, 200]";
ARGS = [exp_path, co2u, co2d];
xtks = ([0, 1000, 2000, 3000])
xrange = [0,3200]
SSTrange = [286, 310]
LHFrange = [50, 220]

orig_co2u = "[200, 300, 400, 600, 800, 1000, 1200, 1400, 1600, 1800]";
orig_co2d = "[1800, 1600, 1400, 1200, 1000, 800, 600, 400, 300, 200, 100]";
orig_path = "20221102_prior/";
include("plot_hysteresis_noLES.jl")
