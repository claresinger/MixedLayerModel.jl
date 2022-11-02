exp_path = "20221101_Cd0.8e-3_α1.08e-3_EIS10_Ex15_ECS3_SW180/";

# create parameters
using MixedLayerModel
par = upCO2();
par.etype = enBal();
par.fttype = co2EIS();
par.rtype = varRad();
par.stype = fixSST();
dt, tmax = 10*24.0, 100.0;

# params
par.decoup_slope = 8;
par.D = 6e-6;

# adjust tunable parameters
par.Cd = 0.8e-3;
par.α_vent = 1.08e-3;
par.EIS0 = 10.0;
par.ECS = 3.0;
par.Eexport = 15.0;
par.SW_b = 180;

# baseline
ARGS = [exp_path]; include("co2_400.jl")

# upsteps
par.stype = varSST();
ARGS = ["200.0", exp_path]; include("co2_upstep.jl")
ARGS = ["300.0", exp_path]; include("co2_upstep.jl")
ARGS = ["400.0", exp_path]; include("co2_upstep.jl")
ARGS = ["600.0", exp_path]; include("co2_upstep.jl")
ARGS = ["800.0", exp_path]; include("co2_upstep.jl")
ARGS = ["1000.0", exp_path]; include("co2_upstep.jl")
ARGS = ["1200.0", exp_path]; include("co2_upstep.jl")
ARGS = ["1400.0", exp_path]; include("co2_upstep.jl")
ARGS = ["1600.0", exp_path]; include("co2_upstep.jl")
ARGS = ["1800.0", exp_path]; include("co2_upstep.jl")

# downsteps
ARGS = ["1800.0", exp_path]; include("co2_downstep.jl")
ARGS = ["1600.0", exp_path]; include("co2_downstep.jl")
ARGS = ["1400.0", exp_path]; include("co2_downstep.jl")
ARGS = ["1200.0", exp_path]; include("co2_downstep.jl")
ARGS = ["1000.0", exp_path]; include("co2_downstep.jl")
ARGS = ["800.0", exp_path]; include("co2_downstep.jl")
ARGS = ["600.0", exp_path]; include("co2_downstep.jl")
ARGS = ["400.0", exp_path]; include("co2_downstep.jl")
ARGS = ["300.0", exp_path]; include("co2_downstep.jl")
ARGS = ["200.0", exp_path]; include("co2_downstep.jl")

# plot
co2u = "[200, 300, 400, 600, 800, 1000, 1200, 1400, 1600, 1800]";
co2d = "[1800, 1600, 1400, 1200, 1000, 800, 600, 400, 300, 200]";
# co2u = "[400, 600, 800, 1000, 1200, 1400, 1600, 1800]";
# co2d = "[1800, 1600, 1400, 1200, 1000, 800, 600, 400]";
ARGS = [exp_path, co2u, co2d]; #include("AGUplot_steadystate.jl")
include("plot_hysteresis.jl")
