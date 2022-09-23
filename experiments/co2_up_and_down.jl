exp_path = "20220923_wrongtropsst_m10_ventMM108_100315_SWb150/";

# create parameters
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using MixedLayerModel
par = upCO2();
par.etype = enBal();
par.fttype = co2EIS();
par.rtype = varRad();
par.stype = fixSST();
dt, tmax = 48.0, 50.0;

# adjust tunable parameters
par.SW_b = 150;
par.decoup_slope = 10;
par.α_vent = 1.08e-3;
par.EIS0 = 10.0;
par.ECS = 3.0;
par.Eexport = 15.0;

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
ARGS = [exp_path, co2u, co2d]; include("AGUplot_steadystate.jl")
include("plot_hysteresis.jl")

# # try last upstep
# ARGS = ["300.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["200.0", exp_path]; include("co2_downstep.jl")
# ARGS = ["100.0", exp_path]; include("co2_downstep.jl")
# co2u = "[400, 600, 800, 1000, 1200, 1400, 1600, 1800]";
# co2d = "[1800, 1600, 1400, 1200, 1000, 800, 600, 400, 300, 200, 100]";
# ARGS = [exp_path, co2u, co2d]; include("AGUplot_steadystate.jl")
# include("plot_hysteresis.jl")