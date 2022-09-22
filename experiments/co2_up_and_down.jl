exp_path = "20220921_EISasTropSST_m5_subLIN_ventMM2_100304_SWb100/";

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
par.SW_a = 120;
par.SW_b = 100;
par.decoup_slope = 5;
par.α_vent = 2.0e-3;
par.EIS0 = 10.0;
par.ECS = 3.0;
par.Eexport = 4.0;

# baseline
ARGS = [exp_path]; include("co2_400.jl")

# upsteps
par.stype = varSST();
# ARGS = ["200.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["300.0", exp_path]; include("co2_upstep.jl")
# ARGS = ["400.0", exp_path]; include("co2_upstep.jl")
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
# ARGS = ["200.0", exp_path]; include("co2_downstep.jl")

# plot
ARGS = [exp_path]; include("AGUplot_steadystate.jl")

# try last upstep
ARGS = ["200.0", exp_path]; include("co2_downstep.jl")
ARGS = ["100.0", exp_path]; include("co2_downstep.jl")