exp_path = "20221126_slope8_α122_V16/";

# create parameters
using MixedLayerModel
par = upCO2();
par.etype = enBal();
par.fttype = co2EIS();
par.rtype = varRad();
par.stype = fixSST();
dt, tmax = 10.0, 100.0; # days

par.τCF = 2.0 # 2 days
par.Hw = 0.2; # 10 days
# adjust tunable parameters
par.decoup_slope = 8;
par.α_vent = 1.22e-3;
# par.Cd = 8e-4;
par.Cd = 5e-4;
par.V = 16.0;
# par.EIS0 = 8.0;
# par.ECS = 1.5;
# par.Eexport = 10.0;
# par.SW_b = 150;

# tuned "20221101_LES_noise5pct_newprior_Nens20_Niter10"
# ϕ_final: [0.0008107894654652039, 0.001140080332952163, 5.0508769710181625, 4.692294040652614, 15.03709820113598, 118.76428605612544]
# par.Cd = 0.81e-3;
# par.α_vent = 1.12e-3;
# par.EIS0 = 5.05;
# par.ECS = 4.69;
# par.Eexport = 15.04;
# par.SW_b = 118.8;

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
ARGS = [exp_path, co2u, co2d];
include("plot_hysteresis.jl")
