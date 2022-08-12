# baseline
include("co2_400.jl")

# upsteps
# ARGS = ["200.0"]; include("co2_upstep.jl")
# ARGS = ["300.0"]; include("co2_upstep.jl")
# ARGS = ["400.0"]; include("co2_upstep.jl")
ARGS = ["600.0"]; include("co2_upstep.jl")
ARGS = ["800.0"]; include("co2_upstep.jl")
ARGS = ["1000.0"]; include("co2_upstep.jl")
ARGS = ["1200.0"]; include("co2_upstep.jl")
ARGS = ["1400.0"]; include("co2_upstep.jl")
# ARGS = ["1500.0"]; include("co2_upstep.jl")
ARGS = ["1600.0"]; include("co2_upstep.jl")
ARGS = ["1800.0"]; include("co2_upstep.jl")
ARGS = ["2000.0"]; include("co2_upstep.jl")

# downsteps
# ARGS = ["1600.0"]; include("co2_downstep.jl")
# # ARGS = ["1500.0"]; include("co2_downstep.jl")
# ARGS = ["1400.0"]; include("co2_downstep.jl")
# ARGS = ["1200.0"]; include("co2_downstep.jl")
# ARGS = ["1000.0"]; include("co2_downstep.jl")
# ARGS = ["800.0"]; include("co2_downstep.jl")
# ARGS = ["600.0"]; include("co2_downstep.jl")
# ARGS = ["400.0"]; include("co2_downstep.jl")
# # ARGS = ["300.0"]; include("co2_downstep.jl")
# ARGS = ["200.0"]; include("co2_downstep.jl")
