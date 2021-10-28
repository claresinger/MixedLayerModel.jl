using MixedLayerModel, Documenter

# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "100"

pages = Any[
    "Home" => "index.md"
    "Mixed Layer Theory" => "theory.md"
    "ODE Solver" => "ode_solver.md"
    "CO``_2`` Perturbation Experiment" => "exp.md"
    # "Fixed SST Results" => "results_fixSST.md"
    # "Slab Ocean Results" => "results_slab.md"
    "APIs" => "library.md"
]

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    collapselevel = 1,
)

makedocs(
    sitename="MixedLayerModel.jl",
    format = format,
    clean = true,
    strict = false,
    modules = [MixedLayerModel],
    pages = pages,
)

deploydocs(
    repo = "github.com/claresinger/MixedLayerModel.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    versions = ["dev" => "dev", "v#.#.#"],
    forcepush = true,
)