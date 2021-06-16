using MixedLayerModel, Documenter

# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "100"

pages = Any[
    "Home" => "index.md"
    "Mixed Layer Theory" => "theory.md"
    "Running an experiment" => "exp.md"
    "Results" => "results.md"
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
    target = "build"
)