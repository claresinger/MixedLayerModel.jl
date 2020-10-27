using Documenter, MixedLayerModel

# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "100"

makedocs(
    sitename="MixedLayerModel.jl",
    doctest = false,
    strict = false,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",),
    modules = [Documenter, MixedLayerModel],
    clean = false,
    pages = Any[
        "Home" => "index.md"
        "Mixed Layer Theory" => "theory.md"
        "Entrainment" => "entrainment.md"
        "Radiation" => "radiation.md"
        "Thermodynamics" => "thermo.md"
        "APIs" => "library.md"
    ]
)

deploydocs(
    repo = "github.com/claresinger/MixedLayerModel.jl.git",
    target = "build"
)