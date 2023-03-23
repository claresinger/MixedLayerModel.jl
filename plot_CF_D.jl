using MixedLayerModel
using MixedLayerModel: g, Cp
using NCDatasets
using Plots

file = "experiments/data/climatology_for_radcool.nc";
ds = Dataset(file, "r");
CFdata = vcat(ds["allsc"]...);
Ddata = vcat((σ*ds["LHF"] ./ ds["dF"])...);

println(σ_of_T(700, 295-(g/Cp)*700))

# α(z,T) = (1 - ϵ − δ*ϵ)/(1 + γ(z,T));
# Ddata = vcat((α.(zb, Tzb) .* (ds["LHF"] .- γ.(zb, Tzb) .* ds["SHF"]) ./ ds["dF"])...);
# Ddata = vcat((α.(zb, Tzb) .* (ds["LHF"] .- γ.(zb, Tzb) .* ds["SHF"]) ./ ds["dF"] .* (700/1200))...);
# χ = 0.15
# Ddata = vcat(( (α.(zb, Tzb) .* (ds["LHF"] .- γ.(zb, Tzb) .* ds["SHF"]) .- (1-χ)*ds["dF"]) ./ (χ*ds["dF"]) .* (700/1200))...);

par = climatology();
par.Dcrit = 1.0;
Dfit = collect(0:0.01:1.5);
CFfit = cloud_fraction_param.(Dfit, Ref(par));

plot(figsize=(8,6), dpi=300, xlabel="\$\\mathcal{D}\$, decoupling parameter", ylabel="Cloud fraction")
scatter!(Ddata, CFdata, color=:grey, ms=5, markerstrokewidth=0, label="ERA5/MODIS/CASSCAD")
plot!(Dfit, CFfit, color=:red, ls=:solid, lw=2, label="CF, Eq. 5")
savefig("D_vs_CF.png")

# add noise to D to reintroduce clouds
# fix sigma to constant (T = 295K)
# make Dcrit some order 1 quantity
