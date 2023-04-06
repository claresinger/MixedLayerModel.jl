using MixedLayerModel
using MixedLayerModel: g, Cp
using NCDatasets
using Plots

file = "experiments/data/climatology_for_radcool.nc";
ds = Dataset(file, "r");
CFdata = vcat(ds["allsc"]...);
Ddata = vcat((σ*ds["LHF"] ./ ds["dF"])...);

par = climatology();
par.Dcrit = 0.9;
Dfit = collect(0:0.01:1.5);
CFfit = cloud_fraction_param.(Dfit, Ref(par));

plot(figsize=(8,6), dpi=300, xlabel="\$\\mathcal{D}\$, decoupling parameter", ylabel="Cloud fraction")
scatter!(Ddata, CFdata, color=:grey, ms=5, markerstrokewidth=0, label="ERA5/MODIS/CASSCAD")
plot!(Dfit, CFfit, color=:red, ls=:solid, lw=2, label="CF, Eq. 5")
savefig("D_vs_CF.png")
