using MixedLayerModel
using Plots
using LaTeXStrings

plot(size=(300,200), dpi=300, legend=false, fontsize=10);

De = collect(range(0,2,length=100));
par = climatology();
CF = cloud_fraction_param.(De, Ref(par));

plot!(De, CF);
xlabel!("Decoupling, D");
ylabel!("Cloud fraction");
ylims!((0,1));

mkpath("./figures/");
savefig("./figures/cf-vs-de.png");