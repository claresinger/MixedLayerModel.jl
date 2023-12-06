using MixedLayerModel
using Plots
using LaTeXStrings

plot(size=(300,200), dpi=300, legend=false, fontsize=10);

S = collect(range(0,2,length=100));
CF = cloud_fraction_S.(S);

plot!(S, CF);
xlabel!("Stability, S");
ylabel!("Cloud fraction");
ylims!((0,1));

mkpath("./figures/");
savefig("./figures/cf-vs-s.png");