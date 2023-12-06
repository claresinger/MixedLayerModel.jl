using MixedLayerModel
using Plots
using LaTeXStrings

plot(size=(300,200), dpi=300, legend=false, fontsize=10);

S = collect(range(0,2,length=100));
m = 10; # tunable parameter for the slope of the CF nonlinearity
S_crit = 0.7; # tunable parameter for the halfway point of CF decrease
CF = 1 .- 0.8 ./ (1 .+ exp.(-m .* (S .- S_crit)));

plot!(S, CF);
xlabel!("Stability, S");
ylabel!("Cloud fraction");
ylims!((0,1));

mkpath("./figures/");
savefig("./figures/cf-vs-s.png");