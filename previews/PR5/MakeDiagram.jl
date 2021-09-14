using Plots
using LaTeXStrings

plot(size=(1000,600), dpi=300, legend=false, axis=nothing, framestyle=:none);

plot!([1,80],[0,0],linecolor="black");
plot!([1,80], [50,50], fillrange = [100,100], color="grey", fillalpha=0.5);
plot!([15,15,5,5],[0,100,100,150],linecolor="black");
plot!([45,45,38,33],[0,100,100,150],linecolor="black");
plot!([62,62,78,62,62],[0,50,100,100,150],linecolor="black");

RH1 = 38 .+ collect(range(0,7,length=50));
RH2 = 45 .+ zeros(50);
RH3 = 28 .- collect(range(0,2,length=5));
RH = vcat(RH1..., RH2..., RH3...);

z1 = collect(range(0,100,length=100));
z2 = collect(range(100,150,length=5));
z = vcat(z1..., z2...);

plot!(RH, z, linecolor="black", linestyle=:dot);

annotate!(0,50, text(L"$z_b$"));
annotate!(0,100, text(L"$z_i$"));

annotate!(17,10, text(L"$h_M$"));
annotate!(36,10, text(L"$RH(z)$"));
annotate!(48,10, text(L"$q_{tM}$"));
annotate!(65,10, text(L"$q_l(z)$"));

annotate!(8,120, text(L"$h_+(z)$"));
annotate!(30,120, text(L"$RH_{+}$"));
annotate!(40,120, text(L"$q_{t+}(z)$"));

mkpath("./figures/");
savefig("./figures/mlm-diagram-with-ql.png");

##################

plot(size=(800,600), dpi=300, legend=false, axis=nothing, framestyle=:none);

plot!([1,50],[0,0],linecolor="black");
plot!([1,50], [50,50], fillrange = [100,100], color="grey", fillalpha=0.5);
plot!([15,15,5,5],[0,100,100,150],linecolor="black");
plot!([45,45,38,33],[0,100,100,150],linecolor="black");

RH1 = 38 .+ collect(range(0,7,length=50));
RH2 = 45 .+ zeros(50);
RH3 = 28 .- collect(range(0,2,length=5));
RH = vcat(RH1..., RH2..., RH3...);

z1 = collect(range(0,100,length=100));
z2 = collect(range(100,150,length=5));
z = vcat(z1..., z2...);

plot!(RH, z, linecolor="black", linestyle=:dot);

annotate!(0,50, text(L"$z_b$"));
annotate!(0,100, text(L"$z_i$"));

annotate!(17,10, text(L"$h_M$"));
annotate!(36,10, text(L"$RH(z)$"));
annotate!(47,10, text(L"$q_{tM}$"));

annotate!(8,120, text(L"$h_+(z)$"));
annotate!(30,120, text(L"$RH_{+}$"));
annotate!(39,120, text(L"$q_{t+}(z)$"));

mkpath("./figures/");
savefig("./figures/mlm-diagram.png");