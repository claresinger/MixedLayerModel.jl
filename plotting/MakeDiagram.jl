using PyPlot
using LaTeXStrings

fig = figure(figsize=(8,6));
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
rcParams["font.size"] = 15;

plot([0,80],[0,0],"k-");
fill_between([0,80],[50,50],[100,100],color="grey",alpha=0.5);
plot([20,20,10,10],[0,100,102,150],"k-");
plot([65,65,50,45],[0,100,102,150],"k-");

RH1 = 58 .+ collect(range(0,7,length=50));
RH2 = 65 .+ zeros(50);
RH3 = 35 .- collect(range(0,2,length=5));
RH = vcat(RH1..., RH2..., RH3...);

z1 = collect(range(0,100,length=100));
z2 = collect(range(102,150,length=5));
z = vcat(z1..., z2...);

plot(RH,z,"k:");

text(-4,50,L"$z_b$");
text(-4,100,L"$z_i$");

text(66,10,L"$q_{tM}$");
text(21,10,L"$h_M$");
text(51,10,L"RH(z)");

text(11,120,L"$h^{ft}(z)$");
text(49,120,L"$q_t^{ft}(z)$");
text(35,120,L"RH$^{ft}$");

axis("off");
box(false);

tight_layout();
savefig("plotting/figures/mlm-diagram.png",dpi=300);