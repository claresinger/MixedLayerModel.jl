push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using MixedLayerModel
using MixedLayerModel: Cp, g, psurf
using Plots

z = 0:1500;
sM = Cp*300;
qM = 10e-3;
D = 6e-6;
ρ0 = rho.(z, temp.(z, sM, qM));
p0 = pres.(z, temp.(z, sM, qM));
sub_func = (psurf .- p0) .* (p0 ./ psurf).^2 ./ (ρ0 .* g);
w_sub = D * sub_func;
ω = D .* (psurf .- p0) .* (p0 ./ psurf).^2;

plot(size=(800,500), layout=(1,2), dpi=200, left_margin = 5Plots.mm);
plot!(w_sub * 1e3, z, xlabel="sub velocity [mm/s]", ylabel="z [m]", subplot=1, legend=false);
plot!(D .* z * 1e3, z, subplot=1);
plot!(ω * (24*60*60) / 100, p0 / 100, xlabel="ω [hPa/day]", ylabel="pressure [hPa]", subplot=2, yflip=true, legend=false);

savefig("subsidence_profile.png");