push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using MixedLayerModel
using FileIO
using Plots
using NCDatasets

par = upCO2();
zft = 3e3;

for Tsurf in 300:315
    T,z = moist_adiabat(Tsurf, zft, par);
end

plot(size=(400,700), dpi=200, left_margin = 5Plots.mm);
for Ts in 295:5:315
    T,z = moist_adiabat(Ts, zft, par);
    p = pres.(z,T) / 100;
    plot!(T, p, yaxis=(:flip, :log, (750,1020)), 
        yticks=(500:100:1000, 500:100:1000), 
        markershape=:circle, markerstrokewidth=0, label=Ts)
end

ds = Dataset("experiments/LES_temp_profiles.nc")
p_LES = ds["pressure"]
T_LES = ds["temperature"]
co2_LES = ds["co2"]
println(T_LES[1,:])
for (i,c) in enumerate(co2_LES)
    if i == 1
        plot!(T_LES[:,i], p_LES, color=:black, markershape=:circle, label="LES")
    else
        plot!(T_LES[:,i], p_LES, color=:black, markershape=:circle, label="")
    end
end

xlims!(280,318)
xlabel!("Temperature (K)")
ylabel!("Pressure (hPa)")
savefig("experiments/figures/moist_adiabats.png");
