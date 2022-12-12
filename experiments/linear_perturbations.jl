using MixedLayerModel
using MixedLayerModel: Cp, g
using Plots
using NCDatasets

include("mlm_solve_funcs.jl")

path = "experiments/figures/20221211_linear_perturb/";
mkpath(path);

N = 11
_p2 = collect(range(98,102,N)) ./ 100;
_p10 = collect(range(90,110,N)) ./ 100;
_p20 = collect(range(80,120,N)) ./ 100;
_p50 = collect(range(50,150,N)) ./ 100;
Xname_arr = ["SST", "V", "D", "RH", "EIS", "CO2"];
Xlabel_arr = ["SST [K]", "V [m/s]", "D [10⁻⁶ s⁻¹]", "RH₊ [%]", "EIS [K]", "CO₂ [ppm]"];
perturb = [_p2, _p10, _p10, _p20, _p20, _p50];
Nvar = length(Xname_arr);
print = false;

# perturb = [[1],[1],[1],[1],[1],[1]];
# Nvar = 1;
# print = true;

X_arr = zeros(length(Xname_arr), N, 2);
CF_arr = zeros(length(Xname_arr), N, 2);
LWP_arr = zeros(length(Xname_arr), N, 2);

for (k,regime) in enumerate(["Sc","Cu"])
    println(regime)
    if regime == "Sc"
        init = 1
        SST0 = 290;
        D0 = 5e-6;
        V0 = 8.0;
        RHft0 = 0.25;
        EIS0 = 12.0;
        CO20 = 400.0;
    end
    if regime == "Cu"
        init = 0
        SST0 = 295;
        D0 = 5e-6;
        V0 = 8.0;
        RHft0 = 0.25;
        EIS0 = 6.0;
        CO20 = 400.0;
    end
    X0_arr = [SST0, V0, D0, RHft0, EIS0, CO20];

    for j in 1:Nvar
        X = X0_arr[j]*perturb[j];
        # println(X)

        # set up parameters
        local par = climatology();
        par.rtype = varRad();
        par.stype = fixSST();
        par.ftype = varFlux();
        par.fttype = fixEIS();
        par.etype = enBal();
        
        # adjust tunable parameters
        par.decoup_slope = 8;
        par.α_vent = 1.69e-3;
        par.Cd = 6e-4; #7.9e-4;

        for (i,Xi) in enumerate(X)            
            par.SST0 = (j == 1) ? Xi : SST0 # (K)
            par.V = (j == 2) ? Xi : V0 # m/s
            par.D = (j == 3) ? Xi : D0 # (1/s)
            par.RHft = (j == 4) ? Xi : RHft0 # (-)
            par.EIS0 = (j == 5) ? Xi : EIS0 # (K)
            par.CO2 = (j == 6) ? Xi : CO20 # (ppm)

            # run to equilibrium
            dt = 3600.0*24.0*5.0;
            tmax = 3600.0*24.0*100.0;
            u0, sol = run_mlm(par, init=init, dt=dt, tspan=(0.0,tmax), quiet=true);
            uf = sol.u[end];
            if print
                println(uf)
                println(calc_LCL(uf))
                println(incloud_LWP(uf, calc_LCL(uf)) * 1e3)
            end

            # save
            X_arr[j,i,k] = Xi;
            CF_arr[j,i,k] = uf[5];
            LWP_arr[j,i,k] = incloud_LWP(uf, calc_LCL(uf));
        end
    end
end

# create output netcdf file
isfile(path*"linear_perturbs.nc") ? rm(path*"linear_perturbs.nc") : "no file"
ds_save = Dataset(path*"linear_perturbs.nc","c")
# Define the dimension index of size N
defDim(ds_save, "index", N)
defVar(ds_save, "index", 1:N, ("index",))
# Define the dimension var of size length(Xname_arr)
defDim(ds_save, "var", length(Xname_arr))
defVar(ds_save, "var", Xname_arr, ("var",))
# Define the dimension of base of size 2
defDim(ds_save, "base", 2)
defVar(ds_save, "base", ["Sc","Cu"], ("base",))
# Define the variables
v = defVar(ds_save, "Xval", X_arr, ("var","index","base"))
v = defVar(ds_save, "CF", CF_arr, ("var","index","base"))
v = defVar(ds_save, "inLWP", LWP_arr, ("var","index","base"))

close(ds_save)






# ENV["GKSwstype"]="nul"

# # plot and save results
# j=1
# p1 = plot(X0_arr[j] .* perturb[j], CF_arr[j,:]*100, lw=2, label="",
#         color=:black,  yguidefontcolor=:black, guidefontsize=16,
#         ylim=ylimCF, ylabel="Cloud fraction [%]",
#         left_margin=15Plots.mm, right_margin=0Plots.mm,
#         bottom_margin=15Plots.mm, top_margin=5Plots.mm)
# plot!(twinx(), X0_arr[j] .* perturb[j], LWP_arr[j,:]*1e3, lw=2, label="",
#         color=:blue, yguidefontcolor=:blue, ytickfontcolor=:blue,
#         ylim=ylimLWP, yformatter=_->"", guidefontsize=16,
#         xlabel=Xlabel_arr[j])
# j=2
# p2 = plot(X0_arr[j] .* perturb[j], CF_arr[j,:]*100, lw=2, label="",
#         color=:black,  yguidefontcolor=:black, guidefontsize=16,
#         ylim=ylimCF, yformatter=_->"",
#         left_margin=0Plots.mm, right_margin=0Plots.mm,
#         bottom_margin=15Plots.mm, top_margin=5Plots.mm)
# plot!(twinx(), X0_arr[j] .* perturb[j], LWP_arr[j,:]*1e3, lw=2, label="",
#         color=:blue, yguidefontcolor=:blue, ytickfontcolor=:blue,
#         ylim=ylimLWP, yformatter=_->"", guidefontsize=16,
#         xlabel=Xlabel_arr[j])
# j=3
# p3 = plot(1e6 .* X0_arr[j] .* perturb[j], CF_arr[j,:]*100, lw=2, label="",
#         color=:black,  yguidefontcolor=:black, guidefontsize=16,
#         ylim=ylimCF, yformatter=_->"",
#         #left_margin=0Plots.mm, right_margin=30Plots.mm,
#         left_margin=0Plots.mm, right_margin=0Plots.mm,
#         bottom_margin=15Plots.mm, top_margin=5Plots.mm)
# plot!(twinx(), 1e6 .* X0_arr[j] .* perturb[j], LWP_arr[j,:]*1e3, lw=2, label="",
#         color=:blue, yguidefontcolor=:blue, guidefontsize=16, ytickfontcolor=:blue,
#         ylim=ylimLWP, yformatter=_->"", #ylabel="In-cloud LWP [g/m²]", 
#         xlabel=Xlabel_arr[j])
# j=4
# p4 = plot(1e2 .* X0_arr[j] .* perturb[j], CF_arr[j,:]*100, lw=2, label="",
#         color=:black,  yguidefontcolor=:black, guidefontsize=16,
#         ylim=ylimCF, yformatter=_->"",#ylabel="Cloud fraction [%]",
#         #left_margin=15Plots.mm, right_margin=0Plots.mm,
#         left_margin=0Plots.mm, right_margin=0Plots.mm,
#         bottom_margin=15Plots.mm, top_margin=5Plots.mm)
# plot!(twinx(), 1e2 .* X0_arr[j] .* perturb[j], LWP_arr[j,:]*1e3, lw=2, label="",
#         color=:blue, yguidefontcolor=:blue, ytickfontcolor=:blue,
#         ylim=ylimLWP, yformatter=_->"", guidefontsize=16,
#         xlabel=Xlabel_arr[j])
# j=5
# p5 = plot(X0_arr[j] .* perturb[j], CF_arr[j,:]*100, lw=2, label="",
#         color=:black,  yguidefontcolor=:black, guidefontsize=16,
#         ylim=ylimCF, yformatter=_->"",
#         left_margin=0Plots.mm, right_margin=0Plots.mm,
#         bottom_margin=15Plots.mm, top_margin=5Plots.mm)
# plot!(twinx(), X0_arr[j] .* perturb[j], LWP_arr[j,:]*1e3, lw=2, label="",
#         color=:blue, yguidefontcolor=:blue, ytickfontcolor=:blue,
#         ylim=ylimLWP, yformatter=_->"", guidefontsize=16,
#         xlabel=Xlabel_arr[j])
# j=6
# p6 = plot(X0_arr[j] .* perturb[j], CF_arr[j,:]*100, lw=2, label="",
#         color=:black,  yguidefontcolor=:black, guidefontsize=16,
#         ylim=ylimCF, yformatter=_->"",
#         left_margin=0Plots.mm, right_margin=30Plots.mm,
#         bottom_margin=15Plots.mm, top_margin=5Plots.mm)
# plot!(twinx(), X0_arr[j] .* perturb[j], LWP_arr[j,:]*1e3, lw=2, label="",
#         color=:blue, yguidefontcolor=:blue, ytickfontcolor=:blue,
#         ylim=ylimLWP, ylabel="In-cloud LWP [g/m²]", guidefontsize=16,
#         xlabel=Xlabel_arr[j])
# # plot(p1, p2, p3, p4, p5, p6, layout=(2,3), link=:y, size=(1200,600), dpi=200);
# plot(p1, p2, p3, p4, p5, p6, layout=(1,6), link=:y, size=(1500,300), dpi=200);
# savefig(path*"linear_perturbations_around_"*regime*"_state.png");
