import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

path = "experiments/figures/20221211_linear_perturb/"
ds = xr.open_dataset(path+"linear_perturbs.nc")

Xname_arr = ["SST", "V", "EIS", "D", "RH", "CO2"]
Xlabel_arr = ["SST [K]", "$U$ [m s$^{-1}$]", "IS [K]", "$D \\times 10^6$ [s$^{-1}$]", 
    "RH$_{\!+}$ [%]", "CO$_2$ [ppmv]"]

Nvar, Ncase = np.size(ds.var), np.size(ds.base)
plt.rcParams.update({"font.size":8})
fig, axes = plt.subplots(len(Xname_arr) // 2, 2, figsize=(6,6), constrained_layout=True)
letter = ["a","b","c","d","e","f"]

for (i,var) in enumerate(Xname_arr):
    for (j,case) in enumerate(["Sc","Cu"]):
        ax = axes.flatten()[i]
        ax.plot(
            ds.sel(base=case, var=var).Xval.values * (1e6 if var=="D" else 1) * (100 if var=="RH" else 1),
            ds.sel(base=case, var=var).CF.values*100,
            color="k",
            linestyle=":" if case=="Cu" else "-",
            label=case,
        )
        ax.set_ylim([0,100])
        ax.set_title(letter[i]+") "+Xlabel_arr[i], loc="left")
        if i == 0:
            ax.legend(loc=3, labelspacing=0.2, borderaxespad=0.3, handlelength=1.5)

        ax2 = ax.twinx()
        ax2.plot(
            ds.sel(base=case, var=var).Xval.values * (1e6 if var=="D" else 1) * (100 if var=="RH" else 1),
            ds.sel(base=case, var=var).inLWP.values*1e3,
            color="b",
            linestyle=":" if case=="Cu" else "-",
            label=case,
        )
        ax2.set_ylim([0,250])
        ax2.spines['right'].set_color("b")
        ax2.yaxis.label.set_color("b")
        ax2.tick_params(axis='y', colors="b")

        if i % 2 == 0:
            ax2.set_yticks([])
        else:
            ax.set_yticks([])

fig.supylabel("Cloud fraction [%]", fontsize=10)
fig.text(1.02, 0.5, "In-cloud liquid water path [g m$^2$]", 
    va="center", rotation="vertical", color="b", fontsize=10)
plt.savefig(path+"linear_perturb.png",dpi=200,bbox_inches="tight")


# Nvar, Ncase = np.size(ds.var), np.size(ds.base)
# fig, axes = plt.subplots(len(Xname_arr), 2, figsize=(6,10))#, constrained_layout=True)
# letter = [["a","b","c","d","e","f"],["g","h","i","j","k","l"]]

# for (i,var) in enumerate(Xname_arr):
#     for (j,case) in enumerate(["Sc","Cu"]):
#         ax = axes[i,j]
#         ax.plot(
#             ds.sel(base=case, var=var).Xval.values * (1e6 if var=="D" else 1) * (100 if var=="RH" else 1),
#             ds.sel(base=case, var=var).CF.values*100,
#             color="k",
#         )
#         ax.set_ylim([0,100])
#         ax.set_title(letter[j][i]+") ", fontweight="bold", loc="left")

#         if i == 0:
#             fig.text(0.32, 1.0, "Stratocumulus", fontsize=15, color="k", ha="center", fontweight="bold")
#             fig.text(0.72, 1.0, "Cumulus", fontsize=15, color="k", ha="center", fontweight="bold")

#         ax2 = ax.twinx()
#         ax2.plot(
#             ds.sel(base=case, var=var).Xval.values * (1e6 if var=="D" else 1) * (100 if var=="RH" else 1),
#             ds.sel(base=case, var=var).inLWP.values*1e3,
#             color="b",
#         )
#         ax2.set_ylim([0,250])

#         if j == 0:
#             ax2.set_yticks([])
#         else:
#             ax2.spines['right'].set_color("b")
#             ax2.yaxis.label.set_color("b")
#             ax2.tick_params(axis='y', colors="b")
#             ax.set_yticks([])
    
#     fig.text(0.52, 1-0.005-(1/6)*(i+1), Xlabel_arr[i], ha="center", fontsize=15)

# fig.supylabel("Cloud fraction [%]", fontsize=15)
# fig.text(1.02, 0.5, "In-cloud liquid water path [g m2]", va="center", rotation="vertical", 
#     fontsize=15, color="b")

# plt.tight_layout(h_pad=2)
# plt.savefig(path+"linear_perturb.png",dpi=200,bbox_inches="tight")
