import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.colors import LinearSegmentedColormap

path = "experiments/figures/20230215_dailytransect_subonly_100days_skip1_1var/"
ds = xr.open_dataset(path+"transect_output_all_LWP.nc")
ds["cf_1var"] = ds.cf_1var.where((ds.cf_1var > 0.05) & (ds.cf_1var < 0.95))
ds = ds.where((ds.zb < 2000) & (ds.zi < 2000) & (ds.zb > 10)) # m

ds["obs_cf_mean"] = ds.obs_cf_mean.mean("time")
ds["obs_cf_std"] = ds.obs_cf_std.mean("time")

###############

fig, axes = plt.subplots(3,2,figsize=(12,6), sharex=True, sharey=False, constrained_layout=True)
plt.rcParams.update({"font.size":15})

vars = ["cf","De","icLWP","dR","zb","zi"]
factor = [100,1,1e3,1,1,1]
colors = ["magenta","firebrick","b","salmon","goldenrod","goldenrod"]
labels = [
    "a) Cloud fraction [%]",
    "b) $\mathscr{D}$, decoupling parameter",
    "c) In-cloud LWP [g m$^{-2}$]",
    "d) Cloud-top radiative cooling [W m$^{-2}$]",
    "e) Cloud base [m]",
    "f) Cloud top [m]",
]

for i,ax in enumerate(axes.flatten()):
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.grid(which="both")
    ax.set_title(labels[i], loc="left")
    # Flatten lon and dR so hexbin sees 1D arrays
    x = np.repeat(ds.lon.values, ds[vars[i]].sizes["time"])
    y = ds[vars[i]].values.flatten()*factor[i]
    hb = ax.hexbin(
        x, y,
        gridsize=15,        # increase for finer resolution
        cmap=LinearSegmentedColormap.from_list("custom",["white", colors[i]]),
        mincnt=1,
        alpha=0.8
    )

    if i == 0:
        ax.plot(ds.lon, ds.obs_cf_mean * 100, lw=3, color="k", label="Observations")
        ax.fill_between(ds.lon, (ds.obs_cf_mean-ds.obs_cf_std)*100, (ds.obs_cf_mean+ds.obs_cf_std)*100, alpha=0.2, color="k")
        ax.plot(ds.lon, ds.cf.mean("time") * 100, lw=3, color="magenta", label="Bulk model")
        # ax.plot(ds.lon, ds.cf_mean.mean("time") * 100, lw=2, color="magenta", ls="--", label="Bulk model, mean forcing")
        ax.set_ylim([0,100])
        ax.legend(loc=4, borderaxespad=0.2, fontsize=14)
    else:
        ax.plot(ds.lon, ds[vars[i]].mean("time")*factor[i], lw=3, color=colors[i])

    if i == 1:
        ax.set_ylim([-1,10])
    if i ==2: 
        ax.set_ylim([0,1000])
    if i == 3:
        ax.set_ylim([0,75])
    if (i == 4) or (i == 5):
        ax.set_ylim([0,1500])

ax.set_xticks([-150, -140, -130, -120])
ax.set_xticklabels(["150°W", "140°W", "130°W", "120°W"])
ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.set_xlim(np.min(ds.lon), np.max(ds.lon))

plt.rcParams.update({"font.size":15})
plt.savefig(path+"daily_transect_means.pdf", dpi=200, bbox_inches="tight")


#######################

fig, axes = plt.subplots(1,1,figsize=(8,4))
plt.rcParams.update({"font.size":15})

ax = axes
ax.tick_params(axis='both', which='major', labelsize=15)
ax.plot(ds.lon, ds.mean("time").cf * 100, lw=3, color="magenta", label="All")
labels = {"sst":"SST", "WS":"$U$", "EIS":"EIS", "D500":"$D_{500}$", "RH500":"RH$_{500}$"}
ls = {"sst":":", "WS":"-", "EIS":"--", "D500":":", "RH500":"-"}
for i,var in enumerate(["sst", "WS", "EIS", "D500", "RH500"]):
    ax.plot(ds.lon, ds.sel(var=var).mean("time").cf_1var * 100, lw=3, ls=ls[var], color="C"+str(i), label=labels[var])
ax.legend(ncol=2, loc=4, borderaxespad=0.2)

ax.set_xlim(np.min(ds.lon), np.max(ds.lon))
ax.set_ylim([0,100])
ax.grid(which="both")

ax.set_xticks([-150, -140, -130, -120])
ax.set_xticklabels(["150°W", "140°W", "130°W", "120°W"])
ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.set_ylabel("Cloud fraction [%]", fontsize=15)

plt.rcParams.update({"font.size":15})
plt.savefig(path+"daily_transect_1var.pdf", dpi=200, bbox_inches="tight")


ds.close()