import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

path = "experiments/figures/20221211_dailytransect_subonly_100days_1var/"
ds = xr.open_dataset(path+"transect_output.nc")

plt.figure(figsize=(10,5))
plt.rcParams.update({"font.size":20})

mean_cf = ds.mean("time").cf * 100
std_cf =  ds.std("time").cf * 100
plt.plot(ds.lon, mean_cf, lw=3, color="magenta")
plt.fill_between(ds.lon, mean_cf-std_cf, 
    mean_cf+std_cf, alpha=0.3, color="magenta")

plt.plot(ds.lon, ds.obs_cf_mean * 100, lw=3, color="k")
plt.fill_between(ds.lon, (ds.obs_cf_mean-ds.obs_cf_std)*100, 
    (ds.obs_cf_mean+ds.obs_cf_std)*100, alpha=0.4, color="k")

plt.plot(ds.lon, ds.cf_mean * 100, "--", lw=3, color="magenta", label="All")

labels = {"sst":"SST", "WS":"WS", "D500":"$D_{500}$", "RH500":"RH$_{500}$", "EIS":"EIS"}
for i,var in enumerate(["sst", "WS", "D500", "RH500", "EIS"]):
    plt.plot(ds.lon, ds.sel(var=var).cf_mean_1var * 100, "--", lw=3, color="C"+str(i), label=labels[var])

plt.legend(fontsize=12)
plt.xticks([-145, -140, -135, -130, -125], 
    ["145 °W", "140 °W", "135 °W", "130 °W", "125 °W"])
plt.xlim(np.min(ds.lon), np.max(ds.lon))
plt.ylabel("Cloud fraction [%]")
plt.ylim([0,90])
plt.savefig(path+"daily_transect.png", dpi=200, bbox_inches="tight")


###########################


fig, axes = plt.subplots(2,1,figsize=(10,8), sharex=True, sharey=True)
plt.rcParams.update({"font.size":20})

ax = axes[0]
mean_cf = ds.mean("time").cf * 100
std_cf =  ds.std("time").cf * 100
ax.plot(ds.lon, mean_cf, lw=3, color="magenta")
ax.fill_between(ds.lon, mean_cf-std_cf, 
    mean_cf+std_cf, alpha=0.3, color="magenta")

ax.plot(ds.lon, ds.obs_cf_mean * 100, lw=3, color="k")
ax.fill_between(ds.lon, (ds.obs_cf_mean-ds.obs_cf_std)*100, 
    (ds.obs_cf_mean+ds.obs_cf_std)*100, alpha=0.4, color="k")

ax.plot(ds.lon, ds.cf_mean * 100, "--", lw=3, color="magenta", label="All")

ax.set_title("a)", loc="left", fontsize=15)

ax = axes[1]
ax.plot(ds.lon, ds.obs_cf_mean * 100, lw=3, color="k", label="Obs")
ax.plot(ds.lon, ds.cf_mean * 100, "--", lw=3, color="magenta", label="All")
labels = {"sst":"SST", "WS":"WS", "D500":"$D_{500}$", "RH500":"RH$_{500}$", "EIS":"EIS"}
for i,var in enumerate(["sst", "WS", "D500", "RH500", "EIS"]):
    ax.plot(ds.lon, ds.sel(var=var).cf_mean_1var * 100, "--", lw=3, color="C"+str(i), label=labels[var])

ax.legend(ncol=2, fontsize=15)
ax.set_xticks([-150, -140, -130, -120], 
    ["150 °W", "140 °W", "130 °W", "120 °W"])
ax.set_xlim(np.min(ds.lon), np.max(ds.lon))
ax.set_ylim([0,90])

ax.set_title("b)", loc="left", fontsize=15)

fig.supylabel("Cloud fraction [%]")
plt.savefig(path+"daily_transect_2panel.png", dpi=200, bbox_inches="tight")

ds.close()
