import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

path = "experiments/figures/20230215_dailytransect_subonly_100days_skip1_1var/"
ds = xr.open_dataset(path+"transect_output_all.nc")
ds["cf_1var"] = ds.cf_1var.where((ds.cf_1var > 0.05) & (ds.cf_1var < 0.95))

###############

fig, axes = plt.subplots(2,1,figsize=(10,8), sharex=True, sharey=True)
plt.rcParams.update({"font.size":15})

ax = axes[0]
ax.tick_params(axis='both', which='major', labelsize=15)
ax.plot(ds.lon, ds.obs_cf_mean * 100, lw=3, color="k", label="CASCCAD Observations")
ax.fill_between(ds.lon, (ds.obs_cf_mean-ds.obs_cf_std)*100, 
    (ds.obs_cf_mean+ds.obs_cf_std)*100, alpha=0.2, color="k")

mean_cf = ds.mean("time").cf * 100
std_cf =  ds.std("time").cf * 100
ax.plot(ds.lon, mean_cf, lw=3, color="magenta", label="Bulk model, 100 day mean")
# ax.fill_between(ds.lon, mean_cf-std_cf, mean_cf+std_cf, alpha=0.3, color="magenta")
ax.plot(ds.lon, ds.cf_mean * 100, lw=2, color="magenta", ls="--", label="Bulk model, mean forcing")

ax.grid(which="both")
ax.set_title("a) Transect comparison", loc="left")
ax.legend(loc=4, borderaxespad=0.2)

ax = axes[1]
ax.tick_params(axis='both', which='major', labelsize=15)
ax.plot(ds.lon, mean_cf, lw=3, color="magenta", label="All")
labels = {"sst":"SST", "WS":"$U$", "EIS":"EIS", "D500":"$D_{500}$", "RH500":"RH$_{500}$"}
ls = {"sst":":", "WS":"-", "EIS":"--", "D500":":", "RH500":"-"}
for i,var in enumerate(["sst", "WS", "EIS", "D500", "RH500"]):
    ax.plot(ds.lon, ds.sel(var=var).mean("time").cf_1var * 100, lw=3, ls=ls[var], color="C"+str(i), label=labels[var])

ax.legend(ncol=2, loc=4, borderaxespad=0.2)
ax.set_xticks([-150, -140, -130, -120], 
    ["150°W", "140°W", "130°W", "120°W"])
ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.set_xlim(np.min(ds.lon), np.max(ds.lon))
ax.set_ylim([0,100])
ax.grid(which="both")
ax.set_title("b) Single variable forcing transect comparison", loc="left")

fig.supylabel("Cloud fraction [%]", x=0.05)
plt.rcParams.update({"font.size":15})
plt.savefig(path+"daily_transect_timemean.png", dpi=200, bbox_inches="tight")

ds.close()

#########################

# plt.figure(figsize=(10,5))
# plt.rcParams.update({"font.size":15})

# mean_cf = ds.mean("time").cf * 100
# std_cf =  ds.std("time").cf * 100
# plt.plot(ds.lon, mean_cf, lw=3, color="magenta")
# plt.fill_between(ds.lon, mean_cf-std_cf, 
#     mean_cf+std_cf, alpha=0.3, color="magenta")

# plt.plot(ds.lon, ds.obs_cf_mean * 100, lw=3, color="k")
# plt.fill_between(ds.lon, (ds.obs_cf_mean-ds.obs_cf_std)*100, 
#     (ds.obs_cf_mean+ds.obs_cf_std)*100, alpha=0.4, color="k")

# plt.plot(ds.lon, ds.cf_mean * 100, "--", lw=3, color="magenta", label="All")

# labels = {"sst":"SST", "WS":"$U$", "D500":"$D_{500}$", "RH500":"RH$_{500}$", "EIS":"EIS"}
# for i,var in enumerate(["sst", "WS", "EIS", "D500", "RH500"]):
#     plt.plot(ds.lon, ds.sel(var=var).cf_mean_1var * 100, "--", lw=3, color="C"+str(i), label=labels[var])

# plt.legend(fontsize=12)
# plt.xticks([-145, -140, -135, -130, -125], 
#     ["145 °W", "140 °W", "135 °W", "130 °W", "125 °W"])
# plt.xlim(np.min(ds.lon), np.max(ds.lon))
# plt.ylabel("Cloud fraction [%]")
# plt.ylim([0,90])
# plt.savefig(path+"daily_transect.png", dpi=200, bbox_inches="tight")


###########################


# fig, axes = plt.subplots(2,1,figsize=(10,8), sharex=True, sharey=True)
# plt.rcParams.update({"font.size":15})

# ax = axes[0]
# ax.plot(ds.lon, ds.obs_cf_mean * 100, lw=3, color="k", label="Observations, CASCCAD")
# ax.fill_between(ds.lon, (ds.obs_cf_mean-ds.obs_cf_std)*100, 
#     (ds.obs_cf_mean+ds.obs_cf_std)*100, alpha=0.2, color="k")

# mean_cf = ds.mean("time").cf * 100
# std_cf =  ds.std("time").cf * 100
# ax.plot(ds.lon, mean_cf, lw=3, color="magenta", label="Prediction, bulk model")
# ax.fill_between(ds.lon, mean_cf-std_cf, 
#     mean_cf+std_cf, alpha=0.3, color="magenta")

# # ax.plot(ds.lon, ds.cf_mean * 100, "--", lw=3, color="magenta", label="All")
# ax.grid(which="both")
# ax.set_title("a) Transect comparison", loc="left")
# ax.legend(loc=4)

# ax = axes[1]
# ax.plot(ds.lon, ds.obs_cf_mean * 100, lw=3, color="k")
# ax.fill_between(ds.lon, (ds.obs_cf_mean-ds.obs_cf_std)*100, 
#     (ds.obs_cf_mean+ds.obs_cf_std)*100, alpha=0.2, color="k")

# ax.plot(ds.lon, ds.cf_mean * 100, "--", lw=3, color="magenta", label="All")
# labels = {"sst":"SST", "WS":"$U$", "EIS":"EIS", "D500":"$D_{500}$", "RH500":"RH$_{500}$"}
# ls = {"sst":":", "WS":"-", "EIS":"--", "D500":":", "RH500":"-"}
# for i,var in enumerate(["sst", "WS", "EIS", "D500", "RH500"]):
#     ax.plot(ds.lon, ds.sel(var=var).cf_mean_1var * 100, lw=3, ls=ls[var], color="C"+str(i), label=labels[var])

# ax.legend(ncol=2)
# ax.set_xticks([-150, -140, -130, -120], 
#     ["150°W", "140°W", "130°W", "120°W"])
# ax.xaxis.set_minor_locator(MultipleLocator(5))
# ax.set_xlim(np.min(ds.lon), np.max(ds.lon))
# ax.set_ylim([0,100])
# ax.grid(which="both")
# ax.set_title("b) Single variable forcing transect comparison", loc="left")

# fig.supylabel("Cloud fraction [%]", x=0.05)
# plt.savefig(path+"daily_transect_2panel.png", dpi=200, bbox_inches="tight")

# ds.close()
