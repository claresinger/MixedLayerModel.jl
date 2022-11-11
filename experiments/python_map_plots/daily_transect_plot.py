import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

path = "../figures/20221110_dailytransect/"
ds = xr.open_dataset(path+"transect_output_Cd4.nc")

plt.figure(figsize=(10,5))

mean_cf = ds.mean("time").cf * 100
std_cf =  ds.std("time").cf * 100
plt.plot(ds.lon, mean_cf, lw=3, color="magenta")
plt.fill_between(ds.lon, mean_cf-std_cf, 
    mean_cf+std_cf, alpha=0.5, color="magenta")

plt.plot(ds.lon, ds.obs_cf_mean * 100, lw=3, color="k")
plt.fill_between(ds.lon, (ds.obs_cf_mean-ds.obs_cf_std)*100, 
    (ds.obs_cf_mean+ds.obs_cf_std)*100, alpha=0.5, color="k")

plt.plot(ds.lon, ds.cf_monthly_mean * 100, "--", lw=3, color="magenta")

plt.xticks([-145, -140, -135, -130, -125, -120], 
    ["145 °W", "140 °W", "135 °W", "130 °W", "125 °W", "120 °W"])
plt.xlim(np.min(ds.lon), np.max(ds.lon))
plt.ylabel("Cloud fraction [%]")
plt.ylim([0,90])
plt.savefig(path+"daily_transect.png", dpi=200, bbox_inches="tight")

ds.close()
