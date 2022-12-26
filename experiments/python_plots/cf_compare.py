import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

## cloud fraction plot
fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize=(10,3), sharex=True, sharey=True,
                subplot_kw={'extent': [-170, -110, 0, 40], 'projection':ccrs.PlateCarree()},
                gridspec_kw={'wspace':0.05})
ds = xr.open_dataset("experiments/data/climatology_for_MLM_std.nc")
ds = ds.sel(month=[6,7,8]).mean("month")
ax = axes[0]
gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
ax.contourf(ds.lon, ds.lat, ds.low * 100, np.linspace(0,100,11), cmap="Greys")
ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')
ax.set_xlim([-160, -110])
ax.set_ylim([10, 40])
ax.set_title("a) Observations, CASCCAD", loc="left")
ds.close()

path = "experiments/figures/20221221_CFdaily_good/"
ds = xr.open_dataset(path+"CF_daily.nc")
print(ds.count())
ds = ds.where((ds.CF >= 0.09) & (ds.CF <= 0.81))
print(ds.count())
# print(ds.min("day").CF)
# print(ds.mean("day").CF)
# print(ds.max("day").CF)
ax = axes[1]
gl = ax.gridlines(draw_labels=True)
gl.left_labels=False
gl.right_labels = False
gl.top_labels = False
h = ax.contourf(ds.lon,ds.lat,ds.mean("day").CF * 100,np.linspace(0,100,11),cmap="Greys")
# x,y = np.meshgrid(ds.lon.values, ds.lat.values)
# ax.plot(x,y,"o",color="yellow")
ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')
ax.set_xlim([-160, -110])
ax.set_ylim([10, 40])
ax.set_title("b) Prediction, bulk model", loc="left")

cax = fig.add_axes([0.35,-0.05,0.3,0.02])
cb = plt.colorbar(h, cax=cax, orientation="horizontal", label="Cloud fraction [%]")
plt.savefig(path+"baseline_cloud_fraction.png", dpi=400, bbox_inches="tight", facecolor="w")
ds.close()
