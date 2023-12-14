import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import rioxarray

## cloud fraction plot
fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(8,5),
                subplot_kw={'extent': [-170, -110, 0, 40], 'projection':ccrs.PlateCarree()},
                gridspec_kw={'wspace':0.05})

path = "experiments/figures/20230223_CFdaily_200day_skip20x10/"
ds = xr.open_dataset(path+"CF_daily.nc")
ds = ds.where((ds.CF >= 0.09) & (ds.CF <= 0.81))
z_nan = ds.mean("day").rename_dims({"lat":"y","lon":"x"}).rename_vars({"lat":"y","lon":"x"}).CF * 100
z_nan.rio.write_nodata(np.nan, inplace=True)
z_nan.rio.write_crs(4326, inplace=True)
z = z_nan.rio.interpolate_na(method="nearest")

gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
h = ax.contourf(ds.lon,ds.lat,z,np.linspace(0,100,11),cmap="Greys")

path = "experiments/figures/20230830_critical_co2_map/"
da = xr.open_dataset(path+"critical_co2_map.nc")
CS = plt.contour(da.lon, da.lat, da.critCO2, levels=np.arange(200,1600,100))
ax.clabel(CS, CS.levels[::3], inline=True, fmt="%.0f", fontsize=14)

ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')
ax.set_xlim([-145, -110])
ax.set_ylim([15, 40])

cax = fig.add_axes([0.35,-0.02,0.3,0.02])
cb = plt.colorbar(h, cax=cax, orientation="horizontal", label="Cloud fraction [%]")
plt.savefig(path+"critical_co2_map.png", dpi=400, bbox_inches="tight", facecolor="w")
ds.close()

# initial sst plot
fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(8,5),
                subplot_kw={'extent': [-170, -110, 0, 40], 'projection':ccrs.PlateCarree()},
                gridspec_kw={'wspace':0.05})
# file = "experiments/data/crit_co2_map_input.nc"
file = "experiments/data/regional_daily_good_BCs_JJA_NEP_subonly.nc"
da = xr.open_dataset(file)
gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
h = plt.contourf(da.lon, da.lat, da.mean("time").sst, levels=np.arange(285,300,1), cmap="Reds", extend="both")
CS = plt.contour(da.lon, da.lat, da.mean("time").sst, levels = [290], colors="k")
ax.clabel(CS, inline=True, fmt="%.0f", fontsize=14)

ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')
ax.set_xlim([-145, -110])
ax.set_ylim([15, 40])

cax = fig.add_axes([0.35,-0.02,0.3,0.02])
cb = plt.colorbar(h, cax=cax, orientation="horizontal", label="SST [K]")
plt.savefig(path+"initial_sst_map.png", dpi=400, bbox_inches="tight", facecolor="w")
da.close()