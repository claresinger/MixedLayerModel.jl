import numpy as np
import xarray as xr
# from metpy.interpolate import cross_section
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

ds = xr.open_dataset("cloud_fraction_compare.nc")
dst = xr.open_dataset("transect_output.nc")

fig = plt.figure(figsize=(10,6), constrained_layout=True)
gs0 = gridspec.GridSpec(2, 2, figure=fig, wspace=1)

# predicted CF
gs00 = gridspec.GridSpecFromSubplotSpec(3, 2, width_ratios=[80, 1], height_ratios=[2,6,2], subplot_spec=gs0[0,0])
ax = fig.add_subplot(gs00[:,0], projection=ccrs.PlateCarree())
gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
cmap = "Blues"
clev = np.linspace(0,100,11)
h = ax.contourf(ds.lon, ds.lat, ds.cf*100, clev, cmap=cmap)
ax.set_title("Predicted CF [%]")
ax.set_title("a)", fontweight="bold", loc="left")
ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')

# observed CF
gs00 = gridspec.GridSpecFromSubplotSpec(3, 2, width_ratios=[80, 1], height_ratios=[1.5,7,1.5], subplot_spec=gs0[1,0])
ax = fig.add_subplot(gs00[:,0], projection=ccrs.PlateCarree())
gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
cmap = "Blues"
clev = np.linspace(0,100,11)
h = ax.contourf(ds.lon, ds.lat, ds.obs_cf*100, clev, cmap=cmap)
ax.set_title("Observed CF [%]")
ax.set_title("b)", fontweight="bold", loc="left")
ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')

# minimum in predicted CF
gs00 = gridspec.GridSpecFromSubplotSpec(6,2, width_ratios=[80,1], subplot_spec=gs0[:,1])
ax = fig.add_subplot(gs00[0:2,0], projection=ccrs.PlateCarree())
gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
cmap = "Blues"
clev = np.linspace(0,100,11)
h = ax.contourf(ds.lon, ds.lat, 100*ds.cf_min, clev, cmap=cmap)
ax.set_title("Minimum predicted CF [%]")
ax.set_title("c)", fontweight="bold", loc="left")
ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')

ax = fig.add_subplot(gs00[2:4,0], projection=ccrs.PlateCarree())
gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
cmap = "Blues"
clev = np.linspace(0,100,11)
h = ax.contourf(ds.lon, ds.lat, 100*ds.cf_max, clev, cmap=cmap)
ax.set_title("Maximum predicted CF [%]")
ax.set_title("d)", fontweight="bold", loc="left")
ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')
cbax = fig.add_subplot(gs00[0:2,1])
cb = plt.colorbar(h,cax=cbax,shrink=0.8)

# normalized difference
ax = fig.add_subplot(gs00[4:6,0], projection=ccrs.PlateCarree())
gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
cmap = "RdBu_r"
clev = np.linspace(-2,2,9)
norm_diff = (ds.cf - ds.obs_cf) / np.sqrt(((ds.cf_max - ds.cf_min)/2)**2 + ds.obs_cf_std**2)
h = ax.contourf(ds.lon, ds.lat, norm_diff, clev, cmap=cmap, extend="both")
ax.set_title("Normalized Difference")
ax.set_title("e)", fontweight="bold", loc="left")
ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')
cbax = fig.add_subplot(gs00[4:6,1])
cb = plt.colorbar(h,cax=cbax,shrink=0.8)

plt.savefig("CF_compare_minmax_nodiff_JJA_NEP.png", dpi=300, bbox_inches="tight", facecolor="w")