import numpy as np
import xarray as xr
# from metpy.interpolate import cross_section
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

ds = xr.open_dataset("cloud_fraction_compare.nc")
dst = xr.open_dataset("transect_output.nc")

fig = plt.figure(figsize=(12,6), constrained_layout=True)
gs0 = gridspec.GridSpec(2, 2, figure=fig, wspace=1)

# data_squeeze = ds.metpy.parse_cf().squeeze()
# latlon = {"Sc":(30,-120), "Cu":(15,-150), "lat":slice(15,30), "lon":slice(-150,-120)}
# transect = cross_section(data_squeeze, latlon['Sc'], latlon['Cu'], steps=20).set_coords(('lat', 'lon'))

# predicted CF
gs00 = gridspec.GridSpecFromSubplotSpec(3, 2, width_ratios=[80, 1], height_ratios=[2,6,2], subplot_spec=gs0[0])
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
# ax.plot(transect.lon, transect.lat, 'r.-')
cbax = fig.add_subplot(gs00[1,1])
cb = plt.colorbar(h,cax=cbax,shrink=0.6)

# range in predicted CF
gs00 = gridspec.GridSpecFromSubplotSpec(3, 2, width_ratios=[80, 1], height_ratios=[2,6,2], subplot_spec=gs0[1])
ax = fig.add_subplot(gs00[:,0], projection=ccrs.PlateCarree())
gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
cmap = "Blues"
clev = np.linspace(0,45,11)
h = ax.contourf(ds.lon, ds.lat, 100*(ds.cf_max - ds.cf_min)/2, clev, cmap=cmap)
ax.set_title("1$\\sigma$ predicted CF [%]")
# clev = np.linspace(0,100,11)
# h = ax.contourf(ds.lon, ds.lat, 100*ds.cf_max, clev, cmap=cmap)
# h2 = ax.contour(ds.lon, ds.lat, 100*ds.cf_min, [20,50], colors='r')
# ax.clabel(h2, h2.levels, inline=True, fmt=lambda x: "{:.0f} %".format(x), fontsize=10)
# ax.set_title("Minimum (maximum) predicted CF [%]")
ax.set_title("b)", fontweight="bold", loc="left")
ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')
# ax.plot(transect.lon, transect.lat, 'r.-')
cbax = fig.add_subplot(gs00[1,1])
cb = plt.colorbar(h,cax=cbax,shrink=0.6)

# observed CF
gs00 = gridspec.GridSpecFromSubplotSpec(3, 2, width_ratios=[80, 1], height_ratios=[1.5,7,1.5], subplot_spec=gs0[2])
ax = fig.add_subplot(gs00[:,0], projection=ccrs.PlateCarree())
gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
cmap = "Blues"
clev = np.linspace(0,100,11)
h = ax.contourf(ds.lon, ds.lat, ds.obs_cf*100, clev, cmap=cmap)
ax.set_title("Observed CF [%]")
ax.set_title("c)", fontweight="bold", loc="left")
ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')
# ax.plot(transect.lon, transect.lat, 'r.-')
cbax = fig.add_subplot(gs00[1,1])
cb = plt.colorbar(h,cax=cbax,shrink=0.6)

# normalized difference
gs00 = gridspec.GridSpecFromSubplotSpec(3, 2, width_ratios=[80, 1], height_ratios=[2,6,2], subplot_spec=gs0[3])
ax = fig.add_subplot(gs00[:,0], projection=ccrs.PlateCarree())
gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
cmap = "RdBu_r"
clev = np.linspace(-2,2,9)
norm_diff = (ds.cf - ds.obs_cf) / np.sqrt(((ds.cf_max - ds.cf_min)/2)**2 + ds.obs_cf_std**2)
h = ax.contourf(ds.lon, ds.lat, norm_diff, clev, cmap=cmap, extend="both")
ax.set_title("Normalized Difference")
ax.set_title("d)", fontweight="bold", loc="left")
ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')
# ax.plot(transect.lon, transect.lat, 'r.-')
cbax = fig.add_subplot(gs00[1,1])
cb = plt.colorbar(h,cax=cbax,shrink=0.6)

# gs00 = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[80, 1], subplot_spec=gs0[3])
# ax = fig.add_subplot(gs00[0])
# ax.grid(color="grey")

# ax.plot(dst.lon, dst.cf*100, "o-", color="green", label="Predicted CF")
# ax.fill_between(dst.lon, dst.cf_min*100, y2=dst.cf_max*100, color='green', alpha=0.3, label="")
# ax.plot(dst.lon, dst.obs_cf*100, "ko-", label="Observed CF")
# ax.set_ylim([0,90])
# # ax.plot(dst.lon, dst.sst, "o-", color="blue", label="SST")
# # ax.plot(dst.lon, dst.Tsurf, "o-", color="cyan", label="Tair")
# # ax.plot(dst.lon, dst.LHF, "o-", color="red", label="LHF")
# ax.set_xlim([-150,-120])
# xticks = ax.get_xticks()[1:-1]
# xtickstr = [str(int(xi*-1))+"°W" for xi in xticks]
# ax.set_xticks(xticks, xtickstr)
# ax.set_ylabel("Cloud fraction [%]")
# ax.set_title("d)", fontweight="bold", loc="left")
# ax.set_title("NEP transect")
# ax.legend(loc=2)

plt.savefig("CF_compare_JJA_NEP.png", dpi=300, bbox_inches="tight", facecolor="w")