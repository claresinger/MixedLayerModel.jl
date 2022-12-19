import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib

def map_plot(ax, data, var):    
    gl = ax.gridlines(draw_labels=True)
    gl.right_labels = False
    gl.top_labels = False
    clev = np.linspace(-50,50,21)
    norm = matplotlib.colors.Normalize()
    cmap = "RdBu_r"
    data = data.where(data.CF0 > 0.5)
    ax.set_facecolor("lightgrey")
    h = ax.contourf(
        data.lon,
        data.lat,
        data.sel(var=var).dLWPdX.values * 1e3,
        clev, norm=norm,
        cmap=cmap, extend="both"
        )
    ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')
    ax.set_xlim([-160, -110])
    ax.set_ylim([15, 40])
    return h,ax,gl

path = "experiments/figures/20221211_perturb_from_mean/"
ds = xr.open_dataset(path+"partial_derivatives.nc")

fig, axes = plt.subplots(nrows = 3, ncols = 2, figsize=(8,6), constrained_layout=True,
                sharex=True, sharey=True,
                subplot_kw={'extent': [-160, -110, 15, 40], 'projection':ccrs.PlateCarree()})

keys = ["sst", "WS", "EIS", "D500", "RH500", "CO2"]
letters = ["a","b","c","d","e","f","g","h"]
titles = ["SST", "$U$", "EIS", "$D$", "RH$_{\!+}$", "CO$_2$"]

plt.rcParams.update({"font.size":10})
cax = fig.add_axes([0.3,-0.05,0.4,0.02])
for i, var in enumerate(keys):
    ax = axes.flatten()[i]
    h,ax,gl = map_plot(ax, ds, var)
    ax.set_title(" "+letters[i]+") "+titles[i], loc="left", y=0.85, fontsize=12)
    if i < 4:
        gl.bottom_labels = False
    if i % 2 == 1:
        gl.left_labels = False
    if i == 0:
        cb = plt.colorbar(h,cax=cax,label="$\\Delta$ All-sky LWP [g m$^{-2}$]", orientation="horizontal")
plt.savefig(path+"box_linear_perturb_JJA_NEP.png", dpi=200, bbox_inches="tight", facecolor="w")
ds.close()

# ########################################################

# ## cloud fraction plot
# fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize=(10,3), sharex=True, sharey=True,
#                 subplot_kw={'extent': [-170, -110, 0, 40], 'projection':ccrs.PlateCarree()},
#                 gridspec_kw={'wspace':0.05})
# ds = xr.open_dataset("experiments/data/climatology_for_MLM_std.nc")
# ds = ds.sel(month=[6,7,8]).mean("month")
# ax = axes[0]
# gl = ax.gridlines(draw_labels=True)
# gl.right_labels = False
# gl.top_labels = False
# ax.contourf(ds.lon, ds.lat, ds.low * 100, np.linspace(0,100,11), cmap="Greys")
# ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')
# ax.set_xlim([-160, -110])
# ax.set_ylim([10, 40])
# ds.close()

# path = "experiments/figures/20221211_perturb_from_mean/"
# ds = xr.open_dataset(path+"partial_derivatives.nc")
# ax = axes[1]
# gl = ax.gridlines(draw_labels=True)
# gl.left_labels=False
# gl.right_labels = False
# gl.top_labels = False
# h = ax.contourf(ds.lon,ds.lat,ds.CF0 * 100,np.linspace(0,100,11),cmap="Greys")
# # x,y = np.meshgrid(ds.lon.values, ds.lat.values)
# # ax.plot(x,y,"o",color="yellow")
# ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')
# ax.set_xlim([-160, -110])
# ax.set_ylim([10, 40])

# cax = fig.add_axes([0.35,-0.05,0.3,0.02])
# cb = plt.colorbar(h, cax=cax, orientation="horizontal", label="Cloud fraction [%]")
# plt.savefig(path+"baseline_cloud_fraction.png", dpi=200, bbox_inches="tight", facecolor="w")
# ds.close()
