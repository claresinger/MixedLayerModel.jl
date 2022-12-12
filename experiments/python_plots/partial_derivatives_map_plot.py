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
        # data.sel(var=var).dCREdX.values,
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
print(ds.dCREdX.min("lon").min("lat"))
print(ds.dCREdX.mean("lon").min("lat"))
print(ds.dCREdX.max("lon").min("lat"))

fig, axes = plt.subplots(nrows = 3, ncols = 2, figsize=(8,6), constrained_layout=True,
                sharex=True, sharey=True,
                subplot_kw={'extent': [-160, -110, 15, 40], 'projection':ccrs.PlateCarree()})
keys = ds["var"].values
print(keys)
letters = ["a","b","c","d","e","f","g","h"]
titles = ["SST", "WS", "D", "RH$_+$", "EIS", "CO$_2$"]


plt.rcParams.update({"font.size":15})
cax = fig.add_axes([0.3,-0.05,0.4,0.02])
for i, var in enumerate(keys):
    ax = axes.flatten()[i]
    h,ax,gl = map_plot(ax, ds, var)
    ax.set_title(" "+letters[i]+") "+titles[i], loc="left", y=0.85, fontweight="bold", fontsize=12)
    if i < 4:
        gl.bottom_labels = False
    if i % 2 == 1:
        gl.left_labels = False
    if i == 0:
        # cb = plt.colorbar(h,cax=cax,label="$\\Delta$SWCRE [W m$^{-2}$]", orientation="horizontal")
        cb = plt.colorbar(h,cax=cax,label="$\\Delta$ All-sky LWP [g m$^{-2}$]", orientation="horizontal")
plt.savefig(path+"box_linear_perturb_JJA_NEP.png", dpi=200, bbox_inches="tight", facecolor="w")

## cloud fraction plot
fig, ax = plt.subplots(nrows = 1, ncols = 1,
                subplot_kw={'extent': [-170, -110, 0, 40], 'projection':ccrs.PlateCarree()})
gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
h = ax.contourf(ds.lon,ds.lat,ds.CF0,np.linspace(0,1,11),cmap="Blues")
# x,y = np.meshgrid(ds.lon.values, ds.lat.values)
# ax.plot(x,y,"o",color="yellow")
ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')
ax.set_xlim([-170, -110])
ax.set_ylim([0, 40])
cb = plt.colorbar(h,shrink=0.6,label="Cloud fraction")
plt.savefig(path+"baseline_cloud_fraction.png", dpi=200, bbox_inches="tight", facecolor="w")

# plot std
# file = "experiments/data/box_BCs_daily_JJA_NEP_subonly.nc"
# file = "experiments/data/box_BCs_dailystats_JJA_NEP_subonly.nc"
# ds = xr.open_dataset(file)

# ocean = ds.sst_mean.notnull()
# mask = ocean.where(ocean == True)
# ds = ds * mask

# var = ["sst", "WS", "D500", "RH500", "EIS"]
# fig, axes = plt.subplots(nrows = 3, ncols = 2, figsize=(8,6), constrained_layout=True,
#         sharex=True, sharey=True,
#         subplot_kw={'extent': [-160, -110, 15, 40], 'projection':ccrs.PlateCarree()})

# for i, v in enumerate(var):
#     ax = axes.flatten()[i]
#     gl = ax.gridlines(draw_labels=True)
#     gl.right_labels = False
#     gl.top_labels = False
#     cmap = "Reds"
#     h = ax.contourf(
#         ds.lon,
#         ds.lat,
#         ds[v+"_std"],# / ds[v+"_mean"],
#         #np.linspace(0,1,10),
#         10,
#         cmap=cmap, extend="max"
#     )
#     ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')
#     ax.set_xlim([-160, -110])
#     ax.set_ylim([15, 40])
#     plt.colorbar(h, ax=ax, shrink=0.6)
#     ax.set_title(v)

# plt.delaxes(axes.flatten()[-1])
# plt.show()