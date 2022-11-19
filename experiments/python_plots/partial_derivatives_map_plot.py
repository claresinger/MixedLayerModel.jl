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
    clev = [-200,-150,-100,-50,-10,-1,-0.1,0.1,1,10,50,100,150,200]
    norm = matplotlib.colors.SymLogNorm(linthresh=1, vmin=-300, vmax=300)
    # clev = 20
    # norm = matplotlib.colors.Normalize()
    cmap = "RdBu_r"
    h = ax.contourf(
        data.lon,
        data.lat,
        data.sel(var=var).dLWPdX.values * 1e3,
        clev, norm=norm,
        cmap=cmap, extend="both"
        )
    h_cf = ax.contour(data.lon, data.lat, data.CF0.values, levels=[0.2,0.4,0.6], linestyles=["-"], colors="k")
    ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')
    ax.set_xlim([-170, -110])
    ax.set_ylim([0, 40])
    # plt.colorbar(h,ax=ax,label="$\\Delta$SWCRE [W m$^{-2}$]",shrink=0.5)
    return h,ax

path = "experiments/figures/20221118_perturb_from_mean/"
ds = xr.open_dataset(path+"partial_derivatives.nc")
print(ds.dCREdX.min("lon").min("lat"))
print(ds.dCREdX.mean("lon").min("lat"))
print(ds.dCREdX.max("lon").min("lat"))

fig, axes = plt.subplots(nrows = 3, ncols = 2, figsize=(8,9), constrained_layout=True,
                sharex=True, sharey=True,
                subplot_kw={'extent': [-170, -110, 0, 40], 'projection':ccrs.PlateCarree()})
keys = ds["var"].values
print(keys)
letters = ["a","b","c","d","e","f","g","h"]
titles = ["SST", "WS", "D", "RH$_+$", "EIS", "CO$_2$"]

plt.rcParams.update({"font.size":15})
cax = fig.add_axes([1.02,0.35,0.015,0.3])
for i, var in enumerate(keys):
    ax = axes.flatten()[i]
    h,ax = map_plot(ax, ds, var)
    ax.set_title(letters[i]+") "+titles[i], fontweight="bold", loc="left")
    if i == 0:
        # cb = plt.colorbar(h,cax=cax,label="$\\Delta$SWCRE [W m$^{-2}$]")
        cb = plt.colorbar(h,cax=cax,label="allsky-$\\Delta$LWP [g m$^{-2}$]")
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