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
    cmap = "RdBu_r"
    # cmax = 200
    # clev = np.linspace(-1*cmax,cmax,6)
    # norm = matplotlib.colors.Normalize()
    # clev = [-2,-0.5,-0.2,-0.05,-0.02,0.02,0.05,0.2,0.5,2]
    # norm = matplotlib.colors.SymLogNorm(linthresh=0.1, vmin=-1, vmax=1)
    clev = [-50, -20, -10, -5, -2, -1, 1, 2, 5, 10, 20, 50]
    norm = matplotlib.colors.SymLogNorm(linthresh=1, vmin=-200, vmax=200)
    h = ax.contourf(
        data.coords["lon"],
        data.coords["lat"],
        data.variables[var],
        clev,norm=norm,
        cmap=cmap,extend="both"
        )
    ax.set_title(data.variables[var].attrs["long_name"])
    ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')
    return h,ax

ds = xr.open_dataset("partial_derivatives.nc")
fig, axes = plt.subplots(nrows = 2, ncols = 3, figsize=(15,6), constrained_layout=True,
                sharex=True, sharey=True, subplot_kw={'projection':ccrs.PlateCarree()})
keys = list(ds.keys())[6:]
print(keys)
letters = ["a","b","c","d","e","f","g","h"]

plt.rcParams.update({"font.size":18})
for i, var in enumerate(keys):
    ax = axes.flatten()[i]
    h,ax = map_plot(ax, ds, var)
    ax.set_title(letters[i]+")", fontweight="bold", loc="left")
cax = fig.add_axes([1.01,0.25,0.015,0.5])
cb = plt.colorbar(h,cax=cax,label="$\\Delta$SWCRE [W m$^{-2}$]")
plt.savefig("box_linear_perturb_JJA_NEP.png", dpi=200, bbox_inches="tight", facecolor="w")