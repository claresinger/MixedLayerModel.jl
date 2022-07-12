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
    cb = plt.colorbar(h,ax=ax,shrink=0.6)
    ax.add_feature(cfeature.LAND, zorder=1, facecolor='black', edgecolor='black')
    return ax

ds = xr.open_dataset("partial_derivatives.nc")
fig, axes = plt.subplots(nrows = 3, ncols = 2, figsize=(8,6), constrained_layout=True,
                sharex=True, sharey=True, subplot_kw={'projection':ccrs.PlateCarree()})
keys = list(ds.keys())[6:]
print(keys)

for i, var in enumerate(keys):
    ax = axes.flatten()[i]
    ax = map_plot(ax, ds, var)
plt.savefig("box_linear_perturb_JJA_NEP.png", dpi=200, bbox_inches="tight", facecolor="w")