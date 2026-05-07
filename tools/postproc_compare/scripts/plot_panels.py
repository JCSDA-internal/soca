"""3-panel plots: Fortran output | C++ output | (Fortran - C++) for each variable.

Run:
    python plot_panels.py
"""
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean  # noqa: F401  (registers cmo.* colormaps)
import xarray as xr

import config as C


CMAPS = {
    "aicen":   "cmo.ice",
    "vicen":   "cmo.deep",
    "vsnon":   "cmo.tempo",
    "Tsfcn":   "cmo.amp",
    "sice":    "cmo.haline",
    "qice":    "cmo.thermal",
    "qsno":    "cmo.thermal",
}


def cmap_for(var):
    for key, cm in CMAPS.items():
        if var.startswith(key):
            return cm
    return "viridis"


def select(ds, var, cat):
    da = ds[var]
    cat_dim = C.get_cat_dim(ds)
    if cat_dim and cat_dim in da.dims:
        da = da.isel({cat_dim: cat})
    return da


def mask_by_aicen(da, aicen):
    return da.where(aicen > 0)


def panel_extent(ax, region):
    if region == "Arctic":
        ax.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
    else:
        ax.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())


def proj_for(region):
    return ccrs.NorthPolarStereo() if region == "Arctic" else ccrs.SouthPolarStereo()


def plot_one(var, cat, ds_f, ds_c, lons, lats):
    a_f = select(ds_f, "aicen", cat)
    a_c = select(ds_c, "aicen", cat)
    d_f = mask_by_aicen(select(ds_f, var, cat), a_f)
    d_c = mask_by_aicen(select(ds_c, var, cat), a_c)
    diff = d_f - d_c

    cmap = cmap_for(var)
    lo = float(np.nanmin([d_f.min().item(), d_c.min().item()]))
    hi = float(np.nanmax([d_f.max().item(), d_c.max().item()]))
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        lo, hi = None, None

    dmax = float(np.nanmax(np.abs(diff))) if np.isfinite(np.nanmax(np.abs(diff))) else 0.0
    if dmax == 0.0:
        dmax = 1e-30  # keep colorbar drawable

    cat_label = f"_cat{cat + 1}" if C.get_cat_dim(ds_f) and C.get_cat_dim(ds_f) in ds_f[var].dims else ""

    for region in C.REGIONS:
        fig, axes = plt.subplots(
            1, 3, figsize=(20, 6),
            subplot_kw={"projection": proj_for(region)},
        )
        fig.suptitle(f"{var}{cat_label} ({region}) — Fortran vs C++", fontsize=14)

        for ax in axes:
            panel_extent(ax, region)
            ax.add_feature(cfeature.LAND, zorder=0)
            ax.add_feature(cfeature.COASTLINE)
            ax.gridlines(draw_labels=False)

        im0 = axes[0].pcolormesh(lons, lats, d_f.squeeze(),
                                 transform=ccrs.PlateCarree(),
                                 cmap=cmap, vmin=lo, vmax=hi)
        axes[0].set_title("Fortran (soca2cice)")
        fig.colorbar(im0, ax=axes[0], orientation="vertical", shrink=0.7)

        im1 = axes[1].pcolormesh(lons, lats, d_c.squeeze(),
                                 transform=ccrs.PlateCarree(),
                                 cmap=cmap, vmin=lo, vmax=hi)
        axes[1].set_title("C++ (PostProcessIce)")
        fig.colorbar(im1, ax=axes[1], orientation="vertical", shrink=0.7)

        im2 = axes[2].pcolormesh(lons, lats, diff.squeeze(),
                                 transform=ccrs.PlateCarree(),
                                 cmap="coolwarm", vmin=-dmax, vmax=dmax)
        axes[2].set_title("Fortran − C++")
        fig.colorbar(im2, ax=axes[2], orientation="vertical", shrink=0.7)

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        out = C.RAW_DIR / f"{var}{cat_label}_{region}.png"
        plt.savefig(out, dpi=110)
        plt.close(fig)
        print(f"saved {out}")


def main():
    grid = xr.open_dataset(C.GRIDSPEC)
    lons = grid["lon"].squeeze().values
    lats = grid["lat"].squeeze().values
    grid.close()

    for var in C.PER_CAT_VARS:
        needed = {"aicen", var}
        ds_f = C.open_squeezed(C.FORTRAN_OUT, variables=needed)
        ds_c = C.open_squeezed(C.CPP_OUT, variables=needed)
        cat_dim = C.get_cat_dim(ds_f)
        cat = C.DEFAULT_CAT if (cat_dim and cat_dim in ds_f[var].dims) else 0
        plot_one(var, cat, ds_f, ds_c, lons, lats)
        ds_f.close()
        ds_c.close()


if __name__ == "__main__":
    main()
