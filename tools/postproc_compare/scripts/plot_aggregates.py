"""Aggregate/derived-field plots: total ice concentration, mean ice thickness,
mean snow thickness — Fortran vs C++ vs (Fortran - C++)."""
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean  # noqa: F401
import xarray as xr

import config as C
from plot_panels import panel_extent, proj_for


def derive(ds):
    cat = C.get_cat_dim(ds)
    aicen = ds["aicen"].sum(dim=cat)
    vicen = ds["vicen"].sum(dim=cat)
    vsnon = ds["vsnon"].sum(dim=cat)
    hice = (vicen / aicen).where(aicen > 0)
    hsno = (vsnon / aicen).where(aicen > 0)
    return {
        "Total Ice Concentration": aicen.where(aicen > 0),
        "Mean Ice Thickness":      hice,
        "Mean Snow Thickness":     hsno,
    }


CMAP_BY_FIELD = {
    "Total Ice Concentration": "cmo.ice",
    "Mean Ice Thickness":      "cmo.deep",
    "Mean Snow Thickness":     "cmo.tempo",
}


def plot_field(name, fortran, cpp, lons, lats):
    diff = fortran - cpp
    lo = float(np.nanmin([fortran.min().item(), cpp.min().item()]))
    hi = float(np.nanmax([fortran.max().item(), cpp.max().item()]))
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        lo, hi = None, None
    dmax = float(np.nanmax(np.abs(diff))) if np.isfinite(np.nanmax(np.abs(diff))) else 0.0
    if dmax == 0.0:
        dmax = 1e-30
    cmap = CMAP_BY_FIELD[name]

    for region in C.REGIONS:
        fig, axes = plt.subplots(
            1, 3, figsize=(20, 6),
            subplot_kw={"projection": proj_for(region)},
        )
        fig.suptitle(f"{name} ({region}) — Fortran vs C++", fontsize=14)
        for ax in axes:
            panel_extent(ax, region)
            ax.add_feature(cfeature.LAND, zorder=0)
            ax.add_feature(cfeature.COASTLINE)
            ax.gridlines(draw_labels=False)

        im0 = axes[0].pcolormesh(lons, lats, fortran.squeeze(),
                                 transform=ccrs.PlateCarree(),
                                 cmap=cmap, vmin=lo, vmax=hi)
        axes[0].set_title("Fortran (soca2cice)")
        fig.colorbar(im0, ax=axes[0], orientation="vertical", shrink=0.7)

        im1 = axes[1].pcolormesh(lons, lats, cpp.squeeze(),
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
        out = C.AGG_DIR / f"{name.replace(' ', '_')}_{region}.png"
        plt.savefig(out, dpi=110)
        plt.close(fig)
        print(f"saved {out}")


def main():
    needed = {"aicen", "vicen", "vsnon"}
    ds_f = C.open_squeezed(C.FORTRAN_OUT, variables=needed)
    ds_c = C.open_squeezed(C.CPP_OUT, variables=needed)
    grid = xr.open_dataset(C.GRIDSPEC)
    lons = grid["lon"].squeeze().values
    lats = grid["lat"].squeeze().values

    f_fields = derive(ds_f)
    c_fields = derive(ds_c)
    for name in f_fields:
        plot_field(name, f_fields[name], c_fields[name], lons, lats)


if __name__ == "__main__":
    main()
