"""Aggregate/derived-field plots: total ice concentration, mean ice thickness,
mean snow thickness — background | analysis | legacy | new | (legacy − new)."""
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean  # noqa: F401
import xarray as xr

import config as C
from plot_panels import panel_extent, proj_for


# Field-by-field display config: (analysis-file variable name, colormap,
# vmax for the shared bkg/analysis/legacy/new colorbar, vmax_abs for the diff).
FIELD_CONFIG = {
    "Total Ice Concentration": ("aice_h", "cmo.ice",  1.0, 0.1),
    "Mean Ice Thickness":      ("hi_h",   "cmo.deep", 4.0, 0.2),
    "Mean Snow Thickness":     ("hs_h",   "cmo.tempo", 0.8, 0.1),
}


def derive_from_cats(ds):
    """Aggregate per-category aicen/vicen/vsnon → conc, mean hice, mean hsno."""
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


def derive_from_analysis(ds):
    """Pull already-aggregated analysis fields, masking by aice_h>0."""
    aice = ds["aice_h"]
    hi   = ds["hi_h"]
    hs   = ds["hs_h"]
    return {
        "Total Ice Concentration": aice.where(aice > 0),
        "Mean Ice Thickness":      hi.where(aice > 0),
        "Mean Snow Thickness":     hs.where(aice > 0),
    }


def plot_field(name, bkg, an, fortran, cpp, lons, lats):
    _, cmap, vmax, dmax = FIELD_CONFIG[name]
    diff = fortran - cpp

    for region in C.REGIONS:
        fig, axes = plt.subplots(
            1, 5, figsize=(28, 6),
            subplot_kw={"projection": proj_for(region)},
        )
        fig.suptitle(f"{name} ({region})", fontsize=14)
        for ax in axes:
            panel_extent(ax, region)
            ax.add_feature(cfeature.LAND, zorder=0)
            ax.add_feature(cfeature.COASTLINE)
            ax.gridlines(draw_labels=False)

        panels = [
            ("Background",          bkg,     cmap, 0.0, vmax),
            ("Analysis",            an,      cmap, 0.0, vmax),
            ("Legacy (soca2cice)",  fortran, cmap, 0.0, vmax),
            ("New (PostProcessIce)", cpp,    cmap, 0.0, vmax),
            ("Legacy − New",        diff, "coolwarm", -dmax, dmax),
        ]
        for ax, (title, data, cm, vmin, vmx) in zip(axes, panels):
            im = ax.pcolormesh(lons, lats, data.squeeze(),
                               transform=ccrs.PlateCarree(),
                               cmap=cm, vmin=vmin, vmax=vmx)
            ax.set_title(title)
            fig.colorbar(im, ax=ax, orientation="vertical", shrink=0.7)

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        out = C.AGG_DIR / f"{name.replace(' ', '_')}_{region}.png"
        plt.savefig(out, dpi=110)
        plt.close(fig)
        print(f"saved {out}")


def main():
    needed_cats = {"aicen", "vicen", "vsnon"}
    needed_an   = {"aice_h", "hi_h", "hs_h"}

    ds_bkg = C.open_squeezed(C.BACKGROUND, variables=needed_cats)
    ds_an  = C.open_squeezed(C.ANALYSIS,   variables=needed_an)
    ds_f   = C.open_squeezed(C.FORTRAN_OUT, variables=needed_cats)
    ds_c   = C.open_squeezed(C.CPP_OUT,    variables=needed_cats)

    grid = xr.open_dataset(C.GRIDSPEC)
    lons = grid["lon"].squeeze().values
    lats = grid["lat"].squeeze().values
    grid.close()

    bkg_fields = derive_from_cats(ds_bkg)
    an_fields  = derive_from_analysis(ds_an)
    f_fields   = derive_from_cats(ds_f)
    c_fields   = derive_from_cats(ds_c)

    for name in FIELD_CONFIG:
        plot_field(name,
                   bkg_fields[name], an_fields[name],
                   f_fields[name], c_fields[name],
                   lons, lats)


if __name__ == "__main__":
    main()
