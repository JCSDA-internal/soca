"""Quantitative summaries of (Fortran - C++) per variable.

Writes:
  - diff_stats/summary.csv     : per-variable max abs diff, RMS diff, nonzero-diff count
  - diff_stats/<var>_hist.png  : histogram of nonzero diffs
"""
import csv
import numpy as np
import matplotlib.pyplot as plt

import config as C


NONZERO_TOL = 1e-12


def diff_array(ds_f, ds_c, var, cat):
    da_f = ds_f[var]
    da_c = ds_c[var]
    cat_dim_f = C.get_cat_dim(ds_f)
    cat_dim_c = C.get_cat_dim(ds_c)
    if cat_dim_f and cat_dim_f in da_f.dims:
        da_f = da_f.isel({cat_dim_f: cat})
    if cat_dim_c and cat_dim_c in da_c.dims:
        da_c = da_c.isel({cat_dim_c: cat})
    return (da_f.values - da_c.values).ravel()


def stats(diff):
    finite = diff[np.isfinite(diff)]
    if finite.size == 0:
        return dict(n=0, n_nonzero=0, max_abs=0.0, rms=0.0, mean=0.0)
    nz_mask = np.abs(finite) > NONZERO_TOL
    return dict(
        n=int(finite.size),
        n_nonzero=int(nz_mask.sum()),
        max_abs=float(np.max(np.abs(finite))),
        rms=float(np.sqrt(np.mean(finite ** 2))),
        mean=float(np.mean(finite)),
    )


def hist(diff, var, cat_label):
    finite = diff[np.isfinite(diff)]
    nz = finite[np.abs(finite) > NONZERO_TOL]
    if nz.size == 0:
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(nz, bins=80)
    ax.set_yscale("log")
    ax.set_xlabel("Fortran − C++")
    ax.set_ylabel("count (log)")
    ax.set_title(f"{var}{cat_label} diff histogram (n_nonzero={nz.size})")
    out = C.DIFFSTATS_DIR / f"{var}{cat_label}_hist.png"
    plt.tight_layout()
    plt.savefig(out, dpi=110)
    plt.close(fig)
    print(f"saved {out}")


def main():
    ds_f = C.open_squeezed(C.FORTRAN_OUT)
    ds_c = C.open_squeezed(C.CPP_OUT)

    rows = []
    for var in C.PER_CAT_VARS:
        cat_dim = C.get_cat_dim(ds_f)
        cats = range(ds_f.dims[cat_dim]) if (cat_dim and cat_dim in ds_f[var].dims) else [0]
        for cat in cats:
            label = f"_cat{cat + 1}" if (cat_dim and cat_dim in ds_f[var].dims) else ""
            d = diff_array(ds_f, ds_c, var, cat)
            s = stats(d)
            s["var"] = var + label
            rows.append(s)
            hist(d, var, label)
            print(f"{var}{label}: max|d|={s['max_abs']:.3e} rms={s['rms']:.3e} "
                  f"n_nonzero={s['n_nonzero']}/{s['n']}")

    csv_path = C.DIFFSTATS_DIR / "summary.csv"
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["var", "n", "n_nonzero", "max_abs", "rms", "mean"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"saved {csv_path}")


if __name__ == "__main__":
    main()
