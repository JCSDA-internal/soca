#!/bin/bash
# concatenate observation files between the letkf observer and solver
# steps. This is temporary, and *should* go away when we go to a halo
# distribution.

cd Data
files=sst.letkf.observer_*.nc
for f in $files; do
    ncks -3 -O -h --mk_rec_dmn nlocs $f -o $f
done
ncrcat -h -O $files sst.letkf.observer.nc
