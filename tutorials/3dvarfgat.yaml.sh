#!/bin/bash
#
# (C) Copyright 2021-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# Creates a 3DVAR FGAT configuration for a 24 hour DA cycle
#

set -e

date_begin=$1 # example: 2018-04-15T00:00:00Z
outeriter=$2
window_begin=$(date '+%Y-%m-%dT%H:%M:%S'Z -d "$date_begin")

cat << EOF > 3dvarfgat.yaml
# common filters used later on
_: &land_mask
  filter: Domain Check
  where:
  - variable: {name: sea_area_fraction@GeoVaLs}
    minvalue: 0.5

cost function:
  cost type: 4D-Var
  window begin: &date_begin $window_begin
  window length: P1D
  analysis variables: &soca_vars  [socn, tocn, ssh, hocn]
  geometry:
    mom6_input_nml: ./inputnml/input.nml
    fields metadata: ./fields_metadata.yml

  model:
    name: PseudoModel
    tstep: PT1H
    state variables: &model_vars [socn, tocn, ssh, hocn, mld, layer_depth]
    states:
EOF

for i in {1..24}
do
    bkg_date=$(date '+%Y-%m-%dT%H:%M:%S'Z -d "$date_begin $(date +%Z) + $i hour")
    pth=PT${i}H
    [[ "$i" -eq "24" ]] && pth=P1D
    cat << EOF >> 3dvarfgat.yaml
    - date: $bkg_date
      basename: ./fcst/
      ocn_filename: ocn.mom6.fc.$window_begin.${pth}.nc
      read_from_file: 1
EOF
done
cat << EOF >> 3dvarfgat.yaml
  background:
    read_from_file: 1
    basename: ./INPUT/
    ocn_filename: MOM.res.nc
    date: &bkg_date $window_begin
    state variables: *model_vars

  background error:
    covariance model: SocaError
    analysis variables: [socn, tocn, ssh]
    date: $window_begin
    bump:
      verbosity: main
      datadir: ./bump
      strategy: specific_univariate
      load_nicas_local: 1
    correlation:
    - name: ocn
      variables: [socn, tocn, ssh]

    variable changes:

    - variable change: VertConvSOCA
      Lz_min: 2.0
      Lz_mld: 1
      Lz_mld_max: 500.0
      scale_layer_thick: 1.5
      input variables: *soca_vars
      output variables: *soca_vars

    - variable change: BkgErrFILT
      ocean_depth_min: 1000 # [m]
      rescale_bkgerr: 1.0
      efold_z: 2500.0       # [m]
      input variables: *soca_vars
      output variables: *soca_vars

    - variable change: BkgErrGODAS
      t_min: 0.25
      t_max: 1.0
      t_dz:  20.0
      t_efold: 500.0
      s_min: 0.0
      s_max: 0.25
      ssh_min: 0.0   # value at EQ
      ssh_max: 0.0   # value in Extratropics
      ssh_phi_ex: 20 # lat of transition from extratropics
      cicen_min: 0.1
      cicen_max: 0.5
      hicen_min: 10.0
      hicen_max: 100.0
      chl_min: 0.001
      chl_max: 30.0
      biop_min: 0.0
      biop_max: 1.0e-6
      input variables: *soca_vars
      output variables: *soca_vars

    - variable change: BalanceSOCA
      dsdtmax: 0.1
      dsdzmin: 3.0e-6
      dtdzmin: 1.0e-6
      nlayers: 1
      input variables: *soca_vars
      output variables: *soca_vars

  observations:
  - obs space:
      name: SeaSurfaceTemp
      obsdatain:  {obsfile: ../obs/sst.nc4}
      obsdataout: {obsfile: ./obs_out/outeriter_${outeriter}/sst.out.nc4}
      simulated variables: [sea_surface_temperature]
    obs operator:
      name: Identity
    obs error:
      covariance model: diagonal
    obs filters:
    - *land_mask

  - obs space:
      name: SeaSurfaceSalt
      obsdatain:  {obsfile: ../obs/sss.nc4}
      obsdataout: {obsfile: ./obs_out/outeriter_${outeriter}/sss.out.nc4}
      simulated variables: [sea_surface_salinity]
    obs operator:
      name: Identity
    obs error:
      covariance model: diagonal
    obs filters:
    - *land_mask

  - obs space:
      name: ADT
      obsdatain:  {obsfile: ../obs/adt.nc4}
      obsdataout: {obsfile: ./obs_out/outeriter_${outeriter}/adt.out.nc4}
      simulated variables: [obs_absolute_dynamic_topography]
    obs operator:
      name: ADT
    obs error:
      covariance model: diagonal
    obs filters:
    - *land_mask

  - obs space:
      name: InsituTemperature
      obsdatain:  {obsfile: ../obs/insitu.T.nc4}
      obsdataout: {obsfile: ./obs_out/outeriter_${outeriter}/insitu.T.out.nc4}
      simulated variables: [sea_water_temperature]
    obs operator:
      name: InsituTemperature
    obs error:
      covariance model: diagonal
    obs filters:
    - *land_mask

  - obs space:
      name: InsituSalinity
      obsdatain:  {obsfile: ../obs/insitu.S.nc4}
      obsdataout: {obsfile: ./obs_out/outeriter_${outeriter}/insitu.S.out.nc4}
      simulated variables: [sea_water_salinity]
    obs operator:
      name: MarineVertInterp
    obs error:
      covariance model: diagonal
    obs filters:
    - *land_mask

variational:
  minimizer:
    algorithm: RPCG
  iterations:
  - geometry:
      mom6_input_nml: ./inputnml/input.nml
      fields metadata: ./fields_metadata.yml
    linear model:
      name: Identity
      increment variables: *soca_vars
      variable change: Identity
      tstep: PT1H
    ninner: 50
    gradient norm reduction: 1e-15
    test: on
    diagnostics:
      departures: ombg
    online diagnostics:
      write increment: true
      increment:
        datadir: 3dvar_out
        date: *bkg_date
        exp: 3dvarfgat.iter$outeriter
        type: incr

output:
  datadir: 3dvar_out
  exp: 3dvarfgat
  type: an

final:
  diagnostics:
    departures: oman
EOF
