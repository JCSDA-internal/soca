#!/bin/bash
set -e
cat <<EOF

#================================================================================
#================================================================================
# da.post.sh
#  post processing of DA and forecast output.
#  Currently this only does compression and moving of the desired files
#================================================================================

EOF

# required environment variables
envar+=()
envar+=("ANA_TIME")   # date and time of analysis (in any valid "date" command format)
envar+=("DA_ANA")     # path to analysis file
envar+=("DA_BKG")     # path to background file
envar+=("DA_OMB")     # path to the O-F file(s)
envar+=("FCST_DIAG")  # path to diagnostic files produced by the forecast run
envar+=("OUT_DIR")    # experiment directory in which to place output compressed files
envar+=("WORK_DIR")   # temporary working directory for this script

# make sure required env vars exist
set +u
for v in ${envar[@]}; do
    if [[ -z "${!v}" ]]; then
	echo "ERROR: env var $v is not set."; exit 1
    fi
    echo " $v = ${!v}"
done
set -u
echo ""

#================================================================================

ymdh=$(date -ud "$ANA_TIME" +%Y%m%d%H)

# move/compress the model diagnostics
# note, since it is an average, it is at 0Z and not 12Z
echo "compressing model diagnostics..."
mkdir -p $OUT_DIR/diag
ncks -O -7 -L 4 --ppc default=4 --ppc u,v,temp,salt=.2 $FCST_DIAG $OUT_DIR/diag/diag.${ymdh:0:8}00.nc &

# move/compress the analysis
echo "compressing the analysis files..."
mkdir -p $OUT_DIR/ana
ncks -O -7 -L 4 --ppc default=.2 -v Temp,Salt $DA_ANA $OUT_DIR/ana/ana.$ymdh.nc &

# move/compress the background
echo "compressing the background files..."
mkdir -p $OUT_DIR/bkg
ncks -O -7 -L 4 --ppc default=.2 -v Temp,Salt,ave_ssh $DA_BKG $OUT_DIR/bkg/bkg.$ymdh.nc &

wait

# move obs space stats
# TODO, should concat these into fewer files
echo "compressing the obs space files ..."
mkdir -p $OUT_DIR/omb
mkdir -p $WORK_DIR/omb
for f in $DA_OMB/*; do
    f2=${f##*/}
    ncks -O -7 -v .*@MetaData,.*@ombg,.*@oman,.*@Obs,.*@EffectiveQC0 $f $WORK_DIR/omb/$f2 &
done
wait
cd $WORK_DIR/omb
tar -caf $ymdh.tgz *.nc
mv $ymdh.tgz $OUT_DIR/omb
