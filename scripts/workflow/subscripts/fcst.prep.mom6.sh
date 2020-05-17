#!/bin/bash
set -e
cat << EOF



#================================================================================
#================================================================================
# fcst.prep.sh
#   MOM6/SIS2 ocean/seaice ensemble member forecast preparation.
#   (Stolen from the Hybrid-GODAS repository)
#   Processes control forcing, ensemble perturbations, and climatological
#   corrections, for the atmospheric forcing
#================================================================================
#================================================================================

EOF

# Prerequisites:
#  * Daily forcing files for the appropriate date range must exist in the
#    $FORC_MEAN_FILE (and $FORC_ENS_FILE locations if doing ensemble DA)
#  * the restart files do NOT need to be ready yet, this step can be done
#     before the previous cycle finished
#
# Results:
#  * ensemble surface forcing files are placed in $FORC_DIR
#
# Required environment variables:
envar=()

envar+=("WORK_DIR")         # The temporary directory in which this script
                            #  will do all of its work

envar+=("ENS_SIZE")         # number of ensemble members
if [[ $ENS_SIZE -gt 1 ]]; then
    envar+=("FORC_ENS_FILE")    # path to the input daily flux file ensemble perturbations
fi

envar+=("FCST_START_TIME")  # Datetime for start of the forecast (YYYYMMDDHH)

envar+=("FCST_END_TIME")    # Datetime for end of the forecast (YYYYMMDDHH)

envar+=("FORC_HR")          # The hour at which forcing is specified (usually 12Z)
FORC_HR=${FORC_HR:-12}

envar+=("FORC_MEAN_FILE")   # path to the input daily flux files

envar+=("FORC_VAR")         # comma separated list of variables to process

envar+=("FORC_VAR_POS")     # comma separated list of variables that should be positive

if [[ ! $FORC_PERTURB -eq 0 ]]; then
    envar+=("FORC_VAR_ENS")   # comma separated list of variables that should have ens
                              # perturbations applied
fi

envar+=("FORC_CORR")        # If =1, a monthly bias correction is done to FORC_MEAN_FILE
FORC_CORR=${FORC_CORR:-0}
if [[ ! $FORC_CORR -eq 0 ]]; then

    envar+=("FORC_CORR_DIR")  # path to the monthly flux bias correction files

    envar+=("FORC_CORR_ADD")  # comma separated list of variables that have additive correction done

    envar+=("FORC_CORR_MUL")  # comma separated list of variables that have multiplicative correction
fi

envar+=("FORC_SKIP_LEAP")   # If == 1, no error is thrown if dataset is missing leap day
FORC_SKIP_LEAP=${FORC_SKIP_LEAP:-1}

envar+=("FORC_PERTURB")         # if =1, ensemble perturbations are pulled from
FORC_PERTURB=${FORC_PERTURB:-1} #  FORC_ENS_FILE and added to FORC_MEAN_FILE for each
                                #  ensemble member.
                                # if =0, FORC_MEAN_FILE is not used at all.

envar+=("FORC_RUNOFF")
if [[ ! $FORC_RUNOFF -eq 0 ]]; then

    envar+=("FORC_RUNOFF_CLIM")
    envar+=("FORC_RUNOFF_PERTURB")
    if [[ ! $FORC_RUNOFF_PERTURB -eq 0 ]]; then
	envar+=("FORC_RUNOFF_VAR")
    fi
fi

envar+=("IC_GEN")
if [[ ! $IC_GEN -eq 0 ]]; then

    envar+=("IC_DIR")
fi

#================================================================================
#================================================================================


# make sure required env vars exist
set +u
for v in ${envar[@]}; do
    if [[ -z "${!v}" ]]; then
	echo "ERROR: env var $v is not set."; exit 1
    fi
    echo " $v = ${!v}"
done
set -u
echo -e "\n================================================================================\n"


# deal with comma separated lists
IFS=',' read -ra forc_var <<< "$FORC_VAR"
IFS=',' read -ra forc_var_pos <<< "$FORC_VAR_POS"
if [[ "$ENS_SIZE" -gt 1 ]]; then
    IFS=',' read -ra forc_var_ens <<< "$FORC_VAR_ENS"
fi
if [[ "$FORC_CORR" -eq 1 ]]; then
  IFS=',' read -ra forc_corr_add <<< "$FORC_CORR_ADD"
  IFS=',' read -ra forc_corr_mul <<< "$FORC_CORR_MUL"
  forc_corr_var="${forc_corr_add[@]} ${forc_corr_mul[@]}"
fi

# surface forcing interpolation method ("bil" or "bic")
interp="bil"

# ------------------------------------------------------------
# Initial COnditions
# ------------------------------------------------------------
if [[ $IC_GEN == 1 ]]; then
    echo "linking initial conditions..."
    echo "------------------------------------------------------------"
    IC_WORK_DIR=$WORK_DIR/ic
    if [[ -e "$IC_WORK_DIR" ]]; then
	echo "WARNING: IC_WORK_DIR already exists, removing:"
	echo " $IC_WORK_DIR"
	rm -rf "$IC_WORK_DIR"
    fi
    mkdir -p "$IC_WORK_DIR"

    for mem in $ENS_LIST; do
	mkdir -p $IC_WORK_DIR/mem_$mem;
	# make sure ic file exists
	ic_file=$IC_DIR/mem_$mem.nc
	if [[ ! -e $ic_file ]]; then
	    echo "ERROR: IC file $ic_file does not exist"
	    exit 1
	fi
	ln -s $ic_file $IC_WORK_DIR/mem_$mem/ic.nc
    done
fi


#------------------------------------------------------------
# Forcing
#------------------------------------------------------------

# setup working directory
#------------------------------------------------------------
FORC_WORK_DIR=$WORK_DIR/forc
if [[ -e "$FORC_WORK_DIR" ]]; then
    echo "WARNING: FORC_WORK_DIR already exists, removing:"
    echo " $FORC_WORK_DIR"
    rm -rf "$FORC_WORK_DIR"
fi
mkdir -p "$FORC_WORK_DIR"
cd $FORC_WORK_DIR
mkdir -p work


# determine the actual start/end dates of the forcing files we need.
#  We are likely to need 1 day before, and 1 day after because the daily forcings
#  are usually specified at 12Z. So a forecast starting at Jan02,0Z will need the
#  Jan01 daily file as well.
#------------------------------------------------------------
forc_start_dy=${FCST_START_TIME:0:8}
forc_start_hr=${FCST_START_TIME:8:12}
if [[ $forc_start_hr -lt $FORC_HR ]]; then
    forc_start_dy=$(date "+%Y%m%d" -d "$forc_start_dy - 1 day")
fi

forc_end_dy=${FCST_END_TIME:0:8}
forc_end_hr=${FCST_END_TIME:8:12}
if [[ $forc_end_hr -gt $FORC_HR ]]; then
    forc_end_dy=$(date "+%Y%m%d" -d "$forc_end_dy + 1 day")
fi

echo ""
echo "preparing surface forcing from $forc_start_dy to $forc_end_dy"
echo " forcing variables are: ${forc_var[@]}"
echo ""

mkdir -p mem_0000


#------------------------------------------------------------
# Create the river runoff mean
#------------------------------------------------------------
if [[ $FORC_RUNOFF == 1 ]]; then
    echo "Generating river runoff files..."
    echo "------------------------------------------------------------"
    cp $FORC_RUNOFF_CLIM mem_0000/runoff.nc
fi


# ------------------------------------------------------------
# Create the combined mean forcing file
# TODO: don't need to do this if we are directly using the ensemble
#  files (instead of using them as just perturbations)
#------------------------------------------------------------
mkdir -p work/mean
echo "Generating forcing mean files..."
echo "------------------------------------------------------------"
for f in ${forc_var[@]}; do

    # does this variable need bias correction applied?
    do_bias=0
    if [[ $FORC_CORR -eq 1 ]]; then
	for v in $forc_corr_var; do
	    [[ "$f" == "$v" ]] && do_bias=1
	done
    fi
    if [[ $do_bias == 1 ]]; then
	bias_op=0
	for v in "${forc_corr_add[@]}"; do [[ "$f" == "$v" ]] && bias_op="add"; done
	for v in "${forc_corr_mul[@]}"; do [[ "$f" == "$v" ]] && bias_op="mul"; done
	if [[ "$bias_op" == 0 ]]; then
	    echo "ERROR: $f needs to be specified in either FORC_CORR_MUL or FORC_CORR_ADD"
	    exit 1
	fi
	echo "$(printf %15s $f) (bias corrected, $bias_op)"
    else
	echo "$(printf %15s $f) "
    fi

    # for each date in the range that needs to be processed for this variable
    files=""
    date_cur=$forc_start_dy
    while [[ $(date -d "$date_cur" +%s) -le $(date -d "$forc_end_dy" +%s) ]];do
	file_in=${FORC_MEAN_FILE//#var#/$f}
	file_in=$(date "+$file_in" -d "$date_cur")
	file_out=$file_in

	# make sure the file exists
	# TODO

	#------------------------------------------------------------
	# if we are bias correcting this variable
	#------------------------------------------------------------
	if [[ $do_bias -eq 1 ]]; then
	    file_out=work/mean/corrected.$f.$date_cur.nc

            # from the monthly climatology correction files
            # we need 3 months to create a file file from which the
            # time interpolated correction is then generated
            # get the time, in seconds, of the 3 desired months
            cf1_s=$(date -d "${date_cur:0:6}01 - 1 month" "+%s" )
            cf2_s=$(date -d "${date_cur:0:6}01" "+%s" )
            cf3_s=$(date -d "${date_cur:0:6}01 + 1 month" "+%s")
            cf1=0; # after the next loop, will contain the files paths for the 3 files
            cf2=0;
            cf3=0
            # for each climatology period available
            for f2 in $FORC_CORR_DIR/*.*.monthly.${f}.nc; do
                f3=${f2##*/}   # remove directory
                f3=${f3#*.}   # remove first bit before the dates
                f3=${f3%%.*}  # remove the bit after the dates
                yr1=${f3%%-*} # starting year of this climatology
                yr2=${f3##*-} # ending year of this climatology
                s1=$(date -d "${yr1}0101" +%s) # starting date of this clim, in seconds
                s2=$(date -d "${yr2}1231" +%s) # ending date of this clim, in seconds
                if [[ $cf1 == 0 ]] || [[ $cf1_s -ge $s1 ]]; then cf1=$f2; fi
                if [[ $cf2 == 0 ]] || [[ $cf2_s -ge $s1 ]]; then cf2=$f2; fi
                if [[ $cf3 == 0 ]] || [[ $cf3_s -ge $s1 ]]; then cf3=$f2; fi
            done

            # generate the interpolated bias correction file
            file_corr=work/mean/bias.$f.$cf1_s.$cf2_s.$cf3_s.nc
            if [[ ! -e $file_corr ]]; then
                cdo -s -L mergetime -setyear,1899 -selmon,12 $cf1 -setyear,1900 $cf2 -setyear,1901 -selmon,1 $cf3 $file_corr
            fi

            # apply the correction
            t=$(date "+1900-%m-%d,12:00:00" -d "$date_cur")
            cdo -s -L $bias_op $file_in -inttime,$t $file_corr $file_out
	fi
	# all done bias correcting this date
	files="$files $file_out"
	date_cur=$(date "+%Y%m%d" -d "$date_cur + 1 day")
    done

    ncrcat -O -L 0 $files mem_0000/$f.nc
    ncatted -O -a axis,time,c,c,T -a calendar,,m,c,gregorian mem_0000/$f.nc
done
echo ""


# ------------------------------------------------------------
# If we are doing an ensemble run, process the ensemble perturbations
# ------------------------------------------------------------
ens_list=""
if [[ "$ENS_SIZE" -gt 1 ]]; then
    ens_list=$(seq -s ' ' -f "%04g" 1 $ENS_SIZE)

    echo ""
    echo "Generating ensemble member forcing files..."
    echo "------------------------------------------------------------"

    if [[ $FORC_RUNOFF == 1 ]]; then
	# generate river runoff files
	echo "River Runoff"
	yr=${date_cur:0:4}
	mn=${date_cur:4:2}
	for m in $ens_list; do
	    d=work/ens/mem_$m
	    mkdir -p $d
 	    mkdir -p mem_$m

	    if [[ ! $FORC_RUNOFF_PERTURB -eq 0 ]]; then
		$ROOT_GODAS_DIR/tools/clim_noise.py -nx 1440 -ny 1080 -month $mn -year $yr -seed $m $d/runoff.noise.nc
		cdo -L add $FORC_RUNOFF_CLIM -mul $FORC_RUNOFF_VAR $d/runoff.noise.nc mem_$m/runoff.nc
		# cdo screws some things up, fix them..
		ncatted -O -a axis,i,d,, -a axis,j,d,, -a cartesian_axis,i,o,c,X -a cartesian_axis,j,o,c,Y \
			-a axis,time,d,, -a cartesian_axis,time,o,c,T -a calendar,time,o,c,noleap \
			-a modulo,time,o,c," " mem_$m/runoff.nc
	    else
		cp $FORC_RUNOFF_CLIM mem_$m/runoff.nc
	    fi
	done
    fi

    # Create the combined atmospheric forcing file for each member
    echo "  Generating combined atmospheric forcing file for each member..."
    for m in $ens_list; do
	d=work/ens/mem_$m
	mkdir -p $d

 	for f in ${forc_var_ens[@]}; do
 	    files=()
 	    date_cur=$forc_start_dy
 	    while [[ $(date -d "$date_cur" +%s) -le $(date -d "$forc_end_dy" +%s) ]];do
 		date_next=$(date "+%F" -d "$date_cur + 1 day")
 		file=${FORC_ENS_FILE//#var#/$f}
 		file=${file//#mem2#/${m: -2}}
 		file=$(date "+$file" -d "$date_cur")
 		files+=("$file")
 		date_cur=$date_next
 	    done
 	    ncrcat ${files[@]} $d/$f.nc
 	done
     done

     # Calculate the mean of the ens files,
     echo "  Generating ensemble forcing mean..."
     mkdir -p work/ens_mean
     for f in ${forc_var_ens[@]}; do
	 ncea work/ens/*/$f.nc -o work/ens_mean/$f.nc
     done

     # generate remap weights
     # (assuming all variables use the same grid
     echo '  Generating "ens->mean" interpolation weights...'
     v=${forc_var_ens[0]}
     cdo -s --no_warnings -L gen${interp},mem_0000/$v.nc work/ens_mean/$v.nc work/remap_weights.nc

     # remap the ens means and calculate mean-ens_mean
     echo '  Calculating "mean - ens_mean"...'
     mkdir -p work/ens_offset
     for f in ${forc_var_ens[@]}; do
 	 cdo -s --no_warnings -L remap,mem_0000/$f.nc,work/remap_weights.nc work/ens_mean/$f.nc tmp.nc
	 cdo -s --no_warnings -L sub mem_0000/$f.nc tmp.nc work/ens_offset/$f.nc
     done

     # generate the final individual ensemble forcing files
     echo "  Generating final ensemble members forcings..."
     for m in $ens_list; do
 	d=mem_$m
 	mkdir -p $d

	for f in ${forc_var[@]}; do
	    is_ens=0
 	    for f2 in ${forc_var_ens[@]}; do
		if [[ $f2 == $f ]]; then
		    is_ens=1
		fi
	    done

	    # if this is a variable that should have an ensemble
	    # perturbation added to it...
	    if [[ $is_ens -ne 0 ]]; then
		cdo -s --no_warnings -L add work/ens_offset/$f.nc -remap,mem_0000/$f.nc,work/remap_weights.nc work/ens/mem_$m/$f.nc $d/$f.nc
 		ncatted -O -a axis,time,c,c,T $d/$f.nc
 		ncatted -O -a calendar,,m,c,gregorian $d/$f.nc
	    else
		# otherwise, just copy the mean file
		cp mem_0000/$f.nc $d/$f.nc
	    fi
	done
     done
fi
ens_list="0000 $ens_list"


# ------------------------------------------------------------
# some fields need to be kept positive, check those
# ------------------------------------------------------------
echo "  Checking that the positiveness of the variables..."
echo "     vars: ${forc_var_pos[@]}"
for m in $ens_list; do
    for v in ${forc_var_pos[@]}; do
	f=mem_$m/$v.nc
	cdo -s --no_warnings -L setrtoc,-1e10,0,0 $f $f.2
	mv $f.2 $f
	ncatted -O -a calendar,,m,c,gregorian $f
    done
done
