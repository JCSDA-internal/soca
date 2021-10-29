#!/bin/bash

# delete the existing log file
output_log=testoutput/${COMPARE_TESTNAME}.log
[[ -f "$output_log" ]] && rm $output_log

# run the executable
echo ""
echo "==============================================================================="
echo "Running test executable"
echo "==============================================================================="
${MPI_CMD} $1 $2 $output_log
e=$?
if [[ $e -gt 0 ]]; then
    echo -e "Failed to run executable. Error code: $e \n"
    exit $e
fi

# run compare, if needed
if [[ $SKIP_COMPARE == "FALSE" ]]; then
    echo ""
    echo "==============================================================================="
    echo "Running compare script"
    echo "==============================================================================="

    $COMPARE_SCRIPT $output_log testref/${COMPARE_TESTNAME}.test ${COMPARE_TOL_F} ${COMPARE_TOL_I}
    e=$?
    if [[ $e -gt 0 ]]; then
        echo -e "Failed in the COMPARE step. Error code: $e \n"
        exit $e
    fi
    echo -e "PASSED \n"
fi
