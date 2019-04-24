#!/usr/bin/env bash
# (C) Copyright 2009-2016 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

# Catch exec errors?
if [ -z "${OOPS_CATCH_EXEC_ERRORS}" ] ; then
  OOPS_CATCH_EXEC_ERRORS=true
fi

# Initialize suffix
SUFFIX=""

# Check if compare.sh has been run with MPI
TESTPE=true
if [ ! -z "${PMIX_RANK}" ] ; then
  # This run uses MPI
  MPI=true

  # Check MPI rank
  if [ ${PMIX_RANK} != 0 ] ; then
    # If MPI rank positive, add a suffix
    SUFFIX="_"${PMIX_RANK}

    # This is not a test PE (only PE 0 will do the comparison)
    TESTPE=false
  fi
else
  # This run does not use MPI
  MPI=false
fi

# Define log and test files
flog=$2${SUFFIX}.log.out

# Run test, check exit status
$1 | tee ${flog}
exit_status=`echo ${PIPESTATUS[0]}`

# Get test lines
if ${TESTPE} ; then
  ftest=$2${SUFFIX}.test.out
  grep 'Test     : ' ${flog} > ${ftest}
fi

# Deal with execution errors
if [ ${exit_status} != 0 ] ; then
  echo "WARNING: exit status is "${exit_status}
  if ${OOPS_CATCH_EXEC_ERRORS} ; then
    exit ${exit_status}
  fi
fi

if ${TESTPE} ; then
  # Get test lines
  ftest=$2${SUFFIX}.test.out
  grep 'Test     : ' ${flog} > ${ftest}

  # Compare files
  NLINES=`diff -s $2 ${ftest} | wc -l`

  # Exit according to n
  if [ ${NLINES} != 1 ] ; then
    diff $2 ${ftest}
    exit 1
  else
    exit 0
  fi
else
  exit 0
fi
