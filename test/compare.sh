#!/usr/bin/env bash
# (C) Copyright 2009-2016 ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


flog=$2.log.out
ftest=$2.test.out

$1 | tee ${flog} && \
grep 'Test     : ' ${flog} > ${ftest} && \
diff -s $2 ${ftest}

