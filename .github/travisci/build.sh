#!/bin/bash
#================================================================================
# (C) Copyright 2019-2020 UCAR
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# This script builds each repo, in the order listed in "$LIB_REPO $MAIN_REPO", if
# "prep.sh" has determined that it needs to be rebuilt.
#
# ccache is used to accelerate CXX files between builds. And upstream repos are
# cached in ${REPO_CACHE}/{repo} and not rebuilt if they have not changed
#================================================================================
set -e

cwd=$(pwd)

# zero out the ccache stats
echo -e "\nzeroing out 'ccache' statistics"
ccache -z

# for each dependency repo...
for repo in $LIB_REPOS; do
    bundle_dir=$cwd/repo.bundle/$repo
    src_dir=$cwd/repo.src/$repo
    build_dir=$cwd/repo.build/$repo
    install_dir=${REPO_CACHE}/$repo

    # set the path to the install dir for subsequent repos to find
    repo_underscore=${repo/-/_}
    typeset "${repo_underscore}_DIR=$install_dir"
    export ${repo_underscore}_DIR

    # do we skip building this repo?
    [[ -e $bundle_dir/skip_rebuild ]] && continue

    # if we are to build
    echo -e "\n"
    echo "************************************************************"
    echo "  $repo"
    echo "************************************************************"
    rm -rf $build_dir
    rm -rf $install_dir
    mkdir -p $build_dir
    cd $build_dir

    # run ecbuild
    build_opt_var=BUILD_OPT_${repo_underscore}
    build_opt="$BUILD_OPT ${!build_opt_var}"
    time ecbuild $src_dir -DCMAKE_INSTALL_PREFIX=${install_dir} -DCMAKE_BUILD_TYPE=${LIB_BUILD_TYPE} \
            -DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DENABLE_TESTS=OFF -DBUILD_TESTING=OFF $build_opt

    # build and install
    time make -j4
    time make install

    # save version info for next time
    cp $bundle_dir/build.version $install_dir/
done


# build the main test repo
echo -e "\n"
echo "************************************************************"
echo "  $MAIN_REPO"
echo "************************************************************"
cd $cwd
src_dir=$cwd/repo.src/${MAIN_REPO}
build_dir=$cwd/repo.build/${MAIN_REPO}
mkdir -p  $build_dir
cd $build_dir

build_opt_var=BUILD_OPT_${MAIN_REPO/-/_}
build_opt=${!build_opt_var}

# valgrind and gprof are mutually exclusive
if [[ "$ENABLE_VALGRIND" == "ON" ]]; then
    $build_opt="$build_opt -DSOCA_TESTS_VALGRIND=ON"
else
    $build_opt="$build_opt -DENABLE_GPROF=ON"
fi


time ecbuild $src_dir -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
       -DCMAKE_BUILD_TYPE=${MAIN_BUILD_TYPE} $build_opt
time make -j4


# how useful was ccache?
echo -e "\nccache statistics:"
ccache -s
