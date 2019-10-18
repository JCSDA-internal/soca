#!/bin/bash
set -e

cwd=$(pwd)
ln -s ${REPO_CACHE} ./repo.install

# zero out the ccache stats
ccache -z

# for each dependency repo, make install, if needed
rebuild=0
for repo in $LIB_REPOS; do
    echo -e "\n"
    echo "************************************************************"
    echo "  $repo"
    echo "************************************************************"
    
    src_dir=$cwd/repo.src/$repo
    build_dir=$cwd/repo.build/$repo
    install_dir=$cwd/repo.install/$repo

    # set the path for the install dir for subsequent repos to find
    typeset "${repo^^}_PATH=$install_dir"
    export ${repo^^}_PATH

    # do we skip building this repo?
    ver_match=0
    if [[ -e $install_dir/build.version ]]; then
	diff -q $src_dir/build.version $install_dir/build.version && ver_match=1
    fi
    [[ $rebuild == 0 && $ver_match == 0 ]] && rebuild=1
    if [[ $rebuild == 0 ]]; then
	echo "Skipping build of $repo"
	continue
    fi

    # if we are to build
    echo "Building $repo ..."
    rm -rf $build_dir
    mkdir -p $build_dir
    cd $build_dir

    # run ecbuild
    build_opt_var=LIB_BUILD_OPT_${repo^^}
    build_opt=${!build_opt_var}
    ecbuild $src_dir -DCMAKE_INSTALL_PREFIX=${install_dir} -DCMAKE_BUILD_TYPE=${LIB_BUILD_TYPE} \
    	    -DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DENABLE_TESTS=OFF -DBUILD_TESTING=OFF $build_opt

    # build and install
    make -j4
    make install

    # save compilation info
    cp $src_dir/build.version $install_dir/    
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
ecbuild $src_dir -DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DCMAKE_BUILD_TYPE=${MAIN_BUILD_TYPE} -DENABLE_GPROF=ON
make -j4


# how useful was ccache?
ccache -s
