#!/bin/bash
#================================================================================
# (C) Copyright 2019-2020 UCAR
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# This script downloads all the repos used by a bundle, and prepares them to be
# compiled individually. (Not how we normally build a bundle, but this allows for
# caching of the upstream repositories).
#
# see ../.travis.yml for a description of the input environment variables that
# need to be defined
#================================================================================
set -e

cwd=$(pwd)

# we assume the MAIN_REPO contains its own bundle
bundle_file=$cwd/repo.src/${MAIN_REPO}/bundle/CMakeLists.txt


# get a list of all repos listed in the bundle
bundle_repos=$(grep "ecbuild_bundle(" $bundle_file | awk '{print $3}')


# prepare a separate bundle for each upstream repo, determine version information,
# and checkout the repo source if needed. (We don't need to checkout, and build
# if the version information matches that in the cached build)
rm -rf repo.build
prev_repo=""
rebuild=0
for repo in $LIB_REPOS $MAIN_REPO; do
    echo ""

    # create a new "bundle" for each repo
    repo_bundle_dir=repo.bundle/$repo
    rm -rf $repo_bundle_dir
    mkdir -p $repo_bundle_dir
    cp $bundle_file $repo_bundle_dir/
    for r in $bundle_repos; do
        if [[ $repo != $r ]]; then
           sed -i "/.* PROJECT $r .*/d" $repo_bundle_dir/CMakeLists.txt
        fi
    done

    # if this is the main test repo, don't download from git, assume
    # that it has already been downloaded by TravisCI
    [[ $repo == $MAIN_REPO ]] && continue

    # determine the branch / git hash tag, based on the branch name in the original bundle
    repo_url=$(grep "PROJECT $repo" $repo_bundle_dir/CMakeLists.txt | awk '{print $5}' | sed 's|"||g')
    repo_branch=$(grep "PROJECT $repo" $repo_bundle_dir/CMakeLists.txt | awk '{print $8}')
    remote_branches=$(git ls-remote -q $repo_url)
    repo_hash=$( echo "$remote_branches" | grep "refs/heads/$repo_branch\$" | awk '{print $1}')

    echo "Repo:    $repo"
    echo "URL:     $repo_url"

    # if this repo is listed in $MATCH_REPOS, see if a repo branch the same name as the main test
    # repo is available, and use that instead
    domatch=0
    for r in $MATCH_REPOS; do [[ $r == $repo ]] && domatch=1 ; done
    if [[ $domatch == 1 ]]; then
        echo -n "* trying to match main test repo branch... "
        line=$(echo "$remote_branches" | grep "refs/heads/$BRANCH\$" || echo "none")
        if [[ "$line" != "none" ]]; then
            echo "match FOUND"
            repo_branch=$BRANCH
            repo_hash=$(echo "$line" | awk '{print $1}')
        else
            echo "match NOT found"
        fi
    fi
    echo "Branch:  $repo_branch"
    echo "githash: $repo_hash"

    # save the version info to a file
    verfile=$repo_bundle_dir/build.version
    rm -rf $verfile
    echo "repo: $repo"                  > $verfile
    echo "git_branch: $repo_branch"    >> $verfile
    echo "git_hash: $repo_hash"        >> $verfile
    echo "build_type: $LIB_BUILD_TYPE" >> $verfile
    echo "dependency_repo: $prev_repo" >> $verfile
    prev_repo=$repo

    # if the repo cache is empty, the version mismatches that in the cache, or
    # an upstream dependency needs to be rebuilt, mark the repo to be rebuilt
    vermatch=0
    verfile_cache=$REPO_CACHE/$repo/build.version
    [[ -e $verfile_cache ]] && diff -q $verfile $verfile_cache && vermatch=1
    [[ $rebuild == 0 && $vermatch == 0 ]] && rebuild=1
    if [[ $rebuild == 0 ]]; then
        touch $repo_bundle_dir/skip_rebuild
        echo "Skipping rebuild"
    fi

    # does the repo need LFS files downloaded?
    skip_lfs=1
    for r in $LFS_REPOS; do [[ $r == $repo ]] && skip_lfs=0 ; done
    [[ $skip_lfs == 0 ]] && echo "Using git-lfs"

    # if LFS files needed, OR rebuild required, download the repo
    if [[ $rebuild == 1 || $skip_lfs == 0 ]]; then
        repo_src_dir=repo.src/$repo
        rm -rf $repo_src_dir
        GIT_LFS_SKIP_SMUDGE=$skip_lfs git clone --recursive -b $repo_branch $repo_url $repo_src_dir
    fi

done
