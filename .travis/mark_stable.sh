#!/bin/bash
#================================================================================
# (C) Copyright 2019 UCAR
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# After automated tests pass, this script is used to create an updated
# "stable branch" tag of the CMakeLists.txt file in the bundle.
#================================================================================
set -e

RELEASE_BRANCH=${RELEASE_BRANCH:-release/stable-nightly}

cwd=$(pwd)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# figure out what git hash is associated with each repo in the bundle.
# Note that there are several places where this repo could exist.
bundle_dir=$cwd/bundle
bundle_repos=$(grep "ecbuild_bundle(" $bundle_dir/CMakeLists.txt | awk '{print $3}')
cp $bundle_dir/CMakeLists.txt $bundle_dir/CMakeLists.new
for r in $bundle_repos; do

    echo ""
    echo "Finding hash tag for $r..."
    verfile_cache=$REPO_CACHE/$r/build.version
    verfile_src=$cwd/repo.bundle/$r/build.version
    hash="none"

    # check if repo exists in the build cache
    echo -n "searching build cache... "
    if [[ -e $verfile_cache ]]; then
        hash=$(grep git_hash $verfile_cache || echo "_ none")
        hash=$(echo "$hash" | awk '{print $2}')
    fi
    [[ "$hash" == "none" ]] && echo "NOT found" || echo "FOUND"

    # check the repo source ( for uncached repos, i.e. the main test repo)
    if [[ "$hash" == "none" ]]; then
        echo -n "searching uncached src.. "
        [[ -e $verfile_src ]] \
            && hash=$(grep git_hash $verfile_src || echo "_ none") \
            && hash=$(echo "$hash" | awk '{print $2}')
        [[ "$hash" == "none" ]] && echo "NOT found" || echo "FOUND"
    fi

    # use ecbuild to find the repo (i.e. check if in the container)
    if [[ "$hash" == "none" ]]; then
        echo -n "searching cmake paths... "
        rm -rf $cwd/repo_hash
        mkdir -p $cwd/repo_hash/build
        cp $SCRIPT_DIR/get_repo_hash.cmake  $cwd/repo_hash/CMakeLists.txt
        cd $cwd/repo_hash/build
        hash=$(ecbuild -DREPO_NAME=${r^^} ../  | grep "git_hash" || echo "none")
        [[ "$hash" != "none" ]] && hash=$(echo "$hash" | awk '{print $3}')
        [[ "$hash" == "none" ]] && echo "NOT found" || echo "FOUND"
    fi

    # otherwise, we couldn't find the repo, either the repo in the bundle wasn't used for
    # the ctests, or something bad happened. Oh well.

    # if a git hash was found, update the bundle with a tagged version
    echo "git_hash: $hash"
    if [[ $hash != "none" ]]; then
        hash=${hash:0:7}
        echo "changing $r to $hash"
        sed -i "s/\(.* PROJECT $r .*\) \(BRANCH\|TAG\) .*/\1 TAG $hash \)/" $bundle_dir/CMakeLists.new
    fi
done


# The following git sorcery ensures that the commit which will be pushed onto the release branch
# is the tagged version of what is in develop, and keeps develop and the current release branch
# as its parents (i.e. it looks like a merge, but it is not a real merge)

# get the git commit hash for the relevant involved branches
cd $bundle_dir
ref_develop=$(git rev-parse HEAD)
git checkout $RELEASE_BRANCH
ref_stable=$(git rev-parse HEAD)
ref_common=$(git merge-base $ref_develop $ref_stable)

# do we need to update the ancestors after merging changes?
# (needs to be done when develop has been updated since the last release tag)
amend=""
if [[ "$ref_common" != "$ref_develop" ]]; then
    git merge -s ours develop -m "temporary branch"
    amend="--amend"
fi

# check in the changes
msg="nightly stable branch $(date +%Y-%m-%d)"
mv CMakeLists.new CMakeLists.txt
git add .
git commit -m "$msg" $amend

# push the changes
git push
