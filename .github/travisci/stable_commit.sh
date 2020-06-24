#!/bin/bash
#================================================================================
# (C) Copyright 2019-2020 UCAR
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# The following git sorcery ensures that the commit which will be pushed onto
# the release branch is the tagged version of what is in develop, and keeps
# develop and the current release branch as its parents (i.e. it looks like a
# merge, but it is not a real merge)
#================================================================================

set -e

echo ""

RELEASE_BRANCH=${RELEASE_BRANCH:-release/stable-nightly}

cwd=$(pwd)

# get the hash for this branch
cd repo.src/$MAIN_REPO
ref_develop=$(git rev-parse HEAD)

# get a modified source url... we have to modify it
# with a credentials token so that we can later push to github
url=$(git remote get-url origin)
url=${url/github.com/"${GH_TOKEN}@github.com"}

# check out release branch in a new directory
# (otherwise, if we do it here, this script may change!)
git clone $url $cwd/repo.release
cd $cwd/repo.release
git checkout $RELEASE_BRANCH || git checkout -b $RELEASE_BRANCH

# do we need to update the ancestors after merging changes?
# (needs to be done when develop has been updated since the last release tag)
amend=""
ref_stable=$(git rev-parse HEAD)
ref_common=$(git merge-base $ref_develop $ref_stable)
if [[ "$ref_common" != "$ref_develop" ]]; then
    git merge -s ours origin/$BRANCH -m "temporary branch"
    amend="--amend"
fi

# replace the CMakeLists with the tagged version
cd bundle
cp $cwd/repo.src/${MAIN_REPO}/bundle/CMakeLists.txt .

# check in the changes
git config --global user.email "travis@travis-ci.org"
git config --global user.name  "TravisCI"
msg="nightly stable  $(date +%Y-%m-%d)"
git add .
git commit -m "$msg" $amend
git push --set-upstream origin $RELEASE_BRANCH