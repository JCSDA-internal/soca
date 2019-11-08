#!/bin/bash
#================================================================================
# (C) Copyright 2019 UCAR
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
git clone ${BUNDLE_URL} bundle.stable
cd bundle.stable

# get the git commit hash for the relevant involved branches 
ref_develop=$(git rev-parse HEAD)
git checkout $RELEASE_BRANCH || git checkout -b $RELEASE_BRANCH
ref_stable=$(git rev-parse HEAD)
ref_common=$(git merge-base $ref_develop $ref_stable)

# do we need to update the ancestors after merging changes?
# (needs to be done when develop has been updated since the last release tag)
amend=""
if [[ "$ref_common" != "$ref_develop" ]]; then
    git merge -s ours develop -m "temporary branch"
    amend="--amend"
fi

# setup git user info
git config --global user.email "travis@travis-ci.org"
git config --global user.name  "TravisCI"
url=${BUNDLE_URL/github.com/"${GH_TOKEN}@github.com"}
git remote add origin-auth $url

# check in the changes
msg="nightly stable  $(date +%Y-%m-%d)"
cat $cwd/bundle/CMakeLists.txt > CMakeLists.txt
git add .
git commit -m "$msg" $amend
git push --set-upstream origin-auth $RELEASE_BRANCH
