#!/bin/bash
set -e

cwd=$(pwd)

# checkout the bundle
rm -rf bundle
git clone $BUNDLE_URL bundle
(cd bundle && git checkout $BRANCH || echo "bundle does not have branch $BRANCH, keeping default branch")
BUNDLE_REPOS=$(grep "ecbuild_bundle(" bundle/CMakeLists.txt | awk '{print $3}')

# prepare the bundle for each repo and checkout the repo source
rm -rf repo.build 
for repo in $MAIN_REPO $LIB_REPOS; do
    echo ""
    
    # setup src directory CMakeLists.txt
    repo_bundle_dir=repo.bundle/$repo
    rm -rf $repo_bundle_dir
    mkdir -p $repo_bundle_dir
    cp $cwd/bundle/CMakeLists.txt $repo_bundle_dir/
    for r in $BUNDLE_REPOS; do
	if [[ $repo != $r ]]; then
	    sed -i "/.* PROJECT $r .*/d" $repo_bundle_dir/CMakeLists.txt
	fi
    done
    	
    # checkout the repos
    repo_src_dir=repo.src/$repo
    rm -rf $repo_src_dir
    url=$(grep "PROJECT $repo" $repo_bundle_dir/CMakeLists.txt | awk '{print $5}' | sed 's|"||g')
    branch=$(grep "PROJECT $repo" $repo_bundle_dir/CMakeLists.txt | awk '{print $8}')
    GIT_LFS_SKIP_SMUDGE=1 git clone -b $branch $url $repo_src_dir
done

# checkout the main repo's branch for each lib, if desired
echo -e "\nChecking out matching repos for \"$MATCH_REPOS\"..."
for repo in $MAIN_REPO $MATCH_REPOS; do
    echo -n "$repo: "
    cd $cwd/repo.src/$repo
    git checkout $BRANCH || echo "$repo: bundle does not have branch $BRANCH, keeping default branch"
done



# determine version information for each repo, to later
# determine if each repo needs to be rebuilt
# TODO add other version info such as docker hash
prev_repo=""
for repo in $LIB_REPOS; do
    cd $cwd/repo.src/$repo
    git_hash=$(git rev-parse HEAD)
    git_branch=$(git symbolic-ref HEAD)    
    rm -rf build.version
    echo "git_branch:$git_branch"     >  build.version
    echo "git_hash:$git_hash"         >> build.version
    echo "build_type:$LIB_BUILD_TYPE" >> build.version
    echo "dependency_repo:$prev_repo" >> build.version    
    prev_repo=$repo
done

# run git-lfs on a subset of the repos
echo -e "\nRunning git-lfs..."
for repo in $LFS_REPOS; do
    echo -n "$repo: "
    cd $cwd/repo.src/$repo
    git lfs fetch
    git lfs checkout
done

