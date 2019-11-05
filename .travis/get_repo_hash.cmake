# This is a simple cmake script that will use ecbuild to search for a given repository, and
# print out the git hash of that repo if found. This is useful for finding the git version
# of a repo that has been preinstalled in a container

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )
include( ecbuild_system )
ecbuild_use_package( PROJECT ${REPO_NAME} REQUIRED)
message(STATUS "git_hash " ${${REPO_NAME}_GIT_SHA1} )
