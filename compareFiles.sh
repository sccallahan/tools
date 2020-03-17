#!/usr/bin/env bash

#######################################################
## Author: Carson Callahan
## Purpose: Check if any two files are identical
## Date: 2019-09-01
## Notes:
#######################################################

## Arguments
file1=$1
file2=$2

## Code
if cmp -s ${file1} ${file2};
then
	echo "Files are the same"
else
	echo "Files are DIFFERENT"
fi