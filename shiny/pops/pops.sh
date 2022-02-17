#!/bin/bash
shopt -s extglob

dir=$@
# The original structures are in two file folders: kappa and lambda

for file in $dir/*.pdb
do
	Rscript --vanilla deltaResidue.R $file $dir ./result
	rm ./result/!(*deltaSASA_rpopsResidue.txt)
#	echo $file
done
