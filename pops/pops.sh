#!/bin/bash
shopt -s extglob

dir=$@
# The directory containing c_region_only pdb structures

for file in $dir/*.pdb
do
	Rscript --vanilla deltaResidue.R $file $dir ./result
#	rm ./result/!(*deltaSASA_rpopsResidue.txt)
#	echo $file
done
