if [ ! -d ../dssp/dssp_results ];then
	mkdir ../dssp/dssp_results
fi

if [ ! -d ../dssp/aln_results ];then
        mkdir ../dssp/aln_results
fi

cd ../dssp/dssp_results
for pdb_file in $(ls ../../pdb_struc/chain_pdb/*.pdb); do mkdssp -i $pdb_file -o $(basename $pdb_file .pdb)".dssp"; done
cd ..
