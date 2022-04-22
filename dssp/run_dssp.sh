if [ ! -d dssp_results ];then
	mkdir dssp_results
fi

if [ ! -d aln_results ];then
        mkdir aln_results
fi

cd dssp_results
for pdb_file in $(ls ../../pdb_struc/chain_pdb/*.pdb); do mkdssp -i $pdb_file -o $(basename $pdb_file .pdb)".dssp"; done
cd ..
