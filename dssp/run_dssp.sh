cd dssp_results
for pdb_file in $(ls ../../pdb_struc/chain_pdb/*.pdb); do mkdssp -i $pdb_file -o $(basename $pdb_file .pdb)".dssp"; done
cd ..
