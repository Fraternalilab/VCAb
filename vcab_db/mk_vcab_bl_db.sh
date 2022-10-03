cd ../seq_db/vcab_db

cp H_seq.fasta h_seq_db
cd h_seq_db
makeblastdb -in H_seq.fasta -dbtype prot
cd ..

cp L_seq.fasta l_seq_db
cd l_seq_db
makeblastdb -in L_seq.fasta -dbtype prot
cd ..

cat H_seq.fasta L_seq.fasta > full_seq_db/all_full_seq.fasta
cd full_seq_db
makeblastdb -in all_full_seq.fasta -dbtype prot
cd ..

cp HV_seq.fasta hv_seq_db
cd hv_seq_db
makeblastdb -in HV_seq.fasta -dbtype prot
cd ..

cp LV_seq.fasta lv_seq_db
cd lv_seq_db
makeblastdb -in LV_seq.fasta -dbtype prot
cd ..

cat HV_seq.fasta LV_seq.fasta > v_seq_db/all_v_seq.fasta
cd v_seq_db
makeblastdb -in all_v_seq.fasta -dbtype prot
cd ..

cd ../../vcab_db

