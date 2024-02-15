#!/usr/bin/env nextflow

/*
 * Defines the pipeline inputs parameters (giving a default value for each for them) 
 * Each of the following parameters can be specified as command line options
 */
import java.time.LocalDate

date = LocalDate.now().format(java.time.format.DateTimeFormatter.ofPattern('yyyyMMdd'))
params.seqres_file = "pdb_seqres_${date}.txt.gz"
params.seqres = "https://ftp.pdbj.org/pub/pdb/derived_data/pdb_seqres.txt.gz"
params.seqres_json = "pdb_protein_seqs_updated_${date}.json"
params.num_result = "num_result"
params.blast_result = "blast_result"
params.result_dir = "result"
params.struct_dir = "pdb_struc"

process download_pdb_seqres {
	/*
	* download the seqres.txt.gz from wwPDB and unzip it
	*/
	publishDir "$baseDir/vcab_db/", mode: 'copy', overwrite: true

	output:
	file "pdb_seqres_${date}.txt"

	input:
	path params.seqres

	// generate the result folders before doing anything
	beforeScript """
        echo ${date}
        if [ ! -d ${baseDir}/vcab_db/${params.num_result} ]; then
                mkdir ${baseDir}/vcab_db/${params.num_result}
        fi
        if [ ! -d ${baseDir}/vcab_db/${params.blast_result} ]; then
                mkdir ${baseDir}/vcab_db/${params.blast_result}
        fi
        if [ ! -d ${baseDir}/vcab_db/${params.result_dir} ]; then
                mkdir ${baseDir}/vcab_db/${params.result_dir}
        fi
        """

	script:
	"""
	wget -O ${params.seqres_file} ${params.seqres}
	gunzip -f ${params.seqres_file}
	"""
}

process anarci_vc {
	/*
	* run ANARCIvc and filter for structures containing both antibody V and C regions
	*/
	debug: 
	true
	
	publishDir "$baseDir/vcab_db/", mode: 'copy', overwrite: true

	output:
	path pdbid_list

	input:
	path seqres_txt

	afterScript """
	cp vnumbering* ${baseDir}/vcab_db/${params.num_result}/
	cp o_cnumbering* ${baseDir}/vcab_db/${params.num_result}/
	"""

	script:
	"""
	#!/usr/bin/env python
	import pandas as pd
	import os
	from anarci_vc import number,run_anarci,chain_type_to_class
	import json
	import sys
	# caution: path[0] is reserved for script path (or '' in REPL)
	sys.path.insert(1, "${workflow.projectDir}/bin")
	from get_paired_seq import collect_all_seq_info, df_column_switch, filter_for_vc

	# Check if the old pdb sequence file existed
	downloaded_pdb_seqs = "${seqres_txt}"
	pdb_seq_fns=[i for i in os.listdir("${baseDir}/vcab_db") if "pdb_seqres_" in i]
	old_pdb_seq_fn=[i for i in pdb_seq_fns if i != downloaded_pdb_seqs]
	old_pdb_seq_fn=old_pdb_seq_fn[0] if old_pdb_seq_fn!=[] else None

	pdb_seq_dict=collect_all_seq_info (downloaded_pdb_seqs,old_pdb_seq_fn)
	with open("${params.seqres_json}", "w") as outfile:
		json.dump(pdb_seq_dict, outfile)

	print ("numbering V&C")
	# The prefix of the output file name (with directory) of the numbering csv file
	vnum_out = f"vnumbering"
	o_cnum_out = f"o_cnumbering"
	cnum_out = f"cnumbering"

	# Transfer the pdb_seq_dict into the list of tuples:[(seqid,seq)], where seqid is "pdbid_chainid"
	pdb_seq_tuple_lst = [(pdbid+"_"+chainid,str(seq)) for pdbid, info in pdb_seq_dict.items() for chainid,seq in info.items()]
	v_result = run_anarci( pdb_seq_tuple_lst, ncpu=1,scheme="imgt", database="ALL", allow=set(["H","K","L"]),assign_germline=True,allowed_species=None,output=True,csv=True,outfile=vnum_out)
	c_result = run_anarci( pdb_seq_tuple_lst, ncpu=1,scheme="imgt_c", database="C_ONLY", allow=set(["H","K","L","H_C1","K_CC","L_CC"]),allowed_species=None,output=True,csv=True,outfile=o_cnum_out)
	
	# Read the numbering csv files
	hvn=pd.read_csv(f"{vnum_out}_H.csv")
	lvn=pd.read_csv(f"{vnum_out}_KL.csv")
	hcn=pd.read_csv(f"{o_cnum_out}_H_C1.csv")
	lcn=pd.read_csv(f"{o_cnum_out}_KL_C.csv")

	# Re-order the columns of c-numbering results
	hcnum=df_column_switch(hcn, [("15B","15A"),("45B","45A")])
	lcnum=df_column_switch(lcn, [("15B","15A"),("45B","45A")])
	hcnum.to_csv(f"{cnum_out}_H_C1.csv",index=False)
	lcnum.to_csv(f"{cnum_out}_KL_C.csv",index=False)

	# Filter for abs with both V and C regions
	hvc,lvc=filter_for_vc (hvn,hcn,lvn,lcn)
	hvc.to_csv("hvc.csv")
	lvc.to_csv("lvc.csv")

	# Get the pdbid of the ab structures containinf V&C H&L
	ab_vc_pdb=generate_pdbid_list (hvc,lvc,out_dir=".")
	pdbid_list = "ab_vc_pdbid_lst.txt"
	"""
}

process download_structs_from_pdb {
	// download structs from pdb
	publishDir "${baseDir}/${params.struct_dir}/full_pdb", mode: 'copy', overwrite: true

	output:
	stdout

	input:
	path pdbid_list

	beforeScript """
	if [ ! -d ${baseDir}/${params.struct_dir}/full_pdb ]; then
		mkdir ${baseDir}/${params.struct_dir}/full_pdb
	fi
	"""

	"""
	mkdir ${params.struct_dir}/
	${params.struct_dir}/pdb_download.sh -f ${pdbid_list} -o ./ -c
	"""
}

process get_paired_ab {
	// identify H and L pairing from the structures
	output:
	stdout

	input:
	path hvc
	path lvc
	path "${params.struct_dir}/full_pdb"

	script:
	"""
   	#!/usr/bin/env python
	import pandas as pd
	import json
	import sys
	# caution: path[0] is reserved for script path (or '' in REPL)
	sys.path.insert(1, "${workflow.projectDir}/bin")
	from get_paired_seq import extract_numbering_from_csv, get_all_pHL_seqs

	hvc=pd.read_csv("${baseDir}/vcab_db/hvc.csv").drop(columns=["Unnamed: 0"])
	lvc=pd.read_csv("${baseDir}/vcab_db/lvc.csv").drop(columns=["Unnamed: 0"])
	pdb_seq_dict=json.load(open("${params.seqres_json}","r"))

	# Reformat the v_result
	n_vresult=extract_numbering_from_csv(f"${params.num_result}/vnumbering_H.csv",f"${params.num_result}/vnumbering_KL.csv")
	ab_vc_pdb = list(set(hvc["pdb"].values)&set(lvc["pdb"].values)) # The intersection of pdbids in both hvc and lvc
	pHL,__=get_all_pHL_seqs(ab_vc_pdb,pdb_seq_dict,n_vresult,hvc,lvc,f"${params.struct_dir}/full_pdb")
	pHL.to_csv("${baseDir}/vcab_db/paired_ab.csv")
	"""
}

process printHead {

	input:
	path outputFile

	script:
	"""
	head -n 10 ${outputFile}
	"""
}


workflow {
	Channel.of(params.seqres).set({seqres_url}) 

	seqres_file = download_pdb_seqres(seqres_url)

	pdbid_list = anarci_vc(seqres_file) 

	download_structs_from_pdb(pdbid_list) | view

	// pair the chains and get a paired_ab.csv file
	Channel.fromPath("${baseDir}/vcab_db/hvc.csv").set({hvc_file})
	Channel.fromPath("${baseDir}/vcab_db/lvc.csv").set({lvc_file})


}
