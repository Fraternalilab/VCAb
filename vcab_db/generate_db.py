# To generate the VCAb database
# Dongjun Guo, Aug.2022


import pandas as pd
import numpy as np
import json
import requests
import re
from Bio.Blast.Applications import NcbiblastpCommandline as ncbiblp

import Bio.PDB
from Bio.PDB import MMCIFParser, PDBParser, NeighborSearch, Selection, StructureBuilder
from Bio.PDB.PDBIO import PDBIO
from Bio.SeqUtils import seq1
from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

import os
import argparse
import time
from collections import Counter

import warnings
warnings.filterwarnings("ignore")

# 0. Get the paired antibody chains
os.system("python get_paired_seq.py")

hvnum=pd.read_csv("./num_result/vnumbering_H.csv")
lvnum=pd.read_csv("./num_result/vnumbering_KL.csv")

############ 1. Collapse & Get antigen chain & Get other information via PDBe API ###########

### 1.1 Get Information from PDBe API & Collapse entries based on coordinate sequence of H&L.

def collapse_by_coordinate_seqs (o_df):
    result= [] #collapsed df: nested list
    df=o_df.copy()
    df=df.sort_values(by=["pdb","Hchain","Lchain"])
    df=df.reset_index(drop=True)

    for i in df.index:
        new=df.loc[i]
        n_pdb,n_Hchain,n_Lchain=df.loc[i]['pdb':'Lchain']
        n_htseq=df.loc[i,"H_coordinate_seq"]
        n_ltseq=df.loc[i,"L_coordinate_seq"]
        if i ==0:
            result.append(new)
        else:
            o_pdb,o_Hchain,o_Lchain=result[-1][0:3]
            #return result[-1]
            o_htseq=result[-1]["H_coordinate_seq"]
            o_ltseq=result[-1]["L_coordinate_seq"]
            if n_pdb==o_pdb and n_htseq==o_htseq and n_ltseq==o_ltseq:
                # If coordinate sequences are the same, only attach the chain names to the old(the last) line
                result[-1]["Hchain"] += f";{n_Hchain}" # add Hchain name
                result[-1]["Lchain"] += f";{n_Lchain}" # add Lchain name
            else:
                # If coordinate sequences are different, attach a seperate new line
                result.append(new)
    result_df=pd.DataFrame(result)
    # Add the "iden_code" column
    result_df["iden_code"]=[result_df.loc[i,'pdb']+'_'+result_df.loc[i,'Hchain'][0]+result_df.loc[i,'Lchain'][0] for i in result_df.index]

    return result_df

def get_info_from_pdbe(pdb):
    """
    return the information from PDBe API:
    1. title
    2. release date
    3. experimental method
    4. resolution
    5. carbohydrate name&chainID
    6. mol_info: molecule type, molecule weight,etc.
    7. PDB species for chains
    """

    result={}

    URL= 'https://www.ebi.ac.uk/pdbe/api/'
    n_pdb=pdb.lower()

    sum_data=requests.get(URL + 'pdb/entry/summary/' + n_pdb)
    if sum_data.status_code==200:
        sum_info=sum_data.json()[n_pdb]

        result["title"]=sum_info[0]["title"]
        result["release_date"]=sum_info[0]["release_date"]
        result["method"]=",".join(sum_info[0]["experimental_method"])
    else:
        result["title"]="error"
        result["release_date"]="error"
        result["method"]="error"

    exp_data=requests.get(URL + 'pdb/entry/experiment/' + n_pdb)
    if exp_data.status_code==200:
        exp_info=exp_data.json()[n_pdb]
        if "resolution" in exp_info[0].keys():
            result["resolution"]=exp_info[0]["resolution"]
        else:
            result["resolution"]="NA"
    else:
        result["resolution"]="error"

    # Get carbohydrate information:
    ligand_data=requests.get(URL+'pdb/entry/carbohydrate_polymer/'+n_pdb)
    if ligand_data.status_code==200:
        sugar_result={}
        ligand_info=ligand_data.json()[n_pdb]
        for l_dict in ligand_info:
            l_name=l_dict['molecule_name']
            l_chain_id=l_dict['in_chains']
            if l_name not in sugar_result.keys():
                sugar_result[l_name]=[]
            sugar_result[l_name]+=l_chain_id
        sugar_result2={k:"".join(v) for k,v in sugar_result.items()}
        sugar_result_str="; ".join([f"{k} (chain_id:{v})" for k,v in sugar_result2.items()])
        result["carbohydrate"]= sugar_result_str


    else:
        result["carbohydrate"]= ""

    # Get mol_info:
    mol_data=requests.get(URL+'pdb/entry/molecules/'+n_pdb)
    mol_type_result={}
    if mol_data.status_code==200:
        mol_info=mol_data.json()[n_pdb]
        for m_dict in mol_info:
            chains=m_dict['in_chains'] # author chainID (the chainID used in biopython)
            #chains_asyms=m_dict["in_struct_asyms"]
            mol_name=m_dict['molecule_name'][0]
            mol_type=m_dict['molecule_type']

            mol_chem_id=""
            mol_weight=np.nan
            if "chem_comp_ids" in m_dict.keys():
                mol_chem_id=m_dict["chem_comp_ids"]

            if "weight" in m_dict.keys():
                mol_weight=m_dict["weight"]

            species=""
            if 'source' in m_dict.keys():
                species_info=m_dict['source']
                if len(species_info)!=0 and 'organism_scientific_name' in species_info[0].keys():
                    species=species_info[0]['organism_scientific_name']

            for c in chains:
                if c not in mol_type_result.keys():
                    mol_type_result[c]=(mol_name,mol_type,mol_chem_id,mol_weight,species)

    result["mol_info"]=mol_type_result

    return result


### 1.2 Get antigen chain
def map_pure_seq_pos_to_aln_pos (pure_pos,aln_seq):
    # return the aln_pos of the corresponding pure_seq
    # aln_seq: the seq from the aln object, which containing the gaps
    # both the numbering of pure_pos and aln_seq starts from 0
    pure_seq_counter=0
    for aln_pos,r in enumerate(aln_seq):
        if r!="-":
            if pure_seq_counter==pure_pos:
                return aln_pos
            pure_seq_counter+=1

def map_aln_pos_to_pure_seq_pos (aln_pos,aln_seq):
    # return the pure_seq_pos, given the aln_pos and aln_seq
    # the numbering of both aln_pos & pure_seq_pos starts from 0
    str_before_aln_pos=aln_seq[0:aln_pos+1]
    pure_str=str_before_aln_pos.replace("-","")
    ## get the pure seq without gap "-" before aln_pos
    ## thus, the residues at the pure_pos given by this function would either be
    ## the residue at aln_pos or be the residue at the the left side of aln_pos
    ## (in case where character at aln_pos is a gap)

    pure_seq_pos=len(pure_str)-1
    return pure_seq_pos

def map_imgt_numbering_to_residue_info(id_code,chain_obj,num_df,pdb_dir):
    """
    Returned a list in this format {imgt_numbering:res_obj}
    :id_code: pdbid_oneChainId, used in num_df
    :args chain: Bio.PDB.chain object
    :args num_df: the df outputed by anarci_c containing the C_numbering results
    :args pdb_dir: the directory of the c_pdb files
    """

    num_info=num_df.loc[num_df["Id"]==id_code]
    pure_num_info0=[(numbering,val.values[0]) for (numbering,val) in num_info.iloc[:,13:].items() if (val.values[0] not in["deleted","-"])]
    # pure_num_info0 removes positions generating gaps(both "deleted" and. "-") in the coor_seq
    pure_num_info={numbering:[pure_pos,res] for pure_pos,(numbering,res) in enumerate(pure_num_info0)}

    #return pure_num_info
    existed_numbering=list(pure_num_info.keys())
    num_seq_frag="".join([i[1] for i in list(pure_num_info.values())])

    coor_seq_info=[res for res in chain_obj if res.resname !="X" and res.id[0]==' ']
    #return coor_seq_info
    coor_seq=seq1(''.join(residue.resname for residue in coor_seq_info))

    # 1. Perform the pairwise alignment
    matrix = matlist.blosum62
    pair_aln=pairwise2.align.globalds(num_seq_frag, coor_seq, matrix,-10,-0.5,penalize_end_gaps=False)[0]
    # use blosum62 as the matrix, gap penalties follow the default value of EMBOSS needle
    # just take the first one as the aln result
    #return pair_aln
    alned_num_seq=pair_aln.seqA
    alned_coor_seq=pair_aln.seqB


    # 2. Convert IMGT_num_pos into pure_pos of coor_seq
    # IMGT_num --> pure_pos (imgt_seq) -->aln_pos (imgt_seq)=aln_pos(coor_seq) -->pure_pos (coor_seq)
    result={}

    for i in existed_numbering:
        try:
            num=int(i)
            insertion=""
        except:
            num=int(i[0:-1])
            insertion=i[-1]

        # pure_pos:
        num_seq_pure_pos=pure_num_info[i][0]


        # aln_positions:
        num_seq_aln_pos=map_pure_seq_pos_to_aln_pos(num_seq_pure_pos,alned_num_seq)

        coor_seq_aln_pos=num_seq_aln_pos
        #test.append(coor_seq_aln_pos)
        #if coor_seq_aln_pos!=None: #and alned_coor_seq[coor_seq_aln_pos]=="-":
            #print (alned_num_seq[num_seq_aln_pos])
            #print (i,pure_num_info[i][1],num_seq_aln_pos,coor_seq_aln_pos,alned_coor_seq[coor_seq_aln_pos])

        if coor_seq_aln_pos==None or alned_coor_seq[coor_seq_aln_pos]=="-":
            # if the residue in the alned position of coor_seq is gap
            # this means this residue is missing in the coordinate sequence
            result[i]=[np.nan,num,insertion]

        else:
            # pure_position:
            coor_seq_pure_pos=map_aln_pos_to_pure_seq_pos(coor_seq_aln_pos,alned_coor_seq)

            res_obj=coor_seq_info[coor_seq_pure_pos]

            result[i]=[res_obj,num,insertion]

    return result

def find_nearby_chains (pdbid,h,l,mol_info,hvnum,lvnum,pdb_dir):
    """
    Get the chain ID of the chains close to the CDR loop fragments (within 7.5 A to the alpha-carbon of CDR)
    : args pdbid: 4-letter pdb_id
    : args h,l: 1-letter chain id for H,L respectively
    : args mol_info: the molecule info from PDBeAPI, derived from the function get_info_from_pdbe
    : args hvnum,lvnum: the numbering result of VH&VL from anarci_vc, stored in csv file, now should be dataframe
    : args pdb_dir: the directory storing the full_pdb (mmcif) file
    """

    # Get the structural_object of the chain
    parser=MMCIFParser()
    structure=parser.get_structure(pdbid,f"{pdb_dir}/{pdbid}.cif")
    hchain_obj=structure[0][h]
    lchain_obj=structure[0][l]

    hid=f"{pdbid}_{h}"
    lid=f"{pdbid}_{l}"

    h_residues=map_imgt_numbering_to_residue_info(hid,hchain_obj,hvnum,pdb_dir)
    l_residues=map_imgt_numbering_to_residue_info(lid,lchain_obj,lvnum,pdb_dir)

    # Extract atoms of CDR loops:
    # IMGT CDR1 position:from 27 to 38, including both ends
    # IMGT CDR2 position:from 56 to 65, including both ends
    # IMGT CDR3 position:from 105 to 117, including both ends
    h_cdr1=[i[0] for i in h_residues.values() if i[1]>=27 and i[1]<=38]
    h_cdr2=[i[0] for i in h_residues.values() if i[1]>=56 and i[1]<=65]
    h_cdr3=[i[0] for i in h_residues.values() if i[1]>=105 and i[1]<=117]

    l_cdr1=[i[0] for i in l_residues.values() if i[1]>=27 and i[1]<=38]
    l_cdr2=[i[0] for i in l_residues.values() if i[1]>=56 and i[1]<=65]
    l_cdr3=[i[0] for i in l_residues.values() if i[1]>=105 and i[1]<=117]

    cdrs=h_cdr1+h_cdr2+h_cdr3+l_cdr1+l_cdr2+l_cdr3

    query_atoms=Selection.unfold_entities(cdrs,'A')
    query_atoms=[atom for atom in query_atoms if atom.name=="CA"]

    # Extract atoms of non-antibody chains
    all_h=[i.split("_")[1] for i in hvnum.loc[map(lambda x:x[0:4]==pdbid,hvnum["Id"])]["Id"].values]
    all_l=[i.split("_")[1] for i in lvnum.loc[map(lambda x:x[0:4]==pdbid,lvnum["Id"])]["Id"].values]
    all_ab_chains=all_h+all_l

    other_chains=[c for c in structure[0].get_chains() if c.id not in all_ab_chains]
    #return [c.id for c in structure[0].get_chains()]
    #return [c.id for c in other_chains]

    # Remove water/glycerol chains & chains with mol_weight less than 50
    other_chains=[c for c in other_chains if ("HOH" not in mol_info[c.id][2]) and ("GOL" not in mol_info[c.id][2]) and (mol_info[c.id][3]>50)]

    searching_res_lst=[r for r in Selection.unfold_entities(other_chains,'R') if r.id[0]!="W"]
    #return searching_res_lst
    searching_atom_lst=Selection.unfold_entities(searching_res_lst,'A')
    #return searching_atom_lst

    if len(searching_atom_lst)!=0:
        # If there are other chains available except the antibody chains
        ns=NeighborSearch(searching_atom_lst)

        nearby_chains={chain for query_atom in query_atoms
                       for chain in ns.search(query_atom.coord, 7.5, 'C')}
        remove_chains=[]

        # If multiple polypeptide antigen chains is catched, keep the one with the longest fragment within 7.5 A.
        peptides=[chain for chain in nearby_chains if mol_info[chain.id][1]=="polypeptide(L)"]
        if len(peptides)>1:
            c_info={}
            for c in peptides:
                nearby_residues={residue for query_atom in query_atoms
                            for residue in ns.search(query_atom.coord, 7.5, 'R') if residue in c.get_residues()}
                c_info[c]=len(nearby_residues)
            remove_chains=[k.id for k,v in c_info.items() if v!=max(c_info.values())]

        return "".join([chain.id for chain in nearby_chains if chain.id not in remove_chains])

    else:
        return ""

def get_antigen_info(pdbid,Hchains,Lchains,mol_info,pdb_dir):
    result=[]
    Hchains=Hchains.split(";")
    Lchains=Lchains.split(";")

    for h,l in zip(Hchains,Lchains):
        # For each H/L chain pair:
        sub_result={}
        ags=find_nearby_chains (pdbid,h,l,mol_info,hvnum,lvnum,pdb_dir) # The chain ID(s) of possible antigens
        #ags=old_find_nearby_chains (pdbid,h,l,hvnum,lvnum,pdb_dir)
        if len(ags)!=0:
            for ag in ags:
                sub_result[ag]=mol_info[ag]
            """
            If there are multiple molecule types in sub_result,
            and one of them is polypeptide, then remove other chains from antigen chain list.
            (Same for polyribonucleotide & polydeoxyribonucleotide)
            """
            mol_types=set([i[1] for i in sub_result.values()])

            if len(mol_types)!=1:
                for i in ["polypeptide(L)","polyribonucleotide","polydeoxyribonucleotide"]:
                    if i in mol_types:
                        sub_result={k:v for k,v in sub_result.items() if v[1]==i}

        if sub_result!={}:
            result.append(sub_result)

    return result

### 1.3 Collect all the information
def gather_pdbeinfo_for_all(o_df,pdb_dir):
    # pdb_dir: the directory of full pdb
    df=o_df.copy()

    title=[]
    date=[]
    method=[]
    resolution=[]
    sugar=[]
    pdb_Hspecies=[]
    pdb_Lspecies=[]

    antigen_chain=[]
    antigen_type=[]
    antigen_description=[]

    error=[]

    for i in df.index:
        pdb=df.loc[i,"pdb"]
        Hchains=df.loc[i,"Hchain"]
        Lchains=df.loc[i,"Lchain"]
        iden_code=df.loc[i,"iden_code"]
        h=Hchains.split(";")[0]
        l=Lchains.split(";")[0]

        # This section is just to try to get info from pdbe multiple times, in case the connection time exceed and everything will be lost
        try:
            pdbinfo=get_info_from_pdbe(pdb)
        except:
            try:
                pdbinfo=get_info_from_pdbe(pdb)
            except:
                pdbinfo=get_info_from_pdbe(pdb)
                print ("connection time exceed: can't try anymore")

        mol_info=pdbinfo["mol_info"]

        title.append(pdbinfo['title'])
        date.append(pdbinfo['release_date'])
        method.append(pdbinfo['method'])
        resolution.append(pdbinfo['resolution'])
        sugar.append(pdbinfo['carbohydrate'])

        # Get species info
        hspecies=""
        lspecies=""
        if h in mol_info.keys():
            hspecies=mol_info[h][4]
        if l in mol_info.keys():
            lspecies=mol_info[l][4]
        pdb_Hspecies.append(hspecies)
        pdb_Lspecies.append(lspecies)

        # Get antigen info
        ag=[] # antigen chain name
        ag_t=[] # antigen chain type
        ag_d=[] # antigen chain description (the mol name)

        try:
            ag_info=get_antigen_info(pdb,Hchains,Lchains,mol_info,pdb_dir)

            if len(ag_info) != 0:
                for d in ag_info:
                    ag.append("".join(list(d.keys())))
                    ag_d+=[i[0] for i in d.values()]
                    ag_t+=[i[1] for i in d.values()]
        except:
            error.append([iden_code])

        ag=";".join(ag)
        ag_d=",".join(set(ag_d))
        ag_t=",".join(set(ag_t))

        antigen_chain.append(ag)
        antigen_type.append(ag_t)
        antigen_description.append(ag_d)



    df['title']=title
    df['release_date']=date
    df['method']=method
    df['resolution']=resolution
    df['carbohydrate']=sugar

    df["H_pdb_species"]=pdb_Hspecies
    df["L_pdb_species"]=pdb_Lspecies

    df["antigen_chain"]=antigen_chain
    df["antigen_type"]=antigen_type
    df["antigen_description"]=antigen_description

    return df,error

############ 2. Identify species and chain types ###########

# 2.1 Functions to perform BLAST

# Generate seq fasta file (VCAb seqs) for the blast purpose
def convert_seq_from_df_to_fasta (df,col_name,out_dir):
    # extract the seq column, then convert it into a fasta file
    file=open(f"{out_dir}/{col_name}.fasta", mode = 'w')
    for i in df.index:
        full_name='>'+df.loc[i,'iden_code']+'-'+col_name
        seq=df.loc[i,col_name]
        if type(seq)==str:
            file.write(full_name+'\n')
            file.write(seq+'\n')
    file.close()

def generate_bl_result (q_seqs,bl_db,bl_out_name,out_dir,out_format=10):
    # for out_format, csv output is 10
    # creat the blastp command line
    bl_cmd=ncbiblp(query=q_seqs,db=bl_db,max_target_seqs=50,outfmt=out_format,out=f"{out_dir}/{bl_out_name}_result.csv") # the output file would be the hit table

    # execute the blastp command line & generate the result.csv (the hit table)
    stdout,stderr=bl_cmd()
    header=['qseqid','sseqid','ident','length','mismatches','gapopen','start_ab','end_ab','start_ref','end_ref','e-value','score']
    bl_result=pd.read_csv(f"{out_dir}/{bl_out_name}_result.csv",names=header)

    # Modify the bl_result
    df=bl_result.copy() # no specific reason, just want to type fewer characters/or in case I want to output the initial bl_result
    df["iden_code"]=list(map(lambda x: x.split("-")[0],df["qseqid"]))

    df["matched_alleles"]=[i.split("|")[0] for i in df["sseqid"]]
    df["species"]=[i.split("|")[1] for i in df["sseqid"]]
    df["chain_type"]=[get_chain_type_names (i) for i in df["matched_alleles"]]

    df.to_csv(f"{out_dir}/{bl_out_name}_result.csv")
    return df

def generate_vbl_result (q_seqs,bl_db,bl_out_name,out_dir,out_format=15):
    # creat the blastp command line
    # for out_format, csv output is 10, json(one-file) is 15
    bl_cmd=ncbiblp(query=q_seqs,db=bl_db,max_target_seqs=50,outfmt=out_format,out=f"{out_dir}/{bl_out_name}_result_{out_format}") # the output file would be the hit table

    # execute the blastp command line & generate the result.csv (the hit table)
    stdout,stderr=bl_cmd()

def filter_bl (o_bl):
    # only keep the top hits for each qseqid, each sseqid(allele)
    # Note: it is different from the get_best_hit_for_each_allele function
    result_lst=[]
    for qseqid,group in o_bl.groupby('qseqid', sort=False):
        for sseqid, a_df in group.groupby('sseqid',sort=False):
            result_lst.append(a_df.iloc[0,:])
    result=pd.DataFrame(result_lst)
    return result

# 2.2 Functions to identify species and chain types
def get_chain_type_names (o_g_name):
    # o_g_name: original gene name, e.g. IGHG1*01
    g_name=o_g_name.split("*")[0]
    if g_name[0:3]=="IGH":
        return "Ig"+g_name[3:]
    elif g_name[0:4]=="IGKC":
        return "kappa"
    elif g_name[0:4]=="IGLC":
        return "lambda"
    elif g_name[0:4]=="IGIC":
        return "iota"
    else:
        return g_name

def species_translator(s_species):
    # translate pdb_species (s_species) into the certain format(all_lower_cases, replace space with "_") which is searchable in blast_species
    result=float('nan')
    if type(s_species)==str:
        result=[]
        s_species_lst=s_species.split(", ") # in case there are multiple species in the listed species.
        for i0 in s_species_lst:
            for i in i0.split(","):
                sub_result=i.replace(" ","_")
                sub_result=sub_result.lower()
                if ("virus" not in sub_result) and ("synthetic" not in sub_result) and ('staphylococcus_aureus' not in sub_result) and ("escherichia_coli" not in sub_result) and ("vector" not in sub_result) and ("other" not in sub_result):
                    result.append(sub_result)
    return result

def include_collapsed_alleles(species,allele,o_c_alleles_dict):
    # Note: the returned string doesn't include the species name
    c_alleles_dict={k.replace(" ","_").lower():v for k,v in o_c_alleles_dict.items()}
    all_alleles=[allele]
    all_alleles+=c_alleles_dict[species][allele]
    return ",".join(all_alleles)

def extract_info_from_json(vbl_json_fn):
    # vbl_json_fn is the file name of the v_blast json file
    # Two purpose:
    ## 1. read the json file, then convert the json to dict
    ## 2. Only take the top hits for each allele

    dct=json.load(open(vbl_json_fn,'r'))
    result_lst=[]
    for query in dct['BlastOutput2']:
        query_id=query['report']["results"]['search']["query_title"]
        total_hits=query['report']["results"]['search']["hits"]

        for hit in total_hits:
            query_info={}
            query_info["query_id"]=query_id
            query_info["subject_id"]=hit["description"][0]["title"]
            query_info.update(hit["hsps"][0]) # only take the top hits for each allele
            query_info_df=pd.DataFrame(query_info,index=[0])

            result_lst.append(query_info_df)

    result=pd.concat(result_lst).reset_index(drop=True)
    result["matched_alleles"]=[i.split("|")[0] for i in result["subject_id"]]
    result["species"]=[i.split("|")[1] for i in result["subject_id"]]
    result["iden_code"]=[i.split("-")[0] for i in result["query_id"]]
    return result

#### Identify the V region species ########
# Possibility: species, humanized, or chimera
def old_extract_CDRs_position (seq,HorL):
    """
    extract cdr loop positions using abnum
    """
    # Extract CDR1, 2 from the given H_seq/L
    # HorL can only be "H" or "L"

    out=[]
    # CDR loops defination (H chain) from IMGT

    cdr1_pat= re.compile(f'{HorL}27 .*\n{HorL}39', re.M|re.S) # CDR1 definition from IMGT
    cdr2_pat= re.compile(f'{HorL}56 .*\n{HorL}66', re.M|re.S) # CDR1 definition from IMGT
    #cdr3_pat= re.compile(f'{HorL}105 .*\n{HorL}118', re.M|re.S) # CDR3 definition from IMGT

    abnum = requests.get('http://www.bioinf.org.uk/abs/abnum/abnum.cgi?plain=1&aaseq=' + seq + '&scheme=-imgt')


    def pat_finder (pat):
        for match in pat.finditer(abnum.content.decode("utf-8", "ignore")):
            ## .decode("utf-8", "ignore") is used to convert bytes-like object to string
            info = match.group(0).split('\n')
            cdr = ''.join([i[-2:].strip() for i in info[:-1]])
            return cdr
    cdr1= pat_finder (cdr1_pat)
    cdr2= pat_finder (cdr2_pat)
    #cdr3= pat_finder (cdr3_pat)

    # Find the CDR positions
    cdr1_s=seq.find(cdr1)
    cdr1_e=cdr1_s+len(cdr1) # residue at end position is not included

    cdr2_s=seq.find(cdr2)
    cdr2_e=cdr2_s+len(cdr2)

    #cdr3_s=seq.find(cdr3)
    #cdr3_e=cdr3_s+len(cdr3)


    return [(cdr1_s,cdr1_e),(cdr2_s,cdr2_e)]

def extract_CDRs_position (pdb,chain,num_df):
    """
    extract cdr loop positions using the numbering result outputed by ANARCI
    :args pdb: pdbid
    :args chain: chain id
    :args num_df_fn: file name of df containing numbering info outputed by ANARCI, hvn or lvn

    """

    out=[]
    # CDR loops defination (H chain) from IMGT

    id_code=f"{pdb}_{chain}"
    num_info=num_df.loc[num_df["Id"]==id_code]

    cdr1= "".join([i for i in list(num_info.loc[:,"27":"38"].values[0]) if i!="-" and i!="deleted"])
    # IMGT CDR1 position:from 27 to 38, including both ends
    cdr2= "".join([i for i in list(num_info.loc[:,"56":"65"].values[0]) if i!="-" and i!="deleted"])
    #IMGT CDR2 position:from 56 to 65, including both ends

    # Find the CDR positions
    seq="".join([i for i in list(num_info.loc[:,"1":].values[0]) if i!="-" and i!="deleted"])
    cdr1_s=seq.find(cdr1)
    cdr1_e=cdr1_s+len(cdr1) # residue at end position is not included

    cdr2_s=seq.find(cdr2)
    cdr2_e=cdr2_s+len(cdr2)

    return [(cdr1_s,cdr1_e),(cdr2_s,cdr2_e)]

def map_pure_seq_pos_to_aln_pos (pure_pos,aln_seq):
    # return the aln_pos of the corresponding pure_seq
    # aln_seq: the seq from the aln object, which containing the gaps
    # both the numbering of pure_pos and aln_seq starts from 0
    pure_seq_counter=0
    for aln_pos,r in enumerate(aln_seq):
        if r!="-":
            if pure_seq_counter==pure_pos:
                return aln_pos
            pure_seq_counter+=1
def find_alternative_vtype (bl_df,b_all_alleles,b_species,bident,bscore,c_alleles):
    """
    find alternative chain type (C region)
    :args b_species & b_all_alleles: the assigned best species and alleles(including the collapsed ones)
    :args bl_df: only contains the filtered df of specific iden_code
    :args bident & bscore: the ident and score of the best hit
    :args c_alleles: collapsed_c_alleles
    """

    result_lst=[]
    s_df=bl_df
    #s_df=bl_df.loc[bl_df["species"]==species] # only contains specific species
    b_all_alleles=b_all_alleles.split(",")
    new_sseqid=[f"{i}|{b_species.capitalize()}" for i in b_all_alleles]
    alternative_hits=s_df.loc[map(lambda a,b,c: (a==bident or b==bscore) and c not in new_sseqid,s_df['identity'],s_df["score"],s_df['subject_id'])]


    if len(alternative_hits)!=0:
        for i in alternative_hits.index:
            a_ident=alternative_hits.loc[i,'identity']
            a_alleles=alternative_hits.loc[i,"matched_alleles"]
            a_species=alternative_hits.loc[i,"species"].lower()
            a_all_alleles=include_collapsed_alleles(a_species,a_alleles,c_alleles)
            result_lst.append(f"{a_species}|{a_all_alleles}: Per.Ident: {a_ident}")
    result=";".join(result_lst)
    return result

def determine_v_species_type (v_bl, iden_code, o_pdb_species, HorL,c_alleles,cut_off,species_limit=None):
    # v_bl is the df transfered from the json file, including only the top hit for each allele.
    # o_pdb_species is the original PDB species, c_alleles: collapsed alleles
    # c_species: a list contaning the species annotation for both CH1 & CL domain
    # Only use this function when PDB_species and blast_species is different

    # Assemble the df contanining the best hit (the one with best identity) for each species
    bl_df=v_bl.loc[v_bl["iden_code"]==iden_code]
    best_hit_total_lst=[] # Contains the hit with best identity for each species
    for species,group in bl_df.groupby("species",sort=False):
        best_hit=group.sort_values(by=["identity","evalue","score"],ascending=[False,True,False]).iloc[0,:]
        best_hit_total_lst.append(best_hit)


    best_hit_total=pd.DataFrame(best_hit_total_lst).sort_values(by=["identity","evalue","score"],ascending=[False,True,False])
    best_hit_total["species"]=[i.lower() for i in best_hit_total["species"]] # change the species name into lower cases
    if species_limit:
        best_hit_total=best_hit_total.loc[best_hit_total["species"]==species_limit].reset_index(drop=True)
    #return best_hit_total


    # Get bl_hits info
    bl_best_hit=best_hit_total.iloc[0,:]
    bl_species=bl_best_hit["species"]
    bl_identity=bl_best_hit["identity"]
    bl_score=bl_best_hit["score"]
    bl_alleles=bl_best_hit["matched_alleles"]
    bl_all_alleles=include_collapsed_alleles(bl_species,bl_alleles,c_alleles)
    bl_alternative_v=find_alternative_vtype (bl_df,bl_all_alleles,bl_species,bl_identity,bl_score,c_alleles)

    # Get pdb_hits info
    pdb_hits_lst=[]
    pdb_species_lst=species_translator(o_pdb_species)
    if type(pdb_species_lst)!=list:
        # pdb_species_lst could be float('nan'), if there is no species info acceessed from PDBeAPI
        return (bl_species,f"{bl_species}|{bl_all_alleles}: Per.Ident: {bl_identity}",bl_alternative_v)

    for i in pdb_species_lst:
        pdb_hit=best_hit_total.loc[best_hit_total["species"]==i]
        pdb_hits_lst.append(pdb_hit)

    if pdb_hits_lst==[]:
        # For cases like "synthetic construc" in pdb species
        return (bl_species,f"{bl_species}|{bl_all_alleles}: Per.Ident: {bl_identity}",bl_alternative_v)

    pdb_total_hits=pd.concat(pdb_hits_lst).sort_values(by=["identity","evalue","score"],ascending=[False,True,False])
    # If the pdb_species not in the bl_species (i.e. pdb_total_hits is empty), then use the blast species
    if len(pdb_total_hits)==0:
        return (bl_species,f"{bl_species}|{bl_all_alleles}: Per.Ident: {bl_identity}",bl_alternative_v)

    pdb_best_hit=pdb_total_hits.iloc[0,:]
    pdb_identity=pdb_best_hit["identity"]
    pdb_species=pdb_best_hit["species"]

    pdb_alleles=pdb_best_hit["matched_alleles"]

    if pdb_species == bl_species:
        return (bl_species,f"{bl_species}|{bl_all_alleles}: Per.Ident: {bl_identity}",bl_alternative_v)
    else:
        # Step 1. Define the first_species
        first_species=pdb_species

        """
        overwrite pdb_species with bl_species when:
        1. bl_identity > pdb_identity+cut_off
        2. bl_identity==100
        #3. bl_species(V region) is the same as bl_species(C region) & bl_identity > pdb_identity
        #   i.e. in this case, the cut_off in 1. is set to 0
        """
        if bl_identity > pdb_identity+cut_off:
            first_species=bl_species
        if bl_identity==100:
            first_species=bl_species
        #if bl_species in c_species:
        #    if bl_identity > pdb_identity:
        #        first_species=bl_species
        if species_limit:
            first_species=species_limit

        first_species_hit=best_hit_total.loc[best_hit_total["species"]==first_species]
        first_species_alleles=first_species_hit["matched_alleles"].values[0]
        first_species_all_alleles=include_collapsed_alleles(first_species,first_species_alleles,c_alleles)

        first_species_identity=first_species_hit["identity"].values[0]
        first_species_score=first_species_hit["score"].values[0]

        first_species_alternative_v=find_alternative_vtype (bl_df,bl_all_alleles,bl_species,bl_identity,bl_score,c_alleles)

        # Step 2. If the first_species is human, check if it is humanized;
        ## if it is not, return the first_species

        if first_species=="homo_sapiens":
            hit_to_be_compared=pd.DataFrame([bl_best_hit,pdb_best_hit])
            human_hit=hit_to_be_compared.loc[hit_to_be_compared["species"]=="homo_sapiens"]
            non_human_hit=hit_to_be_compared.loc[hit_to_be_compared["species"]!="homo_sapiens"]

            human_hit_qseq=human_hit["qseq"].values[0]
            nonhuman_hit_qseq=non_human_hit["qseq"].values[0]
            #pure_seq=human_hit_qseq.replace("-","") # pure_seq for human_qseq and nonhuman_qseq are the same: they are all qseq

            # extract CDR position
            pdbid=iden_code.split("_")[0]
            chainid=iden_code.split("_")[1][0]
            num_df=hvnum
            if HorL=="L":
                chainid=iden_code.split("_")[1][1]
                num_df=lvnum
            cdr1,cdr2=extract_CDRs_position (pdbid,chainid,num_df)

            # Calculate the position of "midline" corresponding to the CDR fragment
            ## human_hits_aln_position:
            hhit_cdr1s=map_pure_seq_pos_to_aln_pos (cdr1[0],human_hit_qseq)
            hhit_cdr1e=map_pure_seq_pos_to_aln_pos (cdr1[1],human_hit_qseq)
            hhit_cdr2s=map_pure_seq_pos_to_aln_pos (cdr2[0],human_hit_qseq)
            hhit_cdr2e=map_pure_seq_pos_to_aln_pos (cdr2[1],human_hit_qseq)

            ## non_human_hits_aln_position:
            nhhit_cdr1s=map_pure_seq_pos_to_aln_pos (cdr1[0],nonhuman_hit_qseq)
            nhhit_cdr1e=map_pure_seq_pos_to_aln_pos (cdr1[1],nonhuman_hit_qseq)
            nhhit_cdr2s=map_pure_seq_pos_to_aln_pos (cdr2[0],nonhuman_hit_qseq)
            nhhit_cdr2e=map_pure_seq_pos_to_aln_pos (cdr2[1],nonhuman_hit_qseq)

            # Calculate the matches% of CDR1&2
            ## human_hits_midline:
            hm=human_hit["midline"].values[0]
            hm_cdr1=hm[hhit_cdr1s:hhit_cdr1e]
            hm_cdr1_match=len([i for i in hm_cdr1 if (i!=" ") and (i != "+")])
            hm_cdr2=hm[hhit_cdr2s:hhit_cdr2e]
            hm_cdr2_match=len([i for i in hm_cdr2 if (i!=" ") and (i != "+")])

            hm_cdr12_match=(hm_cdr1_match+hm_cdr2_match)/(len(hm_cdr1)+len(hm_cdr2))

            ## non_human_hits_midline:
            nhm=non_human_hit["midline"].values[0]
            nhm_cdr1=nhm[nhhit_cdr1s:nhhit_cdr1e]
            nhm_cdr1_match=len([i for i in nhm_cdr1 if (i!=" ") and (i != "+")])
            nhm_cdr2=nhm[nhhit_cdr2s:nhhit_cdr2e]
            nhm_cdr2_match=len([i for i in nhm_cdr2 if (i!=" ") and (i != "+")])
            nhm_cdr12_match=(nhm_cdr1_match+nhm_cdr2_match)/(len(nhm_cdr1)+len(nhm_cdr2))

            if nhm_cdr12_match > hm_cdr12_match:
                #return (f"Humanized:{pdb_species};{bl_species}",f"{pdb_species}|{pdb_alleles}: Per.Ident: {pdb_identity};{bl_species}|{bl_alleles}: Per.Ident: {bl_identity}",bl_alternative_v)
                return (f"Humanized",f"{first_species}|{first_species_alleles}: Per.Ident: {first_species_identity}",bl_alternative_v)
            else:
                # first_species is human in this case
                return (first_species,f"{first_species}|{first_species_alleles}: Per.Ident: {first_species_identity}",first_species_alternative_v)



        else:
            return (first_species,f"{first_species}|{first_species_alleles}: Per.Ident: {first_species_identity}",first_species_alternative_v)

#### Identify the C region species ######
def find_alternative_ctype (bl_df,b_all_alleles,b_species,bident,bscore,c_alleles):
    """
    find alternative chain type (C region)
    :args b_species & b_all_alleles: the assigned best species and alleles(including the collapsed ones)
    :args bl_df: only contains the filtered df of specific iden_code
    :args bident & bscore: the ident and score of the best hit
    :args c_alleles: collapsed_c_alleles
    """

    result_lst=[]
    s_df=bl_df
    #s_df=bl_df.loc[bl_df["species"]==species] # only contains specific species
    b_all_alleles=b_all_alleles.split(",")
    new_sseqid=[f"{i}|{b_species.capitalize()}" for i in b_all_alleles]
    alternative_hits=s_df.loc[map(lambda a,b,c: (a==bident or b==bscore) and c not in new_sseqid,s_df["ident"],s_df["score"],s_df["sseqid"])]


    if len(alternative_hits)!=0:
        for i in alternative_hits.index:
            a_ctype=alternative_hits.loc[i,"chain_type"]
            a_ident=alternative_hits.loc[i,"ident"]
            a_alleles=alternative_hits.loc[i,"matched_alleles"]
            a_species=alternative_hits.loc[i,"species"].lower()
            a_all_alleles=include_collapsed_alleles(a_species,a_alleles,c_alleles)
            result_lst.append(f"{a_ctype}({a_species}|{a_all_alleles}: Per.Ident: {a_ident})")
    result=";".join(result_lst)
    return result

def determine_c_species_type (c_bl, iden_code, o_pdb_species, c_alleles,cut_off,species_limit=None):
    # c_bl is the bl_df contaning only the top_hits for each allele.
    # o_pdb_species is the original species, c_alleles: collapsed alleles
    # Only use this function when PDB_species and blast_species is different

    # Assemble the df contanining the best hit (the one with best identity) for each species
    bl_df=c_bl.loc[c_bl["iden_code"]==iden_code]
    best_hit_total_lst=[] # Contains the hit with best identity,e-value,and score for each species
    for species,group in bl_df.groupby("species",sort=False):
        best_hit=group.sort_values(by=["ident","e-value","score"],ascending=[False,True,False]).iloc[0,:]
        best_hit_total_lst.append(best_hit)


    best_hit_total=pd.DataFrame(best_hit_total_lst).sort_values(by=["ident","e-value","score"],ascending=[False,True,False])

    best_hit_total["species"]=[i.lower() for i in best_hit_total["species"]] # change the species name into lower cases
    if species_limit:
        best_hit_total=best_hit_total.loc[best_hit_total["species"]==species_limit].reset_index(drop=True)

    #return best_hit_total
    # Get bl_hits info
    bl_best_hit=best_hit_total.iloc[0,:]
    bl_species,bl_identity,bl_score,bl_ctype,bl_alleles,bl_vcb=bl_best_hit[["species","ident","score","chain_type","matched_alleles","start_ab"]]

    # Get pdb_hits info
    pdb_hits_lst=[]
    pdb_species_lst=species_translator(o_pdb_species)
    if type(pdb_species_lst)!=list:
        # pdb_species_lst could be float('nan'), if there is no species info acceessed from PDBeAPI
        bl_all_alleles=include_collapsed_alleles(bl_species,bl_alleles,c_alleles)
        bl_alternative_ctype=find_alternative_ctype (bl_df,bl_all_alleles,bl_species,bl_identity,bl_score,c_alleles)
        return (bl_species,f"{bl_ctype}({bl_species}|{bl_all_alleles}: Per.Ident: {bl_identity})", bl_alternative_ctype, bl_vcb)

    for i in pdb_species_lst:
        pdb_hit=best_hit_total.loc[best_hit_total["species"]==i]
        pdb_hits_lst.append(pdb_hit)

    if pdb_hits_lst==[]:
        # For cases like "synthetic construc" in pdb species
        bl_all_alleles=include_collapsed_alleles(bl_species,bl_alleles,c_alleles)
        bl_alternative_ctype=find_alternative_ctype (bl_df,bl_all_alleles,bl_species,bl_identity,bl_score,c_alleles)
        return (bl_species,f"{bl_ctype}({bl_species}|{bl_all_alleles}: Per.Ident: {bl_identity})", bl_alternative_ctype, bl_vcb)

    pdb_total_hits=pd.concat(pdb_hits_lst).sort_values(by=["ident","e-value","score"],ascending=[False,True,False])
    # If the pdb_species not in the bl_species (i.e. pdb_total_hits is empty), then use the blast species
    if len(pdb_total_hits)==0:
        bl_all_alleles=include_collapsed_alleles(bl_species,bl_alleles,c_alleles)
        bl_alternative_ctype=find_alternative_ctype (bl_df,bl_all_alleles,bl_species,bl_identity,bl_score,c_alleles)
        return (bl_species,f"{bl_ctype}({bl_species}|{bl_all_alleles}: Per.Ident: {bl_identity})", bl_alternative_ctype, bl_vcb)

    pdb_best_hit=pdb_total_hits.iloc[0,:]
    pdb_species=pdb_best_hit["species"]
    pdb_identity=pdb_best_hit["ident"]

    if pdb_species == bl_species:
        bl_all_alleles=include_collapsed_alleles(bl_species,bl_alleles,c_alleles)
        bl_alternative_ctype=find_alternative_ctype (bl_df,bl_all_alleles,bl_species,bl_identity,bl_score,c_alleles)
        return (bl_species,f"{bl_ctype}({bl_species}|{bl_all_alleles}: Per.Ident: {bl_identity})", bl_alternative_ctype, bl_vcb)
    else:
        # Define the first_species (pdb_species by default, unless the Id% for bl_species is 10% higher)
        first_species=pdb_species
        """
        overwrite pdb_species with bl_species when:
        1. bl_identity > pdb_identity+cut_off
        2. bl_identity==100
        """
        if bl_identity > pdb_identity+cut_off:
            first_species=bl_species
        if bl_identity==100:
            first_species=bl_species

        first_species_hit=best_hit_total.loc[best_hit_total["species"]==first_species]
        first_species_identity=first_species_hit["ident"].values[0]
        first_species_score=first_species_hit["score"].values[0]
        first_species_ctype=first_species_hit["chain_type"].values[0]

        first_species_alleles=first_species_hit["matched_alleles"].values[0]
        first_species_all_alleles=include_collapsed_alleles(first_species,first_species_alleles,c_alleles)
        first_species_vcb=first_species_hit["start_ab"].values[0]
        first_species_alternative_ctype=find_alternative_ctype (bl_df,first_species_all_alleles,first_species,first_species_identity,first_species_score,c_alleles)

        return (first_species,f"{first_species_ctype}({first_species}|{first_species_all_alleles}: Per.Ident: {first_species_identity})",first_species_alternative_ctype, first_species_vcb)

def add_all_domain_species_type (o_df,vhbl,vlbl,hbl,lbl,c_vh,c_vl,c_ch,c_cl,c_v_alleles,c_c_alleles):
    # vh_bl,vl_bl,hbl,lbl: filtered blast result for vh,vl,ch,cl.
    # c_vh,c_vl,c_ch,c_cl: cut-off for vh,vl,ch,cl
    # c_v_alleles,c_c_alleles: collapsed v alleles & c alleles
    df=o_df.copy()

    vh_assigned_species=[]
    vh_assigned_allele=[]
    vl_assigned_species=[]
    vl_assigned_allele=[]

    ch_assigned_species=[]
    ch_assigned_allele=[]
    cl_assigned_species=[]
    cl_assigned_allele=[]

    h_vcb=[]
    l_vcb=[]
    al_htype=[]
    al_ltype=[]

    al_vhtype=[]
    al_vltype=[]

    dropped_index=[]

    for i in df.index:
        iden_code=df.loc[i,"iden_code"]
        h_pdb_species=df.loc[i,'H_pdb_species']
        l_pdb_species=df.loc[i,'L_pdb_species']
        h_seq=df.loc[i,"H_seq"]
        l_seq=df.loc[i,"L_seq"]

        try:
        #if True:
            s_ch_species,s_htype,s_al_htype,s_hvcb=determine_c_species_type (hbl,iden_code, h_pdb_species,c_c_alleles,c_ch)
            s_cl_species,s_ltype,s_al_ltype,s_lvcb=determine_c_species_type (lbl,iden_code, l_pdb_species,c_c_alleles,c_cl)

            s_vh_species,s_vh_alleles,s_al_vhtype=determine_v_species_type (vhbl, iden_code, h_pdb_species, "H",c_v_alleles,c_vh)
            s_vl_species,s_vl_alleles,s_al_vltype=determine_v_species_type (vlbl, iden_code, l_pdb_species, "L",c_v_alleles,c_vl)



            vh_assigned_species.append(s_vh_species)
            vh_assigned_allele.append(s_vh_alleles)
            vl_assigned_species.append(s_vl_species)
            vl_assigned_allele.append(s_vl_alleles)

            ch_assigned_species.append(s_ch_species)
            ch_assigned_allele.append(s_htype)
            al_htype.append(s_al_htype)
            h_vcb.append(s_hvcb)

            cl_assigned_species.append(s_cl_species)
            cl_assigned_allele.append(s_ltype)
            al_ltype.append(s_al_ltype)
            l_vcb.append(s_lvcb)

            al_vhtype.append(s_al_vhtype)
            al_vltype.append(s_al_vltype)

        except:
            print (iden_code)
            dropped_index.append(i)

    df=df.drop(index=dropped_index)
    df["VH_species"]=vh_assigned_species
    df["VL_species"]=vl_assigned_species
    df["HC_species"]=ch_assigned_species
    df["LC_species"]=cl_assigned_species

    df["VH_allele"]=vh_assigned_allele
    df["VL_allele"]=vl_assigned_allele
    df["Htype"]=ch_assigned_allele
    df["Ltype"]=cl_assigned_allele
    df["Alternative_VH_allele"]=al_vhtype
    df["Alternative_VL_allele"]=al_vltype
    df["Alternative_Htype"]=al_htype
    df["Alternative_Ltype"]=al_ltype


    df["H_seq_VC_boundary"]=h_vcb
    df["L_seq_VC_boundary"]=l_vcb

    return df


############ 3. Identify structural coverage ###########
# Collect the domain information for the constant region from the title of reference sequences
def extract_domain_information_from_imgt_fasta(fn, ifH):
    #fn is a fasta file downloaded directly from IMGT, covering all the sequences from all the species.
    # ifH is an integer, indicating if CH or CL will be extracted, with 1 being yes (H chain), 0 being no (L chain)

    # Step1. Convert the fasta tilte into df containing the sequence information.
    all_title=[]
    for record in SeqIO.parse(fn,"fasta"):
        title_lst=record.description.split("|")

        all_title.append(title_lst)
    fasta_titles=pd.DataFrame(all_title,columns=["IMGT_accession","allele_name","species","functionality","domains","st_end_position","nt_length","codon_start","9","10","11","aa_length","13","partial","15","16"]).drop(columns=["16"])
    ### extract the "pure_species" info, then replace the space in "pure_species" with "_"
    fasta_titles["pure_species"]=[i.split("_")[0].replace(" ","_") for i in fasta_titles["species"]]
    n_fasta_titles=fasta_titles.loc[(fasta_titles["functionality"]=="F") | (fasta_titles["functionality"]=="(F)")]
    f_fasta_titles1=n_fasta_titles.loc[map(lambda x:("M" not in x) and ("-LIKE" not in x), n_fasta_titles["domains"])]
    f_fasta_titles=f_fasta_titles1.loc[map(lambda x: x[0:2]=="IG",f_fasta_titles1["allele_name"])]



    # Now, the left are: V-region, D-region, J-region, C-region, and other CH domains
    ## only include C_genes:
    c_titles=f_fasta_titles.loc[map(lambda x: (x!="V-REGION") and (x!="D-REGION") and (x!="J-REGION"),f_fasta_titles["domains"])]
    c_titles["if_H"]=[int(i[0:3]=="IGH") for i in c_titles["allele_name"]]

    # Filter only for the HC or LC sequence information
    flt_c_titles=c_titles.loc[c_titles["if_H"]==ifH]

    domains={}
    result=[]
    for i in flt_c_titles.index:
        a_name=flt_c_titles.loc[i,"allele_name"]
        species=flt_c_titles.loc[i,"pure_species"].lower()
        allele_name=f"{species}|{a_name}"

        if_partial=flt_c_titles.loc[i,"partial"]# if the allele sequence only contains a fragment
        domain_name=flt_c_titles.loc[i,"domains"]
        domain_seq_len=int(flt_c_titles.loc[i,"aa_length"][0:-3])

        if allele_name not in domains.keys():
            d_counter=domain_seq_len
            domains[allele_name]={domain_name:[(1,domain_seq_len),d_counter,if_partial]}

            result.append([1,domain_seq_len,allele_name,domain_name])
        else:
            last_val=list(domains[allele_name].values())[-1]
            previous_d_counter=last_val[1]
            d_counter=domain_seq_len
            domains[allele_name][domain_name]=[(1+previous_d_counter,domain_seq_len+previous_d_counter),d_counter+previous_d_counter,if_partial]
            result.append([1+previous_d_counter,domain_seq_len+previous_d_counter,allele_name,domain_name])
    result_df=pd.DataFrame(result,columns=["a","b","q","dom"])
    return domains,result_df


imgt_original="../seq_db/ref_db/all_species_ref_seq_imgt.fasta"
h_domains,h_domains_df=extract_domain_information_from_imgt_fasta(imgt_original,1)
l_domains,l_domains_df=extract_domain_information_from_imgt_fasta(imgt_original,0)

with open("h_domains.json","w") as outfile:
    json.dump(h_domains,outfile)
with open("l_domains.json","w") as outfile:
    json.dump(l_domains,outfile)

h_domains_df.to_csv("h_domains_info.csv")
l_domains_df.to_csv("l_domains_info.csv")

def aligned_to (start,end,d_name,allele,domains=h_domains):
    # return the aligned coverage for certain domain in specific Ig
    # start & end: start and end positions of the aligned region in the ref_seq in the blast aln result
    # d_name:domain name ('CH1',etc.)
    # allele: allele_name (e.g. homo_sapiens|IGHG1*01)

    d_region=domains[allele][d_name][0] #domain boundary

    a=d_region[0] # d_region start
    b=d_region[1] # d_region end
    overlapped_region=()
    overlapped_perc=0
    if b<start or end<a:
        return 'no overlap'
    else:
        ordered=sorted([a,b,start,end])
        med_small=ordered[1]
        med_large=ordered[2]

        overlapped_region=(med_small,med_large)
        overlapped_perc_domains=((med_large-med_small)/(b-a))*100
        overlapped_perc_protein=((med_large-med_small)/(end-start))*100
        return [d_name,overlapped_perc_domains,overlapped_perc_protein]

def get_struc_cov_coor_VCB (iden_code,ctype_info,df,h_domains):
    # Note: df must be the filtered (not filtered should also be fine) bl result of the true seqs (coordinate seqs)
    # return: 1. align_info; 2.Structural Coverage; 3.VC_Boundary of the coor sequence
    # df: the bl_result of the coordinate sequence
    # iden_code: pdb_HL
    # ctype_info: the chain type information directly from vcab["Htype","Ltype"],e.g.IgG1(mus_musculus|IGHG1*02,IGHG1*03: Per.Ident: 99.01)

    ig=ctype_info.split("(")[0] #ig: the isotype/light chain type of the antibody
    allele_name=ctype_info.split("(")[1].split(":")[0].split(",")[0] # The allele name with species & without collapsed alleles, format consistant with h_domains
    p_a_name=allele_name.split("|")[1] # The pure allele name without species
    species=allele_name.split("|")[0].capitalize() # The pure species name with the first letter capitalized.
    reordered_allele=f"{p_a_name}|{species}" # change the format of allele name to make it consistant with the format in blast table


    true_hit_info=df.loc[(df["iden_code"]==iden_code)&(df["species"]==species)]
    # if select the corresponding alleles in the "Htype"/"Ltype",use:
    #true_hit_info=df.loc[(df["iden_code"]==iden_code)&(df["sseqid"]==reordered_allele)]
    #hit_allele=allele_name # the matched allele of the isotype
    ## However, Abs like "5dk3_BA" wouldn't be picked up as "full antibody"
    best_true_hit_info=true_hit_info.iloc[0,]

    p_hit_allele=best_true_hit_info["matched_alleles"]
    hit_species=best_true_hit_info["species"].lower()
    hit_allele=f"{hit_species}|{p_hit_allele}"


    # Get the true_VCB info
    true_VCB=best_true_hit_info["start_ab"]

    if p_a_name[0:3]=="IGH":
        # if the chain is heavy chain
        # Get the align_info
        s,e=best_true_hit_info["start_ref"],best_true_hit_info["end_ref"]
        s_t,e_t=best_true_hit_info["start_ab"],best_true_hit_info["end_ab"]#start and end point of the aligned region in the target antibody


        regions=[i for i in h_domains[hit_allele].keys() if "M" not in i]


        align_info=[(s_t,e_t),(s,e)]
        displayed_align_info=[]

        for r in regions:
            subresult=aligned_to (s,e,r,hit_allele)
            if type(subresult)!= str:
                if subresult[1:] !=[0,0]:
                    perc1=round(subresult[1],2)
                    perc2=round(subresult[2],2)
                    sub_string=f"{subresult[0]}(covers {perc1}% of {subresult[0]}, accounts for {perc2}% of the H_coor_seq)"
                    displayed_align_info.append(sub_string)
                    align_info.append(subresult)

        # Get the struc_cov
        ig_d=[i for i in h_domains[hit_allele].keys() if "M" not in i]
        # the domain list for the ref_seq
        pro_d=[a[0] for a in align_info[2:]]
        # the domain list for this antibody
        domain_info="; ".join(displayed_align_info)
        if true_VCB > 70:
            if pro_d[0]!="CH1":
                struc_cov = 'not classified'
                # mostly: structure contains only V region and is forced to align with C region
                # other cases: 1za6 (CH2 deleted)
            elif ig_d==pro_d:
                struc_cov = 'full antibody'
            elif ("CH3" not in pro_d) and ("CH3-CHS" not in pro_d):
                struc_cov = 'Fab'
            else:
                struc_cov = 'not classified'

        else:
            struc_cov = 'check V/C Annotation' # Probably V region only
    else:
        # if the chain is light chain
        align_info,domain_info,struc_cov=("","","")

    return align_info, domain_info, struc_cov, true_VCB

def add_struc_cov(o_df,htbl,ltbl):
    df=o_df.copy()

    aln_info=[]
    domain_info=[]
    struc_cov=[]
    hvcb_true=[]
    lvcb_true=[]

    for i in df.index:
        iden_code=df.loc[i,"iden_code"]
        this_h_type=df.loc[i,"Htype"]
        this_l_type=df.loc[i,"Ltype"]

        try:
        #if True:
            this_aln_info,this_domain_info,this_struc_cov,this_ht_vcb=get_struc_cov_coor_VCB (iden_code,this_h_type,htbl,h_domains)
            __,__,__,this_lt_vcb=get_struc_cov_coor_VCB (iden_code,this_l_type,ltbl,l_domains)
        except:
            print (f"warning:struc_cov:{iden_code}")
            this_aln_info,this_struc_cov,this_ht_vcb,this_lt_vcb,this_domain_info=("unidentified","unidentified","unidentified","unidentified","unidentified")

        aln_info.append(this_aln_info)
        domain_info.append(this_domain_info)
        struc_cov.append(this_struc_cov)
        hvcb_true.append(this_ht_vcb)
        lvcb_true.append(this_lt_vcb)

    df["align_info"]=aln_info
    df["Domains in HC"]=domain_info
    df["Structural Coverage"]=struc_cov

    df["H_coordinate_seq_VC_boundary"]=hvcb_true
    df["L_coordinate_seq_VC_boundary"]=lvcb_true

    return df

############ 4. Cut mmcif files & generate PDB files ###########
def keep_chain (model,chains):
    chain_to_remove = []
    for c in model:
        if c.id not in chains:
            chain_to_remove.append(c.id)
    for c in chain_to_remove:
        model.detach_child(c)
    for c in chains:
        if len(c)!=1:
            # The length of chain id must be 1 if later the model is written out and saved
            model[c].id=c[0]
    return model

def keep_constant_domain (chain,v_c):
    # the h&l constant region starts at h_vc-1,l_vc-1 (since the numbering of python starts at 0)
    # this is why I use "c+1<v_c"
    # since I use "c" to count residues, not the residue number written in pdb, true_V_C_boundary is used for cutting C region
    residue_to_remove = []
    for c,r in enumerate(chain):
        if (c+1<v_c) or (r.id[0]!=' '): #r.id is in this format(' ',num,' '), for normal amino acid, r.id[0] is ' '
            residue_to_remove.append(r.id)
    for residue in residue_to_remove:
        chain.detach_child(residue)
    chain.detach_parent()

    if len(chain.id)!=1:
        chain.id=chain.id[0]
    return chain

def generate_C_pdb (pdb_c,h,l,h_vc,l_vc,in_dir,out_dir):
    pdb=pdb_c[0:4]



    parser = MMCIFParser()
    structure=parser.get_structure(pdb, f'{in_dir}/{pdb}.cif')
    model=structure[0]

    model_id=[c.id for c in model.get_chains()]
    if h not in model_id:
        # for chains with ID III, MMM,etc
        h=h*3
    if l not in model_id:
        l=l*3

    f_struc=Bio.PDB.Structure.Structure(f'{pdb_c}_C_struc')
    f_model=Bio.PDB.Model.Model(f'{pdb_c}_C')

    out_h=keep_constant_domain (model[h],h_vc)
    out_l=keep_constant_domain (model[l],l_vc)

    # In the pdb-file to be generated, heavy chain first, then light chain
    f_model.add(out_h)
    f_model.add(out_l)
    f_struc.add(f_model)

    io=PDBIO()
    io.set_structure(f_struc)
    io.save(f'{out_dir}/{pdb_c}_C.pdb')

def generate_C_pdb_total (df,in_dir,out_dir):
    error=[]
    for i in df.index:
        pdb_c=df.loc[i,'iden_code']
        h=df.loc[i,'Hchain'].split(";")[0]
        l=df.loc[i,'Lchain'].split(";")[0]
        
        if os.path.exists(f"{out_dir}/{pdb_c}_C.pdb"):
            continue

        try:
            h_vc=df.loc[i,'H_coordinate_seq_VC_boundary']
            l_vc=df.loc[i,'L_coordinate_seq_VC_boundary']
            generate_C_pdb (pdb_c,h,l,h_vc,l_vc,in_dir,out_dir)

        except:
            error.append(df.loc[i,:])
    return pd.DataFrame(error,columns=df.keys())

def generate_chain_pdb_total (df,in_dir,out_dir):
    error=[]
    for i in df.index:
        pdb_c=df.loc[i,'iden_code']

        if os.path.exists(f"{out_dir}/{pdb_c}.pdb"):
            #print (pdb_c)
            continue

        try:
            pdb=pdb_c[0:4]
            h=df.loc[i,'Hchain'].split(";")[0]
            l=df.loc[i,'Lchain'].split(";")[0]

            parser = MMCIFParser()
            structure=parser.get_structure(pdb, f'{in_dir}/{pdb}.cif')
            model=structure[0]

            model_id=[c.id for c in model.get_chains()]
            if h not in model_id:
                # for chains with ID III, MMM,etc
                h=h*3
            if l not in model_id:
                l=l*3

            chains=[h,l]
            new_model=keep_chain (model,chains)

            f_struc=Bio.PDB.Structure.Structure(f'{pdb_c}_struc')
            f_struc.add(new_model)

            io=PDBIO()
            io.set_structure(f_struc)
            io.save(f'{out_dir}/{pdb_c}.pdb')


        except:
            error.append(df.loc[i,:])
    return pd.DataFrame(error,columns=df.keys())

################### 5. ADD domain_seqs and PDB_VC_BOUNDARY (In order to specify the VC_Boundary in the 3D viewer of the webserver) ##############
def get_V_C_seq(df):
    all_hv=[]
    all_hc=[]
    all_lv=[]
    all_lc=[]

    for i in df.index:

        h_vc = df.loc[i,'H_seq_VC_boundary']
        l_vc = df.loc[i,'L_seq_VC_boundary']
        h_seq = df.loc[i,'H_seq']
        l_seq = df.loc[i,'L_seq']

        hv_seq = h_seq[0:h_vc]
        hc_seq = h_seq[h_vc-1:]
        lv_seq = l_seq[0:l_vc]
        lc_seq = l_seq[l_vc-1:]

        all_hv.append(hv_seq)
        all_hc.append(hc_seq)
        all_lv.append(lv_seq)
        all_lc.append(lc_seq)

    result_df=df.copy()
    result_df['HV_seq']=all_hv
    result_df['HC_seq']=all_hc
    result_df['LV_seq']=all_lv
    result_df['LC_seq']=all_lc

    return result_df

def read_coor_seq_and_the_start_position_of_chains (pdb_file_name,h,l,pdb_dir):
    # Read the coordinate sequence and the start positions of the H-L chains
    # The start position in the C region-only files would be the position of VCBoundary in PDB files
    # pdb_file_name: the pdb file name without ".pdb", is the pdb file where you are going to extract sequence from. e.g. "7c2l", "7c2l_HL", "7c2l_HL_C"
    # h, l: the chain id for H and L chains
    # pdb_dir: the directory of the pdb files

    parser = PDBParser()

    try:
        structure=parser.get_structure(pdb_file_name, f'{pdb_dir}/{pdb_file_name}.pdb')
        model=structure[0]
        try:
            hchain=model[h]
            lchain=model[l]
        except:
            # e.g.: The Hchain in vcab table is "H", but the chain name in pdb file is "h", this would cause an error
            h=h.swapcase()
            l=l.swapcase()
            hchain=model[h]
            lchain=model[l]

        hchain_start = list(hchain)[0].id[1]
        lchain_start = list(lchain)[0].id[1]
        # Read true C seq
        ## get the seqs for all the chains in the pdb file
        chains = {chain.id:seq1(''.join(residue.resname for residue in chain)).replace("X","") for chain in structure.get_chains()}

        h_seq = chains[h]
        l_seq = chains[l]
    except:
        # file not found
        hchain_start = "pdb not found"
        lchain_start = "pdb not found"
        h_seq = ""
        l_seq = ""
    return hchain_start,lchain_start,h_seq,l_seq

def read_pdb_VC_Boundary_and_C_true_seqs (df,pdb_dir):
    # df: the vcab db
    # pdb_dir: the directory of the pdb files CONTAINING C REGION ONLY
    pdb_h_vcb=[]
    pdb_l_vcb=[]
    hc_true_seqs=[]
    lc_true_seqs=[]
    error=[]
    result_df=df.copy()

    parser = PDBParser()
    for i in df.index:
        iden_code=df.loc[i,'iden_code']

        h=df.loc[i,'Hchain'][0]
        l=df.loc[i,'Lchain'][0]

        hchain_start,lchain_start,h_seq,l_seq = read_coor_seq_and_the_start_position_of_chains (f"{iden_code}_C",h,l,"../pdb_struc/c_pdb")

        pdb_h_vcb.append(hchain_start)
        pdb_l_vcb.append(lchain_start)
        hc_true_seqs.append(h_seq)
        lc_true_seqs.append(l_seq)

    result_df["pdb_H_VC_Boundary"]=pdb_h_vcb
    result_df["pdb_L_VC_Boundary"]=pdb_l_vcb

    result_df['HC_coordinate_seq']=hc_true_seqs
    result_df['LC_coordinate_seq']=lc_true_seqs

    error_df=pd.DataFrame(error,columns=df.keys())

    return result_df,error_df

########################## 6. Calculate the disulfide bond ##########################
def calc_dis_between_cys (iden_code,c_pdb_dir):
    # cat_chain: concatenated chain object of the H and L chain
    parser=PDBParser()
    structure=parser.get_structure(iden_code,f"{c_pdb_dir}/{iden_code}_C.pdb")
    h,l=iden_code.split("_")[1]
    h_obj=structure[0][h]
    l_obj=structure[0][l]

    cat_chain=Bio.PDB.Chain.Chain("n")
    for res in h_obj:
        cat_chain.add(res)
    for res in l_obj:
        id_lst=list(res.id)
        res.id=(id_lst[0],id_lst[1]+3000,id_lst[2])
        cat_chain.add(res)

    result=[]
    for counter,res1 in enumerate(cat_chain):
        res1_name=res1.resname
        for res2 in list(cat_chain)[counter+1:]:
            # to exlude the cases like CYS155(H)-CYS211(H),CYS211(H)-CYS155(H)
            res2_name=res2.resname

            if res1_name=="CYS" and res2_name=="CYS":
                diff_vector=res1["CA"].coord-res2["CA"].coord
                dist=np.linalg.norm(diff_vector)

                if dist<=7.5:
                    res1_id=res1.id[1]
                    res2_id=res2.id[1]

                    chain1=h
                    chain2=h
                    if res1_id>3000:
                        chain1=l
                        res1_id-=3000
                    if res2_id>3000:
                        chain2=l
                        res2_id-=3000

                    display_dist='%.2f' % round(dist, 2)
                    result.append(f"{res1_name}{res1_id}({chain1})-{res2_name}{res2_id}({chain2}):{display_dist}")
    return ", ".join(result)

def add_disulfide_info(df,c_pdb_dir):
    disulfide_info=[]
    for i in df.index:
        iden_code=df.loc[i,"iden_code"]
        try:
            disulfide_val=calc_dis_between_cys (iden_code,c_pdb_dir)
        except:
            disulfide_val=""
        disulfide_info.append(disulfide_val)
    df["disulfide_bond"]=disulfide_info
    return df

def assign_species_summary(o_df,vhbl,vlbl,c_v_alleles,c_vh=8,c_vl=8):
    df=o_df.copy()
    df["Alternative_VH_allele"]=[i if type(i)==str else "" for i in df["Alternative_VH_allele"]]
    df["Alternative_VL_allele"]=[i if type(i)==str else "" for i in df["Alternative_VL_allele"]]
    df["Alternative_Htype"]=[i if type(i)==str else "" for i in df["Alternative_Htype"]]
    df["Alternative_Ltype"]=[i if type(i)==str else "" for i in df["Alternative_Ltype"]]
    species=[]
    for i in df.index:

        x,y,a,b=df.loc[i,'VH_species'],df.loc[i,'VL_species'],df.loc[i,'HC_species'],df.loc[i,'LC_species']
        if len(set([x,y,a,b]))==1:
            species.append(x.title())
        elif len(set([x,y,a,b]))==2 and "Humanized" not in f"{x}{y}{a}{b}":
            #s_info=",".join(sorted(list(set([x,y,a,b]))))
            #counts=Counter([x,y,a,b])
            #majority=[(k,v) for k,v in counts.items() if v==3]
            if x==y and a==b and y!=a:
                species.append(f"Chimera")
            else:
                # Check if the two species have the same Percentage identity

                # get the alternative chain type(including species) information
                x_al,y_al,a_al,b_al=df.loc[i,"Alternative_VH_allele"],df.loc[i,"Alternative_VL_allele"],df.loc[i,"Alternative_Htype"],df.loc[i,"Alternative_Ltype"]
                # get the alternative species for domains
                x_al_s=[i.split("|")[0] for i in x_al.split(";") if i !=""]
                y_al_s=[i.split("|")[0] for i in y_al.split(";") if i !=""]
                a_al_s=[i.split("(")[1].split("|")[0] for i in a_al.split(";") if i !=""]
                b_al_s=[i.split("(")[1].split("|")[0] for i in b_al.split(";") if i !=""]


                common_species=list(set.intersection(set(x_al_s+[x]),set(y_al_s+[y]),set(a_al_s+[a]),set(b_al_s+[b])))
                if len(common_species)==1:
                    species.append(common_species[0].title())
                else:
                    species.append(f"ambiguous species: please check the species annotation for domains")


        elif ("Humanized" in x or "Humanized" in y):
            #vh_info=x.split(":")[-1].split(";")
            #vl_info=y.split(":")[-1].split(";")
            #s_info=",".join(sorted(list(set(vh_info+vl_info+[a,b]))))
            if a=="homo_sapiens" and b=="homo_sapiens":
                species.append(f"Humanized")
            else:
                # Can't be humanized:
                # Re-assign the V domain species
                c_other_species=[i for i in [a,b] if i != "homo_sapiens"]
                c_other_species_1=[i.split("_")[0] for i in c_other_species]

                iden_code=df.loc[i,"iden_code"]
                h_pdb_species=df.loc[i,"H_pdb_species"]
                l_pdb_species=df.loc[i,"L_pdb_species"]

                """
                Sometimes for the V domains identified as "humanized" (in the format of "Humanized: human, other species"),
                if the "other species" is also present in the C domains,
                it is more likely that this antibody is not "humanized" (because C domains should be human, not other species).
                Cases like this should be corrected.
                """
                if "Humanized" in x:
                    vh_other_species=[i for i in x.split(":")[1].split(";") if i != "homo_sapiens"][0]
                    if (vh_other_species in c_other_species) or (vh_other_species.split("_")[0] in c_other_species_1):
                        df.loc[i,'VH_species']=vh_other_species

                        s_vh_species,s_vh_alleles,s_al_vhtype=determine_v_species_type (vhbl, iden_code, h_pdb_species, "H",c_v_alleles,c_vh,vh_other_species)
                        df.loc[i,"VH_allele"]=s_vh_alleles
                        df.loc[i,"Alternative_VH_allele"]=s_al_vhtype

                        print (f"vh_species is changed:{iden_code}: from {x} to {vh_other_species}")
                        x=vh_other_species


                if "Humanized" in y:
                    vl_other_species=[i for i in y.split(":")[1].split(";") if i != "homo_sapiens"][0]
                    if (vl_other_species in c_other_species)  or (vl_other_species.split("_")[0] in c_other_species_1):
                        df.loc[i,'VL_species']=vl_other_species

                        s_vl_species,s_vl_alleles,s_al_vltype=determine_v_species_type (vlbl, iden_code, l_pdb_species, "L",c_v_alleles,c_vl,vl_other_species)
                        df.loc[i,"VL_allele"]=s_vl_alleles
                        df.loc[i,"Alternative_VL_allele"]=s_al_vltype

                        print (f"vl_species is changed:{iden_code}: from {y} to {vl_other_species}")
                        y=vl_other_species

                # Give species summary
                if len(set([x,y,a,b]))==1:
                    species.append(x.title())
                elif len(set([x,y,a,b]))==2 and "Humanized" not in f"{x}{y}{a}{b}":
                    #s_info=",".join(sorted(list(set([x,y,a,b]))))
                    if x==y and a==b and y!=a:
                        species.append(f"Chimera")
                    else:
                        # Check if the two species have the same Percentage identity

                        # get the alternative chain type(including species) information
                        x_al,y_al,a_al,b_al=df.loc[i,"Alternative_VH_allele"],df.loc[i,"Alternative_VL_allele"],df.loc[i,"Alternative_Htype"],df.loc[i,"Alternative_Ltype"]
                        # get the alternative species for domains
                        x_al_s=[i.split("|")[0] for i in x_al.split(";") if i !=""]
                        y_al_s=[i.split("|")[0] for i in y_al.split(";") if i !=""]
                        a_al_s=[i.split("(")[1].split("|")[0] for i in a_al.split(";") if i !=""]
                        b_al_s=[i.split("(")[1].split("|")[0] for i in b_al.split(";") if i !=""]

                        common_species=list(set.intersection(set(x_al_s+[x]),set(y_al_s+[y]),set(a_al_s+[a]),set(b_al_s+[b])))
                        if len(common_species)==1:
                            species.append(common_species[0].title())
                        else:
                            species.append(f"ambiguous species: please check the species annotation for domains")
                else:
                    species.append(f"ambiguous species: please check the species annotation for domains")


        else:
            species.append(f"ambiguous species: please check the species annotation for domains")

    df["Species"]=species

    return df

def further_classify_ambiguous(df,vh_bl,ch_bl,vl_bl,cl_bl):
    """
    1. Replace domain species annotation when:
        1.1. 3/4 domain species are the same (Domain species of four domains: A:A:A:B (3:1));
        and
        1.2. In the blast hits for the minority domain, Perc.Ident difference for majority & minority
        species is <=8.
    2. Then reassign the species summary annotation if the four domain species annotations are
    the same after modification

    df: the dataframe after the assignment of both domain species and species summary
    vh_bl,ch_bl vl_bl,cl_bl: the filtered blast table
    """
    #a_df=df.loc[df["Species"]=='ambiguous species: please check the species annotation for domains']
    df=df.copy()
    def macaca_domains(df,d_species):
        # replace the "macaca_fascicularis" with "macaca_mulatta"
        df[d_species]=["macaca_mulatta" if i=="macaca_fascicularis" else i for i in df[d_species].values]
    macaca_domains(df,"VH_species")
    macaca_domains(df,"VL_species")
    macaca_domains(df,"HC_species")
    macaca_domains(df,"LC_species")

    for i in df.index:
        check_if_3= 3 not in Counter(df.loc[i,["VH_species","VL_species","HC_species","LC_species"]].values).values()

        if df.loc[i,"Species"]!='ambiguous species: please check the species annotation for domains' or check_if_3:
            continue
            # skip the non_ambiguous ones and the ambiguous ones which doesn't have domain_species 3:1

        iden_code=df.loc[i,"iden_code"]

        vh_s,vl_s,ch_s,cl_s=df.loc[i,["VH_species","VL_species","HC_species","LC_species"]].values
        domain_species={"vh":vh_s,"vl":vl_s,"ch":ch_s,"cl":cl_s}

        species_count=Counter(domain_species.values())
        majority_species=list(species_count.keys())[list(species_count.values()).index(3)]
        minority_species=list(species_count.keys())[list(species_count.values()).index(1)]
        if "Humanized" in minority_species:
            continue
            print (iden_codes)

        minority_domain=list(domain_species.keys())[list(domain_species.values()).index(minority_species)]

        # Extract the blast hits in the minority_domain for the species: majority_species & minority_species
        bl_dfs={"vh":vh_bl,"vl":vl_bl,"ch":ch_bl,"cl":cl_bl}
        bl_df=bl_dfs[minority_domain].rename(columns={
                                                      "ident":"identity",
                                                      'qseqid':'query_id',
                                                      'sseqid':'subject_id',
                                                      'e-value':"evalue"})
        # CHECK: if this needs to be changed:
        this_bl_df=bl_df.loc[bl_df["iden_code"]==iden_code]
        this_bl_df["species"]=[i.lower() for i in this_bl_df["species"].values]
        #if True:
        try:
            major_hit=this_bl_df.loc[this_bl_df["species"]==majority_species].sort_values(by=["identity","evalue","score"],ascending=[False,True,False]).iloc[0,:]
            minor_hit=this_bl_df.loc[this_bl_df["species"]==minority_species].sort_values(by=["identity","evalue","score"],ascending=[False,True,False]).iloc[0,:]

            major_ident=major_hit["identity"]
            minor_ident=minor_hit["identity"]

            if abs(major_ident-minor_ident)<=8:

                """
                replace the minority species with majority species when the Difference in Perc.Ident for
                the two hits corresponding to these two species is smaller than 8.

                """

                # domain name: VH/VL/HC/LC
                d_name=minority_domain.upper() if minority_domain[0]=="v" else minority_domain[1].upper()+"C"

                allele=major_hit["subject_id"].split("|")[0]
                allele_info=f"{majority_species}|{allele}: Per.Ident: {major_ident}"
                if minority_domain[0]=="c":
                    chaintype=major_hit["chain_type"]
                    allele_info=f"{chaintype}({allele_info})"

                df.loc[i,d_name+"_species"]=majority_species
                df.loc[i,d_name+"_allele"]=allele_info
                print (f"{iden_code} is changed: domain: {d_name}: from {minority_species} to {majority_species}")
        except:
            print ("errors:",iden_code,minority_domain)

        # Change species annotation if four domain species are the same after domain species modification:
        x,y,a,b=df.loc[i,'VH_species'],df.loc[i,'VL_species'],df.loc[i,'HC_species'],df.loc[i,'LC_species']
        if len(set([x,y,a,b]))==1:
            df.loc[i,"Species"]=x.title()



    return df

def combine_fastas(rf,cf):
    """
    For the update of seq_db/vcab_db:
    check for the files presented in both result folder and current folder
    1. Combine the common files in rf and cf into one, stored in rf.
    :args rf: result folder: stores the old data
    :args cf: current folder: stores the updated data(only the newly added part)
    returns nothing
    """
    common_files=[fn for fn in os.listdir(rf) if os.path.exists(f"{cf}/{fn}")]
    for f in common_files:
        if ".fasta" in f:
            #print (f)

            os.system(f"cat {rf}/{f} {cf}/{f} > {cf}/combined_{f}")
            os.system(f"mv {cf}/combined_{f} {rf}/{f}")

############################ Apply the functions ############################
if __name__=="__main__":
    # Used variables:
    pHL_fn="paired_ab.csv"

    # The directories of VCAb sequences:
    hfqseqs="../seq_db/vcab_db/H_seq.fasta"
    lfqseqs="../seq_db/vcab_db/L_seq.fasta"

    htqseqs="../seq_db/vcab_db/H_coordinate_seq.fasta"
    ltqseqs="../seq_db/vcab_db/L_coordinate_seq.fasta"

    h_ref_db="../seq_db/ref_db/ch_db/all_species_unique_CH_alleles.fasta"
    l_ref_db="../seq_db/ref_db/cl_db/all_species_unique_CL_alleles.fasta"

    vh_ref_db="../seq_db/ref_db/vh_db/all_species_unique_VH_alleles.fasta"
    vl_ref_db="../seq_db/ref_db/vl_db/all_species_unique_VL_alleles.fasta"

    full_pdb_dir="../pdb_struc/full_pdb/"



    # 1. generate the collapsed paired HL antibodies containing both V and C regions
    pHL=pd.read_csv(pHL_fn).drop(columns=["Unnamed: 0"]).reset_index(drop=True)
    if "H_coor_seq" in list(pHL.columns):
        pHL=pHL.rename(columns={"H_coor_seq":"H_coordinate_seq","L_coor_seq":"L_coordinate_seq"})
    cHL1=collapse_by_coordinate_seqs (pHL) # collapse the entries based on the coordinate sequences
    cHL1.to_csv("collapsed_paired_ab.csv")
    #print ("test_this_file")
    #cHL1=pd.read_csv("collapsed_paired_ab.csv").drop(columns=["Unnamed: 0"])
    print ("collapsing completed")
    cHL,error=gather_pdbeinfo_for_all(cHL1,full_pdb_dir)
    cHL.to_csv("final_collapsed_paired_ab.csv")
    print("pdbe info collected")


    #with open(f'./pdbe_info_not_collected_lst.txt','w') as output:
    #    output.write(','.join(error))

    # 2. identify chain types
    ## generate the fasta files
    print ("Identifing chain types")
    convert_seq_from_df_to_fasta(cHL,'H_seq',"../seq_db/vcab_db")
    convert_seq_from_df_to_fasta(cHL,'H_coordinate_seq',"../seq_db/vcab_db")

    convert_seq_from_df_to_fasta(cHL,'L_seq',"../seq_db/vcab_db")
    convert_seq_from_df_to_fasta(cHL,'L_coordinate_seq',"../seq_db/vcab_db")

    ## Do the blast:
    print ("start author_seq blast")
    hbl=generate_bl_result (hfqseqs,h_ref_db,"h_seq_bl","./blast_result")
    lbl=generate_bl_result (lfqseqs,l_ref_db,"l_seq_bl","./blast_result")

    print ("start author_seq V blast")
    vhbl=generate_vbl_result (hfqseqs,vh_ref_db,"vh_seq_bl","./blast_result",15)
    vlbl=generate_vbl_result (lfqseqs,vl_ref_db,"vl_seq_bl","./blast_result",15)

    print ("start coor_seq blast")
    hcoorbl=generate_bl_result (htqseqs,h_ref_db,"h_coordinate_seq_bl","./blast_result")
    lcoorbl=generate_bl_result (ltqseqs,l_ref_db,"l_coordinate_seq_bl","./blast_result")

    ## Do the cross-blast in order to find possible domain-swapped antibodies
    print ("start cross blast")
    l_coor_cross_bl=generate_bl_result (ltqseqs,h_ref_db,"cross_bl_l_coor_seq","./blast_result")
    h_coor_cross_bl=generate_bl_result (htqseqs,l_ref_db,"cross_bl_h_coor_seq","./blast_result")

    l_author_cross_bl=generate_bl_result (lfqseqs,h_ref_db,"cross_bl_l_author_seq","./blast_result")
    h_author_cross_bl=generate_bl_result (hfqseqs,l_ref_db,"cross_bl_h_author_seq","./blast_result")

    ## Generate the filtered_blast result, which only contaning top hits for each pair of qseqid and sseqid (allele)
    print ("filter blast result")
    # For test:
    #hbl=pd.read_csv("./blast_result/h_seq_bl_result.csv").drop(columns=["Unnamed: 0"])
    #lbl=pd.read_csv("./blast_result/l_seq_bl_result.csv").drop(columns=["Unnamed: 0"])
    flt_hbl=filter_bl (hbl)
    flt_lbl=filter_bl (lbl)
    flt_hbl.to_csv("./blast_result/flt_bl_result/flt_h_seq_bl_result.csv")
    flt_lbl.to_csv("./blast_result/flt_bl_result/flt_l_seq_bl_result.csv")


    print ("filter V blast result")
    flt_vhbl=extract_info_from_json("./blast_result/vh_seq_bl_result_15")
    flt_vhbl.to_csv("./blast_result/flt_bl_result/vh_bl_15.csv")
    flt_vlbl=extract_info_from_json("./blast_result/vl_seq_bl_result_15")
    flt_vlbl.to_csv("./blast_result/flt_bl_result/vl_bl_15.csv")

    collapsed_v_alleles=json.load(open("../seq_db/ref_db/all_species_collapsed_v_alleles.json","r"))
    collapsed_c_alleles=json.load(open("../seq_db/ref_db/all_species_collapsed_c_alleles.json","r"))

    print ("identifying isotype & species")
    # For test:
    #flt_vhbl=pd.read_csv("./blast_result/flt_bl_result/vh_bl_15.csv").drop(columns=["Unnamed: 0"])
    #flt_vlbl=pd.read_csv("./blast_result/flt_bl_result/vl_bl_15.csv").drop(columns=["Unnamed: 0"])

    #flt_hbl=pd.read_csv("./blast_result/flt_bl_result/flt_h_seq_bl_result.csv").drop(columns=["Unnamed: 0"])
    #flt_lbl=pd.read_csv("./blast_result/flt_bl_result/flt_l_seq_bl_result.csv").drop(columns=["Unnamed: 0"])

    cHL=pd.read_csv("final_collapsed_paired_ab.csv").drop(columns=["Unnamed: 0"])
    iso_vcab=add_all_domain_species_type (cHL,flt_vhbl,flt_vlbl,flt_hbl,flt_lbl,8,8,8,8,collapsed_v_alleles,collapsed_c_alleles)
    iso_vcab.to_csv("iso_vcab.csv")

    print ("identifying struc_cov")
    # For test:
    #iso_vcab=pd.read_csv("iso_vcab.csv").drop(columns=["Unnamed: 0"])
    hcoorbl=pd.read_csv("./blast_result/h_coordinate_seq_bl_result.csv").drop(columns=["Unnamed: 0"])
    lcoorbl=pd.read_csv("./blast_result/l_coordinate_seq_bl_result.csv").drop(columns=["Unnamed: 0"])

    vcab_iso_cov=add_struc_cov(iso_vcab,hcoorbl,lcoorbl)
    vcab_iso_cov.to_csv("iso_cov_vcab.csv")

    print ("cut mmcif files")
    generate_chain_pdb_total (vcab_iso_cov,"../pdb_struc/full_pdb/","../pdb_struc/chain_pdb/")
    generate_C_pdb_total (vcab_iso_cov,"../pdb_struc/full_pdb/","../pdb_struc/c_pdb/")

    print ("adding pdb_vcb")
    vcab_iso_cov=pd.read_csv("iso_cov_vcab.csv").drop(columns=["Unnamed: 0"])

    vcab_pdb_vcb,err=read_pdb_VC_Boundary_and_C_true_seqs (vcab_iso_cov,"../pdb_struc/c_pdb/")
    if (len(err)>0):
        err.to_csv("./unusual_cases/failed_to_add_pdb_vcb.csv")

    # Filtering out antibodies with abnormal structural coverage
    f_vcab=vcab_pdb_vcb.loc[map(lambda x: x not in ["not classified","check V/C Annotation","unidentified"],vcab_pdb_vcb['Structural Coverage'])]
    print ("adding V C seq")
    ff_vcab=get_V_C_seq(f_vcab)
    ff_vcab=ff_vcab.reset_index(drop=True)
    ff_vcab.to_csv("ff_vcab.csv")

    # For test:
    #ff_vcab=pd.read_csv("ff_vcab.csv").drop(columns=["Unnamed: 0"])

    print ("adding disulfide")
    vcab0=add_disulfide_info(ff_vcab,"../pdb_struc/c_pdb/")
    vcab0.to_csv("disul_vcab.csv")

    print("adding species")
    # For test:
    #vcab0=pd.read_csv("disul_vcab.csv").drop(columns=["Unnamed: 0"])
    vcab1=assign_species_summary(vcab0,flt_vhbl,flt_vlbl,collapsed_v_alleles,c_vh=8,c_vl=8)
    vcab=further_classify_ambiguous(vcab1,flt_vhbl,flt_hbl,flt_vlbl,flt_lbl)
    vcab.to_csv("./vcab.csv")

    # 3. Generate files for the shiny app:
    # 3.1. Generate POPSComp results
    print ("Going through POPSComp analysis...")
    os.chdir("../pops")
    os.system("sh pops.sh ../pdb_struc/c_pdb/") # PDB structures with C region only are inputted for POPSComp analysis
    os.system("sh total_pops.sh ../pdb_struc/chain_pdb/")
    os.chdir("../vcab_db")

    # Change this part!!!!
    print ("generating files for shiny app")
    convert_seq_from_df_to_fasta(vcab,'H_seq',"../seq_db/vcab_db")
    convert_seq_from_df_to_fasta(vcab,'H_coordinate_seq',"../seq_db/vcab_db")

    convert_seq_from_df_to_fasta(vcab,'L_seq',"../seq_db/vcab_db")
    convert_seq_from_df_to_fasta(vcab,'L_coordinate_seq',"../seq_db/vcab_db")

    convert_seq_from_df_to_fasta (vcab,'HV_seq',"../seq_db/vcab_db")
    convert_seq_from_df_to_fasta (vcab,'LV_seq',"../seq_db/vcab_db")
    combine_fastas("../seq_db/vcab_db/fasta","../seq_db/vcab_db")

    os.system ("sh mk_vcab_bl_db.sh")



    os.system("python cal_angles_new.py")
    #os.system("python cal_interface_matrix_new.py")
