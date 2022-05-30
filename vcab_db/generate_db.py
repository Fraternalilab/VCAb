# To generate the VCAb database
# Dongjun Guo, Apr.2022

import pandas as pd
import numpy as np
import json
import requests
from Bio.Blast.Applications import NcbiblastpCommandline as ncbiblp

import Bio
from Bio.PDB import StructureBuilder
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.SeqUtils import seq1
from Bio import SeqIO

import os
import argparse
import time

import warnings
warnings.filterwarnings("ignore")

# Print iterations progress
# This function (printProgressBar()) is from https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters:
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()

######################## 1.1 Filter SAbDab/Get seq info ################################################
def sabdab_filter (df):
    #filter out the pdbs contain not-paired Hchain or Lchain (H chain only or L chain only)
    # df: initial SAbDab table directly downloaded from SAbDab
    H_only=[]
    L_only=[]

    excluded_position=[]

    for i in range(len(df)):
        pdb,Hchain,Lchain=list(df.iloc[i][0:3])

        if type(Hchain) != str:
            L_only.append(list(df.iloc[i]))
            excluded_position.append(i)
        if type(Lchain) != str:
            H_only.append(list(df.iloc[i]))
            excluded_position.append(i)

    HL = df.drop(excluded_position)
    H_df=pd.DataFrame(H_only,columns=df.keys())
    L_df=pd.DataFrame(L_only,columns=df.keys())

    new_HL=HL.sort_values(["pdb", "Hchain"], ascending = (True, True)).reset_index(drop=True)
    return new_HL, H_df, L_df


# download the author-seq directly from PDBe
# The part to incorporate the seq_fragments to generate coordinate_seqs is deleted (in the comment), instead we extract coordinate sequence directly from pdb file
def get_seq_from_pdbe(pdb,c):
    # c: chain id
    URL = 'https://www.ebi.ac.uk/pdbe/api/'
    n_pdb=pdb.lower()
    chainid=c
    mol_data=requests.get(URL + 'pdb/entry/molecules/' + n_pdb)
    #structural_coverage=requests.get(URL+'pdb/entry/polymer_coverage/'+n_pdb+'/chain/'+chainid)

    #if mol_data.status_code == 200 and structural_coverage.status_code == 200:
    if mol_data.status_code == 200:
        # test if API is working
        mol_info=mol_data.json()[n_pdb]
        #coverage_info=structural_coverage.json()[n_pdb]['molecules'][0]['chains'][0]['observed']
        #coverage=[[i['start']['residue_number']-1,i['end']['residue_number']]for i in coverage_info]

        for i in range(len(mol_info)):
            # to find the molecular information(mol_name) of specific chain, because mol_info contains info of multiple entity
            chain=mol_info[i]['in_chains']
            if chainid in chain:
                try:
                    fseq=mol_info[i]['pdb_sequence'] #seq[start-1:end]
                except KeyError as e:
                    fseq="mol: Error - sequence can't be accessed by API"
                    print (f"This sequence can't be downloaded via PDBe API. PDB: {pdb}, ChainID:{chainid}. \n Error message: KeyError:{e}. \n Go to https://www.ebi.ac.uk/pdbe/api/doc/ for details")

                #check if sequence begins with (PCA):
                if '(PCA)' in fseq:
                    fseq=fseq.replace('(PCA)','')
                    #coverage=[[i['start']['residue_number']-1,i['end']['residue_number']-1]for i in coverage_info]
                    # since '(PCA)' is at the start of the sequence, the ending point for the coverage needs to -1

                #sseq=[fseq[i[0]:i[1]] for i in coverage]
                #coor_seq=''.join(sseq)
                ## return the full_seq, coordinate_seq
                #return fseq,coor_seq
                return fseq
        # Use recursion to solve this problem: due to the messy pdb annotation, the chain id could be "HHH", thus the sequence can't be fetched using "H" as the id
        if len(c)>3: # recursion: break condition
            return 'mol: Error - %s' % (mol_data.status_code)
        else:
            new_c=c*3 # due to the messy pdb annotation, some chain is labeled as "H" in sabdab, but actually the chain id is "HHH"
            return get_seq_from_pdbe(pdb,new_c)

    else:
        return 'mol: Error - %s' % (mol_data.status_code)


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

# Add The author_seq(seq) & coordinate_seq  to the sabdab
def add_seq_to_sabdab (o_df):
    # add the following sequence: H_seq, H_coordinate_seq, L_full_seq, L_coordinate_seq
    df=o_df.copy()
    df=df.reset_index(drop=True)
    H_seq_lst=[]
    H_coor_seq_lst=[]
    L_seq_lst=[]
    L_coor_seq_lst=[]

    printProgressBar(0, len(df), prefix = 'Downloading/Extracting sequence:', suffix = 'Complete', length = 50)
    for i in df.index:
        pdb=df.loc[i,'pdb']
        h=df.loc[i,'Hchain'][0]# Take only the first chain name
        l=df.loc[i,'Lchain'][0]
        H_seq=get_seq_from_pdbe(pdb,h)
        L_seq=get_seq_from_pdbe(pdb,l)
        __,__,H_coor_seq,L_coor_seq=read_coor_seq_and_the_start_position_of_chains (pdb,h,l,"../pdb_struc/full_pdb")

        H_seq_lst.append(H_seq)
        H_coor_seq_lst.append(H_coor_seq)
        L_seq_lst.append(L_seq)
        L_coor_seq_lst.append(L_coor_seq)
        printProgressBar(i + 1, len(df), prefix = 'Downloading/Extracting sequence:', suffix = 'Complete', length = 50)

    df['H_seq']=H_seq_lst
    df['H_coordinate_seq']=H_coor_seq_lst
    df['L_seq']=L_seq_lst
    df['L_coordinate_seq']=L_coor_seq_lst

    return df

def collapse_by_coordinate_seqs (df):
    result= [] #collapsed df: nested list

    for i in df.index:
        new=list(df.loc[i])
        n_pdb,n_Hchain,n_Lchain=df.loc[i]['pdb':'Lchain']
        n_htseq=df.loc[i,"H_coordinate_seq"]
        n_ltseq=df.loc[i,"L_coordinate_seq"]
        if i ==0:
            result.append(new)
        else:
            o_pdb,o_Hchain,o_Lchain=result[-1][0:3]
            o_htseq=result[-1][-3]
            o_ltseq=result[-1][-1]
            if n_pdb==o_pdb and n_htseq==o_htseq and n_ltseq==o_ltseq:
                # If coordinate sequences are the same, only attach the chain names to the old(the last) line
                result[-1][1] += n_Hchain # add Hchain name
                result[-1][2] += n_Lchain # add Lchain name
            else:
                # If coordinate sequences are different, attach a seperate new line
                result.append(new)
    result_df=pd.DataFrame(result,columns=df.keys()) # the output would be DataFrame

    # Add the "iden_code" column
    result_df["iden_code"]=[result_df.loc[i,'pdb']+'_'+result_df.loc[i,'Hchain'][0]+result_df.loc[i,'Lchain'][0] for i in result_df.index]

    return result_df

def filter_out_seq_error (df):
    filtered_df=df.copy()
    unusual_index=[] # store the index of unusual cases
    unusual=[]
    for i in df.index:
        row=df.loc[i]
        for item in row:
            if type(item)==str:
                if "mol: Error" in item:
                    unusual.append(row)
                    unusual_index.append(i)
                    break

    filtered_df=filtered_df.drop(filtered_df.index[unusual_index])
    filtered_df=filtered_df.reset_index(drop=True)

    unusual_df=pd.DataFrame(unusual)

    return filtered_df, unusual_df

# Generate seq fasta file for the blast purpose
# IMPORTANT: for update purpose later: might need to modify the file name if we don't want to overwrite the previous file
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

######################## 1.2 Generate vcab (identify isotype, light chain type, structural coverage) except the VCB in the PDB file ########################
chain_type_names={'IGHG1':"IgG1",
                 'IGHG2':"IgG2",
                 'IGHG3':"IgG3",
                 'IGHG4':"IgG4",
                  'IGHA1':"IgA1",
                  'IGHA2':"IgA2",
                  'IGHM':"IgM",
                  'IGHD':"IgD",
                  'IGHE':"IgE",
                'IGLC7':"lambda",
                 'IGKC':"kappa",
                 'IGLC6':"lambda",
                 'IGLC1':"lambda",
                 'IGLC2':"lambda",
                 'IGLC3':"lambda"
                 }

# Collect the domain information
def extract_domain_information_from_imgt_fasta(fn):
    domains={}
    for record in SeqIO.parse(fn,"fasta"):
        title_lst=record.description.split("|")

        allele_name=title_lst[1]
        if_partial=title_lst[13]# if the allele sequence only contains a fragment
        domain_name=title_lst[4]
        domain_seq_len=int(title_lst[11][0:-3])

        if allele_name not in domains.keys():
            d_counter=domain_seq_len
            domains[allele_name]={domain_name:[(1,domain_seq_len),d_counter,if_partial]}
        else:
            last_val=list(domains[allele_name].values())[-1]
            previous_d_counter=last_val[1]
            d_counter=domain_seq_len
            domains[allele_name][domain_name]=[(1+previous_d_counter,domain_seq_len+previous_d_counter),d_counter+previous_d_counter,if_partial]
    return domains

imgt_h_original="../seq_db/ref_db/all_alleles/H_chains/imgt_original.fasta"
imgt_l_original="/Users/dongjung/Downloads/vcab/seq_db/ref_db/all_alleles/L_chains/imgt_original.fasta"

h_domains=extract_domain_information_from_imgt_fasta(imgt_h_original)
l_domains=extract_domain_information_from_imgt_fasta(imgt_l_original)

def generate_bl_result (q_seqs,bl_db,bl_out_name,out_dir):
    # creat the blastp command line
    bl_cmd=ncbiblp(query=q_seqs,db=bl_db,outfmt=10,out=f"{out_dir}/{bl_out_name}_result.csv") # the output file would be the hit table

    # execute the blastp command line & generate the result.csv (the hit table)
    stdout,stderr=bl_cmd()
    header=['qseqid','sseqid','ident','length','mismatches','gapopen','start_ab','end_ab','start_ref','end_ref','e-value','score']
    bl_result=pd.read_csv(f"{out_dir}/{bl_out_name}_result.csv",names=header)

    # Modify the bl_result
    df=bl_result.copy() # no specific reason, just want to type fewer characters/or in case I want to output the initial bl_result
    df["iden_code"]=list(map(lambda x: x.split("-")[0],df["qseqid"]))

    try:
        df["chain_type"]=list(map(lambda x: chain_type_names[x.split("*")[0]],df["sseqid"])) # In the format of "IgG1"
    except:
        df["gene"]=list(map(lambda x: x.split("*")[0],df["sseqid"])) # This is for the blast result for V region
    df["matched_alleles"]=list(df["sseqid"]) # In the format of "IgG1"
    df.to_csv(f"{out_dir}/{bl_out_name}_result.csv")
    return df



def found_special_cases_ab (vhbl,vlbl,hcbl,lcbl,cross_bl,df):
    # To found:
    # 1. the possible domain_exchanged antibody
    # 2. entries with sequence possibly from other species
    # df: vcab, cross_bl: blast heavy against light in order to found possibly dom_exchanged antibody

    best_vhbl=extract_hit_bl_result (vhbl,"ident")
    best_vlbl=extract_hit_bl_result (vlbl,"ident")
    best_hcbl=extract_hit_bl_result (hcbl,"ident")
    best_lcbl=extract_hit_bl_result (lcbl,"ident")

    c_cross_bl=extract_hit_bl_result (cross_bl,"")

    result=[]
    for i in df.index:
        iden_code=df.loc[i,"iden_code"]
        vh_ident=best_vhbl.loc[best_vhbl["iden_code"]==iden_code]
        vl_ident=best_vlbl.loc[best_vlbl["iden_code"]==iden_code]
        hc_ident=best_hcbl.loc[best_hcbl["iden_code"]==iden_code]
        lc_ident=best_lcbl.loc[best_lcbl["iden_code"]==iden_code]

        c_cross_ident=c_cross_bl.loc[c_cross_bl["iden_code"]==iden_code]

        sub_result=""
        if len(vh_ident)!=0 and len(vl_ident)!=0:
            vh_ident_val=vh_ident["ident"].item()
            vl_ident_val=vl_ident["ident"].item()
            if vh_ident_val < 80 or  vl_ident_val < 80:
                sub_result=f"V region sequence identity to human reference is low, indicating sequence possibly from other species (VH:{vh_ident_val}, VL:{vl_ident_val})"
        else:
            sub_result="No significant match found to human reference V sequence, indicating sequence possibly from other species"

        if len(hc_ident)!=0 and len(lc_ident)!=0:
            hc_ident_val=hc_ident["ident"].item()
            lc_ident_val=lc_ident["ident"].item()
            if hc_ident_val < 70 or lc_ident_val < 70:
                sub_result=f"C region sequence identity to human reference is low, indicating sequence possibly from other species (CH:{hc_ident_val}, CL:{lc_ident_val})"
        else:
            sub_result="No significant match found to human reference C sequence, indicating sequence possibly from other species"

        if len(c_cross_ident)!=0:
            if c_cross_ident["ident"].item() >50:
                sub_result="possibly domain_exchanged antibody"

        result.append(sub_result)
    df["special_cases"]=result
    return df

def include_collapsed_alleles(allele,c_alleles_dict):
    all_alleles=[allele]
    all_alleles+=c_alleles_dict[allele]
    return ",".join(all_alleles)

def get_best_hit_for_each_allele(df):
    # get the best hit for each allele in blast result
    result_lst=[]
    for name,sub_df in df.groupby('matched_alleles',sort=False):
        result_lst.append(sub_df.iloc[0,:])
    return pd.DataFrame(result_lst)

def get_chain_type_VCB (iden_code,df,collapsed_alleles):
    # df: the bl result of the full seqs, could be the bl results of H or L chain
    # note: for L chain: the returned type is the LSubtype
    # return 1. chain type of the seq; 2. alternative chain type of the seq (the same score with the best hit); 3.the VC_Boundary of the sequence (author-submitted sequence)
    hit_info=df.loc[df["iden_code"]==iden_code]
    hit_info=hit_info.reset_index(drop=True)

    # Best match:
    # Assign the chain type based on alignment length & identity:
    best_hit_for_each_allele=get_best_hit_for_each_allele(hit_info)
    hit_info=best_hit_for_each_allele.sort_values(by=["ident"],ascending=False)

    best_hit_df=hit_info.iloc[0,]
    best_chain_type=best_hit_df["chain_type"]

    allele=best_hit_df["matched_alleles"]
    ident=best_hit_df["ident"]
    all_alleles=include_collapsed_alleles(allele,collapsed_alleles)
    best_hit_info=f"{best_chain_type}({all_alleles}: Per. Ident: {ident})"

    best_score=best_hit_df["score"]
    VCB=best_hit_df["start_ab"]

    # Alternative match:
    alternative_hit_df=hit_info.loc[hit_info["score"]==best_score]
    alternative_hit_info=""
    #return alternative_hit_df
    for i in alternative_hit_df.index:
        if i == 0:
            # skip the best results
            continue
        elif len(alternative_hit_df) != 0:
            a_chain_type=alternative_hit_df.loc[i,"chain_type"]
            a_allele=alternative_hit_df.loc[i,"matched_alleles"]
            a_all_alleles=include_collapsed_alleles(a_allele,collapsed_alleles)
            a_ident=alternative_hit_df.loc[i,"ident"]
            if a_chain_type != best_chain_type:
                # exclude the cases with the same chain type but different alleles
                alternative_hit_info+=f"{a_chain_type}({a_all_alleles}: Per. Ident: {a_ident}) "

    return best_hit_info,alternative_hit_info,VCB

def aligned_to (start,end,d_name,allele,domains=h_domains):
    # return the aligned coverage for certain domain in specific Ig
    # start & end: start and end positions of the aligned region in the ref_seq in the blast aln result
    # d_name:domain name ('CH1',etc.)
    # allele: allele_name (e.g. IGHG1*01)

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

def get_struc_cov_coor_VCB (iden_code,ig,df,h_domains):
    # Note: df must be the bl result of the true seqs (coordinate seqs)
    # return: 1. align_info; 2.Structural Coverage; 3.VC_Boundary of the coor sequence
    # df: the bl_result of the coordinate sequence
    # iden_code: pdb_HL; ig: the isotype/light chain type of the antibody

    true_hit_info=df.loc[df["iden_code"]==iden_code]
    best_true_hit_info=true_hit_info.iloc[0,]
    hit_allele=best_true_hit_info["matched_alleles"] # the matched allele of the isotype

    # Get the true_VCB info
    true_VCB=best_true_hit_info["start_ab"]

    if ig in ["IgG1","IgG2","IgG3","IgG4","IgA1","IgA2","IgM","IgD","IgE"]:
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

def get_v_c_seqs (seq, vcb, v_or_c):
    # v_or_c can only be "v" or "c"
    # vcb: the boundary between v region and c region
    if v_or_c.lower()=="v":
        return seq[0:vcb]
    else:
        return seq[vcb-1:]

def get_carbohydrate_info_from_pdbe(pdb):
    URL= 'https://www.ebi.ac.uk/pdbe/api/'
    n_pdb=pdb.lower()
    ligand_data=requests.get(URL+'pdb/entry/carbohydrate_polymer/'+n_pdb)
    if ligand_data.status_code==200:
        result={}
        ligand_info=ligand_data.json()[n_pdb]
        for l_dict in ligand_info:
            l_name=l_dict['molecule_name']
            l_chain_id=l_dict['in_chains']
            if l_name not in result.keys():
                result[l_name]=[]
            result[l_name]+=l_chain_id
        result2={k:"".join(v) for k,v in result.items()}
        result_str="; ".join([f"{k} (chain_id:{v})" for k,v in result2.items()])
        return result_str


    else:
        return ""


### FUNCTIONS DESIGNED FOR DOMAIN_SWAPPED_AB: ###
def extract_domain_from_dom_swapped_ab(iden_code,bl,crbl,h_domains,l_domains):
    # bl,crbl should belongs to the same chain sequence
    # bl: straight blast to the corresponding ref profile (e.g. H chain seq to H ref seqs)
    # crbl: cross blast to the other ref profile (e.g. H chain seq to L ref seqs)
    # both bl and crbl should be the collpased blast result. So, the len(bl_info)==len(crbl_info)==1
    bl_info=bl.loc[bl["iden_code"]==iden_code]
    crbl_info=crbl.loc[crbl["iden_code"]==iden_code]

    focus_bl=[i for i in [bl_info,crbl_info] if len(i) !=0]
    focus_bl=[i for i in focus_bl if i["ident"].item()>50] # Only focus on the blast hit with identity larger than 50

    total_aln_info=[]
    displayed_align_info=[]
    total_chain_type_info={} # in the format of{chain_type:ident}, at most one H hit one L hit
    for hit in focus_bl:
        s_ab=hit["start_ab"].item()
        e_ab=hit["end_ab"].item()

        s_ref=hit["start_ref"].item()
        e_ref=hit["end_ref"].item()

        hit_allele=hit["matched_alleles"].item()
        chain_type=hit["chain_type"].item()
        total_chain_type_info[chain_type]=[hit_allele,hit["ident"].item()]
        domains=h_domains
        if chain_type in ["kappa","lambda"]:
            domains=l_domains
        regions=[i for i in domains[hit_allele].keys() if "M" not in i]

        align_info=[(s_ab,e_ab),(s_ref,e_ref)]

        for r in regions:
            subresult=aligned_to (s_ref,e_ref,r,hit_allele,domains)
            if type(subresult)!= str:
                if subresult[1:] !=[0,0]:
                    perc1=round(subresult[1],2)
                    perc2=round(subresult[2],2)
                    domain_name=subresult[0]
                    domain_name+=f"-{chain_type}"
                    sub_string=f"{domain_name}(covers {perc1}% of {subresult[0]}, accounts for {perc2}% of the H_coor_seq)"
                    displayed_align_info.append(sub_string)
                    align_info.append(subresult)

        total_aln_info.append(align_info)

    # pro_d: extract only the domain names
    # vcb: VCBoundary
    pro_d=[]
    total_vcb=[]
    for aln in total_aln_info:
        sub_d=[i[0] for i in aln[2:]]
        pro_d +=sub_d

        vcb=aln[0][0]
        total_vcb.append(vcb)

    h_d=[i for i in pro_d if i != 'C-REGION']
    l_d=[i for i in pro_d if i not in h_d]

    domain_info="; ".join(displayed_align_info)
    final_vcb=total_vcb
    if len(total_vcb)!=0:
        final_vcb=min(total_vcb)
    # take the min value of the aln_start_ab_pos as the VCBoundary, if the seq mapped to both H and l ref

    return focus_bl,total_aln_info,h_d,l_d,domain_info,final_vcb,total_chain_type_info

def annotations_for_dom_swapped_ab (iden_code,hfbl,hf_crbl,lfbl,lf_crbl,htbl,ht_crbl,ltbl,lt_crbl,h_domains,l_domains,c_h_alleles,c_l_alleles):

    # Part1. Identify chain types:
    __,__,__,__,__,hf_vcb,h_ctype=extract_domain_from_dom_swapped_ab(iden_code,htbl,ht_crbl,h_domains,l_domains)
    __,__,__,__,__,lf_vcb,l_ctype=extract_domain_from_dom_swapped_ab(iden_code,ltbl,lt_crbl,h_domains,l_domains)
    htype,ltype,lsubtype,alter_htype,alter_ltype=("","","","","")

    def mergeDictionary(dict_1, dict_2):
        dict_3 = {**dict_1, **dict_2}
        for key, value in dict_3.items():
            if key in dict_1 and key in dict_2:
                # take the one with the highest val (in this case highest identity) as the value
                   dict_3[key] = max([dict_1[key],dict_2[key]])
        return dict_3

    total_ctype=mergeDictionary(h_ctype,l_ctype)

    if len(total_ctype) >1:
        # If len(total_ctype <=1), that means h_ctype or l_ctype is empty, mostly because we force the seq only containing V region to align with C refs
        htype_dict={k:v for k,v in total_ctype.items() if k not in ["kappa","lambda"]}
        ltype_dict={k:v for k,v in total_ctype.items() if k not in htype_dict.keys()}


        ##in case chain_type_dictionary contains multiple chain_types
        ### (at most two, since total_ctype contains 2 H hit and 2 L hit at most)
        #### (that is because h_ctype & l_ctype contains 1 H hit and 1 L hit at most)
        htype_max_ident=max(list(map(lambda x:x[1],list(htype_dict.values()))))
        ltype_max_ident=max(list(map(lambda x:x[1],list(ltype_dict.values()))))
        f_htype_dict={k:v for k,v in htype_dict.items() if v[1]==htype_max_ident}
        f_ltype_dict={k:v for k,v in ltype_dict.items() if v[1]==ltype_max_ident}

        f_htype_tuple=list(f_htype_dict.items())[0]
        f_ltype_tuple=list(f_ltype_dict.items())[0]
        h_allele=f_htype_tuple[1][0]
        l_allele=f_ltype_tuple[1][0]
        all_h_allele=include_collapsed_alleles(h_allele,c_h_alleles)
        all_l_allele=include_collapsed_alleles(l_allele,c_l_alleles)
        htype=f"{f_htype_tuple[0]}({all_h_allele}: Per. Ident: {f_htype_tuple[1][1]})"
        ltype=f"{f_ltype_tuple[0]}({all_l_allele}: Per. Ident: {f_ltype_tuple[1][1]})"
        lsubtype=f_ltype_tuple[1][0].split("*")[0]

        ##alternative chain types:
        alter_htype_dict={k:v for k,v in htype_dict.items() if k!=f_htype_tuple[0]}
        alter_ltype_dict={k:v for k,v in ltype_dict.items() if k!=f_ltype_tuple[0]}
        alter_htype=""
        alter_ltype=""
        if len(alter_htype_dict)!=0:
            alter_h_allele=list(alter_htype_dict.values())[0][0]
            all_alter_h_allele=include_collapsed_alleles(alter_h_allele,c_h_alleles)
            alter_htype=f"{list(alter_htype_dict.keys())[0]}({all_alter_h_allele}: Per. Ident: {list(alter_htype_dict.values())[0][1]})"

    #Part2. Identify the structural coverage and the coor_vcb of the domain_swapped antibody

    h_focus_bl,h_aln_info,h_hd,h_ld,h_domain_info,ht_vcb,ht_ctype=extract_domain_from_dom_swapped_ab(iden_code,htbl,ht_crbl,h_domains,l_domains)
    l_focus_bl,l_aln_info,l_hd,l_ld,l_domain_info,lt_vcb,lt_ctype=extract_domain_from_dom_swapped_ab(iden_code,ltbl,lt_crbl,h_domains,l_domains)


    total_hd=list(set(h_hd+l_hd))
    total_ld=list(set(h_ld+l_ld))

    total_hd.sort()
    total_ld.sort()

    hit_allele_dict=mergeDictionary(ht_ctype,lt_ctype)
    ht_type_dict={k:v for k,v in hit_allele_dict.items() if k not in ["kappa","lambda"]}
    hit_allele=list(ht_type_dict.items())[0][1][0]
    ig_d=[i for i in h_domains[hit_allele].keys() if "M" not in i]
    ig_d.sort()


    hchain,lchain=iden_code.split("_")[1]
    aln_info={}
    aln_info[hchain]=h_aln_info
    aln_info[lchain]=h_aln_info

    domain_info="SPECIAL_CASE:"
    domain_info+=f"Chain {hchain}: {h_domain_info}. "
    domain_info+=f"Chain {lchain}: {l_domain_info}. "

    # Find the structural coverage
    struc_cov=""

    if len(h_focus_bl)==0 or len(l_focus_bl)==0:
        struc_cov='check V/C Annotation'
    else:
        if ht_vcb <=70 or lt_vcb <=70:
            struc_cov='check V/C Annotation'
        else:
            if total_hd[0]!="CH1":
                struc_cov = 'not classified'
                # mostly: structure contains only V region and is forced to align with C region
                # other cases: 1za6 (CH2 deleted)

            elif ig_d==total_hd:
                struc_cov='full antibody'
            elif ("CH3" not in total_hd) and ("CH3-CHS" not in total_hd):
                struc_cov='Fab'
            else:
                struc_cov='not classified'
    return htype,ltype,lsubtype,alter_htype,alter_ltype,hf_vcb,lf_vcb,aln_info,domain_info,struc_cov,ht_vcb,lt_vcb
###------------###

def generate_final_db (o_df,hfbl,lfbl,htbl,ltbl,hf_crbl,lf_crbl,ht_crbl,lt_crbl,h_domains,l_domains,c_h_alleles,c_l_alleles):
    #o_df: the df table from sabdab, labeled with if_domain_swapped
    # hfbl & lfbl: the bl_results of the seqs of H & L chains
    # htbl & ltbl: the bl_results of the coordinate seqs of H & L chains
    # return the final results of vcab (except the pdb_VCB)

    # Note: when fetch the align_info/struc_cov/true_vcb, the bl_result of true seqs are used
    # for some abs, due to the inconsistant of the author_seqs (seq) and the coordinate_seqs (true_seq),
    # the bl_result of true seqs is not available
    # Thus, in these cases, the align_info/struc_cov/true_vcb would be assigned as unidentified
    # And these unidentified abs would be filtered out.

    df=o_df.copy()
    Htype=[]
    Ltype=[]
    LSubtype=[]
    alter_Htype=[]
    alter_Ltype=[]
    sugar_info=[]

    aln_info=[]
    domain_info=[]
    struc_cov=[]

    hvcb_full=[]
    lvcb_full=[]
    hvcb_true=[]
    lvcb_true=[]

    # I include these sequences just for the convenience to build the db
    # for the v seqs to allow the user search in the shiny app:
    hv_full_seq=[]
    lv_full_seq=[]


    for i in df.index:
        iden_code=df.loc[i,"iden_code"]
        pdb_code=df.loc[i,"pdb"]
        carbohydrate=get_carbohydrate_info_from_pdbe(pdb_code)
        if_special_case=df.loc[i,"special_cases"]

        this_h_info,this_l_info,this_l_subtype,this_alter_h_info,this_alter_l_info,this_hf_vcb,this_lf_vcb,this_aln_info,this_domain_info,this_struc_cov,this_ht_vcb,this_lt_vcb=[""]*12
        if if_special_case=="True" or if_special_case==True:
            this_h_info,this_l_info,this_l_subtype,this_alter_h_info,this_alter_l_info,this_hf_vcb,this_lf_vcb,this_aln_info,this_domain_info,this_struc_cov,this_ht_vcb,this_lt_vcb=annotations_for_dom_swapped_ab (iden_code,hfbl,hf_crbl,lfbl,lf_crbl,htbl,ht_crbl,ltbl,lt_crbl,h_domains,l_domains,c_h_alleles,c_l_alleles)
        else:

            try:
                this_h_info,this_alter_h_info,this_hf_vcb=get_chain_type_VCB (iden_code,hfbl,c_h_alleles)
                this_l_info,this_alter_l_info,this_lf_vcb=get_chain_type_VCB (iden_code,lfbl,c_l_alleles)
            except:
                this_h_info,this_alter_h_info,this_hf_vcb=("","","")
                this_l_info,this_alter_l_info,this_lf_vcb=("","","")
                print (iden_code)

            this_h_type=this_h_info.split("(")[0]
            this_l_type=this_l_info.split("(")[0]
            try:
                this_l_subtype=this_l_info.split("(")[1].split("*")[0]
            except:
                print(f"error:line504:{iden_code} \n {this_l_info} \n {this_l_type}")
                this_l_subtype=""

            try:
                this_aln_info,this_domain_info,this_struc_cov,this_ht_vcb=get_struc_cov_coor_VCB (iden_code,this_h_type,htbl,h_domains)
                __,__,__,this_lt_vcb=get_struc_cov_coor_VCB (iden_code,this_l_type,ltbl,h_domains)
            except:
                this_aln_info,this_struc_cov,this_ht_vcb,this_lt_vcb,this_domain_info=("unidentified","unidentified","unidentified","unidentified","unidentified")


        Htype.append(this_h_info)
        Ltype.append(this_l_info)
        LSubtype.append(this_l_subtype)

        alter_Htype.append(this_alter_h_info)
        alter_Ltype.append(this_alter_l_info)
        sugar_info.append(carbohydrate)

        hvcb_full.append(this_hf_vcb)
        lvcb_full.append(this_lf_vcb)

        aln_info.append(this_aln_info)
        domain_info.append(this_domain_info)
        struc_cov.append(this_struc_cov)
        hvcb_true.append(this_ht_vcb)
        lvcb_true.append(this_lt_vcb)

    df["Htype"]=Htype
    df["Ltype"]=Ltype
    df["LSubtype"]=LSubtype
    df["Alternative_Htype"]=alter_Htype
    df["Alternative_Ltype"]=alter_Ltype
    df["carbohydrate_polymer"]=sugar_info

    df["align_info"]=aln_info
    df["Domains in HC"]=domain_info
    df["Structural Coverage"]=struc_cov

    df["H_seq_VC_boundary"]=hvcb_full
    df["L_seq_VC_boundary"]=lvcb_full
    df["H_coordinate_seq_VC_boundary"]=hvcb_true
    df["L_coordinate_seq_VC_boundary"]=lvcb_true

    return df

# Filter out the unusual cases (mostly no C region in the true seq)
def find_unusual_cases (df):
    filtered_df=df.copy()
    unusual_index=[] # store the index of unusual cases
    unusual=[]
    for i in df.index:
        row=df.loc[i]
        for item in row:
            if item == "unidentified type" or item == "unidentified":
                unusual.append(row)
                unusual_index.append(i)
                break

    filtered_df=filtered_df.drop(filtered_df.index[unusual_index])
    filtered_df=filtered_df.reset_index(drop=True)

    unusual_df=pd.DataFrame(unusual)

    return filtered_df, unusual_df


######################### 2.1 CUT THE PDB FILES  ######################################################
def keep_chain (model,chains):
    chain_to_remove = []
    for c in model:
        if c.id not in chains:
            chain_to_remove.append(c.id)
    for c in chain_to_remove:
        model.detach_child(c)
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
    return chain

def generate_C_pdb (pdb_c,h_vc,l_vc,in_dir,out_dir):
    pdb=pdb_c[0:4]
    h=pdb_c[5]
    l=pdb_c[6]



    parser = PDBParser()
    structure=parser.get_structure(pdb, f'{in_dir}/{pdb}.pdb')
    model=structure[0]

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
        try:
            pdb_c=df.loc[i,'iden_code']
            h_vc=df.loc[i,'H_coordinate_seq_VC_boundary']
            l_vc=df.loc[i,'L_coordinate_seq_VC_boundary']

            if os.path.exists(f"{out_dir}/{pdb_c}_C.pdb"):
                continue
            generate_C_pdb (pdb_c,h_vc,l_vc,in_dir,out_dir)

        except:
            error.append(df.loc[i,:])
    return pd.DataFrame(error,columns=df.keys())

def generate_chain_pdb_total (df,in_dir,out_dir):
    error=[]
    for i in df.index:
        try:
            pdb_c=df.loc[i,'iden_code']

            pdb=pdb_c[0:4]
            h=pdb_c[5]
            l=pdb_c[6]

            if os.path.exists(f"{out_dir}/{pdb_c}.pdb"):
                print (pdb_c)
                continue

            parser = PDBParser()
            structure=parser.get_structure(pdb, f'{in_dir}/{pdb}.pdb')
            model=structure[0]

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


######################### 2.2 ADD PDB_VC_BOUNDARY (In order to specify the VC_Boundary in the 3D viewer of the webserver)  ######################################################
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

######################### 2.3 Separate AUTHOR SEQS into V and C REGIONS  ######################################################
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

def generate_pdbid_list (df,out_dir):
    pdbid=list(set(df['pdb']))
    with open(f'{out_dir}/pdbid_lst.txt','w') as output:
        output.write(','.join(pdbid))

########################## 3.3 Calculate the disulfide bond ##########################
def calc_dis_between_cys (iden_code,c_pdb_dir):
    # cat_chain: concatenated chain object of the H and L chain
    parser=Bio.PDB.PDBParser()
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

def extract_hit_bl_result (bl,standard=""):
    # extract the best (the first one by the default order given by blast) hit in bl_result
    # bl: the blast result dataframe to be extracted
    # Note: the difference between this function and "extract_hit_bl_result_for_shiny":
    ## in this function, the chain type would not be needed. This is mainly for the domain extraction step in domain-swaaped antibodies.
    ## The reason why "extract_hit_bl_result_for_shiny" need the chain type is because when generating the sequence coverage plot in shiny app, we want to compare the author_seq and coor_seq with the same allele
    ## Because some coor_seq lacking some fragment compared with author_seq, this minor difference might lead to different hit alleles when we blast them against the ref_alleles

    # standard can only be "default" or "ident"
    c_bl_lst=[]
    iden_code=list(set(bl["iden_code"]))
    for i in iden_code:
        # extract the first hit as the best hit:
        sub=pd.DataFrame(bl.loc[bl["iden_code"]==i].iloc[0,]).T
        # extract the one with highest identity as the best hit:
        if standard=="ident":
            sub=bl.loc[bl["iden_code"]==i].sort_values(by=["ident"], ascending=False)
            sub=pd.DataFrame(sub.iloc[0,]).T
        c_bl_lst.append(sub)
    c_bl=pd.concat(c_bl_lst)

    return c_bl

def extract_hit_bl_result_for_shiny (bl_df,horltype,new_bl_name,df):
    # extract only the bl_result of the hit chain type (when the hit chain type is known)
    # bl_df: the blast result dataframe to be extracted
    # horltype can only be "Htype" or "Ltype"
    new_bl_lst=[]
    for i in df.index:
        iden_code=df.loc[i,"iden_code"]
        ctype=df.loc[i,horltype]
        allele_info=ctype.split("(")[1]
        allele=allele_info.split(",")[0]

        sub_bl_df=pd.DataFrame(bl_df.loc[(bl_df["iden_code"]==iden_code)&(bl_df["matched_alleles"]==allele)].iloc[0,]).T

        new_bl_lst.append(sub_bl_df)
    new_bl=pd.concat(new_bl_lst)

    new_bl.to_csv(f"{new_bl_name}.csv")

########################## Apply the functions #######################################
# The directory to run this is the directory of this python file
parser = argparse.ArgumentParser(description="Generate the VCAb database")
parser.add_argument("--in_file",help="the input tsv file containing paired H and L chains from antibody PDBs, such as the csv file downloaded from SAbDab")
args = parser.parse_args()

# 1.1
sabdab=pd.read_csv(args.in_file,sep='\t')
#sub_sabdab=sabdab.iloc[0:50,:] # For test

paired,H_df,L_df=sabdab_filter(sabdab)

# download the pdb files
print ("Downloading PDB structures...")
generate_pdbid_list(paired,"../pdb_struc/") # out_dir: ../pdb_struc
os.system("sh ../pdb_struc/pdb_download.sh -f ../pdb_struc/pdbid_lst.txt -o ../pdb_struc/full_pdb/ -p")

paired_seq=add_seq_to_sabdab (paired) # add sequences to the Abs with paired chains
# For testing: show the acquired df now
#paired_seq.to_csv("paired.csv")
filtered_paired_seq,seq_err=filter_out_seq_error (paired_seq) # filter out abs with seqs can't downloaded from PDBe (ususally because of the wrong annotation of chain)
## For testing: show the acquired df now
#filtered_paired_seq.to_csv("paired.csv")
seq_err=seq_err.reset_index(drop=True)
seq_err.to_csv("./unusual_cases/seq_err.csv")

print ("Generating VCAb database ...")
cHL=collapse_by_coordinate_seqs (filtered_paired_seq)
# For testing: show the acquired df now
#cHL.to_csv("cHL_paired.csv")
#cHL=pd.read_csv("cHL_paired.csv").drop(columns="Unnamed: 0")


# generate the fasta files
convert_seq_from_df_to_fasta(cHL,'H_seq',"../seq_db/vcab_db")
convert_seq_from_df_to_fasta(cHL,'H_coordinate_seq',"../seq_db/vcab_db")

convert_seq_from_df_to_fasta(cHL,'L_seq',"../seq_db/vcab_db")
convert_seq_from_df_to_fasta(cHL,'L_coordinate_seq',"../seq_db/vcab_db")

# 1.2 Identify isotypes and struc_coverage
hfqseqs="../seq_db/vcab_db/H_seq.fasta"
lfqseqs="../seq_db/vcab_db/L_seq.fasta"

htqseqs="../seq_db/vcab_db/H_coordinate_seq.fasta"
ltqseqs="../seq_db/vcab_db/L_coordinate_seq.fasta"

h_ref_db="../seq_db/ref_db/all_alleles/H_chains/unique_alleles.fasta"
l_ref_db="../seq_db/ref_db/all_alleles/L_chains/unique_light_alleles.fasta"

vh_ref_db="../seq_db/ref_db/all_alleles/v_gene/heavy/unique_ighv.fasta"
vl_ref_db="../seq_db/ref_db/all_alleles/v_gene/light/unique_vl.fasta"

# Do the blast:
hbl=generate_bl_result (hfqseqs,h_ref_db,"h_seq_bl","./blast_result")
lbl=generate_bl_result (lfqseqs,l_ref_db,"l_seq_bl","./blast_result")

vhbl=generate_bl_result (hfqseqs,vh_ref_db,"vh_seq_bl","./blast_result")
vlbl=generate_bl_result (lfqseqs,vl_ref_db,"vl_seq_bl","./blast_result")

hcoorbl=generate_bl_result (htqseqs,h_ref_db,"h_coordinate_seq_bl","./blast_result")
lcoorbl=generate_bl_result (ltqseqs,l_ref_db,"l_coordinate_seq_bl","./blast_result")

# Do the cross-blast in order to find possible domain-swapped antibodies
l_coor_cross_bl=generate_bl_result (ltqseqs,h_ref_db,"cross_bl_l_coor_seq","./blast_result")
h_coor_cross_bl=generate_bl_result (htqseqs,l_ref_db,"cross_bl_h_coor_seq","./blast_result")

l_author_cross_bl=generate_bl_result (lfqseqs,h_ref_db,"cross_bl_l_author_seq","./blast_result")
h_author_cross_bl=generate_bl_result (hfqseqs,l_ref_db,"cross_bl_h_author_seq","./blast_result")

#cHL_dom_swap_labeled=found_dom_exchanged_ab(l_coor_cross_bl,cHL)
cHL_dom_swap_labeled=found_special_cases_ab(vhbl,vlbl,hbl,lbl,l_coor_cross_bl,cHL)
# For test:
cHL_dom_swap_labeled.to_csv("special_cases_labeled.csv")

"""# For test:
hbl=pd.read_csv("./blast_result/h_seq_bl_result.csv").drop(columns=["Unnamed: 0"])
lbl=pd.read_csv("./blast_result/l_seq_bl_result.csv").drop(columns=["Unnamed: 0"])
hcoorbl=pd.read_csv("./blast_result/h_coordinate_seq_bl_result.csv").drop(columns=["Unnamed: 0"])
lcoorbl=pd.read_csv("./blast_result/l_coordinate_seq_bl_result.csv").drop(columns=["Unnamed: 0"])

l_coor_cross_bl=pd.read_csv("./blast_result/cross_bl_l_coor_seq_result.csv").drop(columns=["Unnamed: 0"])
h_coor_cross_bl=pd.read_csv("./blast_result/cross_bl_h_coor_seq_result.csv").drop(columns=["Unnamed: 0"])
l_author_cross_bl=pd.read_csv("./blast_result/cross_bl_l_author_seq_result.csv").drop(columns=["Unnamed: 0"])
h_author_cross_bl=pd.read_csv("./blast_result/cross_bl_h_author_seq_result.csv").drop(columns=["Unnamed: 0"])"""

# Get the collapsed bl_result:
htbl=extract_hit_bl_result(hcoorbl)
ltbl=extract_hit_bl_result(lcoorbl)

hfbl=extract_hit_bl_result(hbl)
lfbl=extract_hit_bl_result(lbl)

ht_crbl=extract_hit_bl_result (h_coor_cross_bl)
lt_crbl=extract_hit_bl_result (l_coor_cross_bl)

hf_crbl=extract_hit_bl_result (h_author_cross_bl)
lf_crbl=extract_hit_bl_result (l_author_cross_bl)

"""# For test:
cHL_dom_swap_labeled=pd.read_csv("cHL_dom_swap_labeled.csv").drop(columns=["Unnamed: 0"])
cHL_dom_swap_labeled=cHL_dom_swap_labeled.rename(columns={"domain_swapped_ab":"special_cases"})"""
collapsed_h_alleles=json.load(open("../seq_db/ref_db/all_alleles/H_chains/collapsed_h_alleles.json","r"))
collapsed_l_alleles=json.load(open("../seq_db/ref_db/all_alleles/L_chains/collapsed_l_alleles.json","r"))
total_vcab=generate_final_db (cHL_dom_swap_labeled,hfbl,lfbl,htbl,ltbl,hf_crbl,lf_crbl,ht_crbl,lt_crbl,h_domains,l_domains,collapsed_h_alleles,collapsed_l_alleles)
# For test:
total_vcab.to_csv("test_vcab.csv")

# Extract unusual cases:
vcab,unusual=find_unusual_cases (total_vcab)
f_vcab=vcab.loc[(vcab["Structural Coverage"]!="check V/C Annotation") & (vcab["Structural Coverage"]!="not classified")] # exclude unusual entries
sus_v_c=vcab.loc[vcab["Structural Coverage"]=="check V/C Annotation"] # These entries are probably Abs only containing V region.
struc_cov_unclassified=vcab.loc[vcab["Structural Coverage"]=="not classified"]
sus_v_c=sus_v_c.reset_index(drop=True)
unusual=unusual.reset_index(drop=True)
struc_cov_unclassified=struc_cov_unclassified.reset_index(drop=True)
sus_v_c.to_csv("./unusual_cases/suspicious_v_c_annotation_in_SAbDab.csv")
unusual.to_csv("./unusual_cases/unusual.csv")
struc_cov_unclassified.to_csv("./unusual_cases/struc_cov_not_classified.csv")


# 2.1.CUT THE PDB FILES (also add the pdb_VCB to vcab)
generate_chain_pdb_total (f_vcab,"../pdb_struc/full_pdb/","../pdb_struc/chain_pdb/")
generate_C_pdb_total (f_vcab,"../pdb_struc/full_pdb/","../pdb_struc/c_pdb/")

# 2.2 ADD PDB_VC_BOUNDARY & true seqs of C region TO VCAB
vcab_pdb_vcb,err=read_pdb_VC_Boundary_and_C_true_seqs (f_vcab,"../pdb_struc/c_pdb/")
if (len(err)>0):
    err.to_csv("./unusual_cases/failed_to_add_pdb_vcb.csv")
# 2.3 ADD full seqs of V and C
ff_vcab=get_V_C_seq(vcab_pdb_vcb)
ff_vcab=ff_vcab.reset_index(drop=True)
ff_vcab=add_disulfide_info(ff_vcab,"../pdb_struc/c_pdb/")
ff_vcab.to_csv("new_vcab.csv")

# Collapse bl_result to make the file smaller and the shiny app to load the file faster:
extract_hit_bl_result_for_shiny (hbl,"Htype","./blast_result/best_h_seq_bl_result",ff_vcab)
extract_hit_bl_result_for_shiny (lbl,"Ltype","./blast_result/best_l_seq_bl_result",ff_vcab)

extract_hit_bl_result_for_shiny (hcoorbl,"Htype","./blast_result/best_h_coordinate_seq_bl_result",ff_vcab)
extract_hit_bl_result_for_shiny (lcoorbl,"Ltype","./blast_result/best_l_coordinate_seq_bl_result",ff_vcab)

# 3. Generate files for the shiny app:
# 3.1. Generate POPSComp results
print ("Going through POPSComp analysis...")
os.chdir("../pops")
os.system("sh pops.sh ../pdb_struc/c_pdb/") # PDB structures with C region only are inputted for POPSComp analysis
os.chdir("../vcab_db")

# 3.2. Generate BLAST databases
# generate fasta files
print ("Generating BLAST databases...")
#ff_vcab=pd.read_csv("new_vcab.csv").drop(columns="Unnamed: 0")
convert_seq_from_df_to_fasta(ff_vcab,'HV_seq',"../seq_db/vcab_db")
convert_seq_from_df_to_fasta(ff_vcab,'LV_seq',"../seq_db/vcab_db")

#Important: rewrite the H_seq.fasta and L_seq.fasta.
# Previous files containing sequences in cHL are used to identify isotype & structural coverage.
# Now we need to build the database only containing sequences from ff_vcab. Some entries are excluded from cHL, because they are find to only include V region, etc.
convert_seq_from_df_to_fasta(ff_vcab,'H_seq',"../seq_db/vcab_db")
convert_seq_from_df_to_fasta(ff_vcab,'L_seq',"../seq_db/vcab_db")
os.system("cat ../seq_db/vcab_db/HV_seq.fasta ../seq_db/vcab_db/LV_seq.fasta > ../seq_db/vcab_db/all_v_seq.fasta") # Combining two fasta files together to generate the fasta file containing all sequences of V region
os.system("cat ../seq_db/vcab_db/H_seq.fasta ../seq_db/vcab_db/L_seq.fasta > ../seq_db/vcab_db/all_full_seq.fasta") # Combining two fasta files together to generate the fasta file containing all sequences of V + C region (full sequence)
# make BLAST database. Note: BLAST should be installed on command line
os.system("makeblastdb -in ../seq_db/vcab_db/H_seq.fasta -dbtype prot")
os.system("makeblastdb -in ../seq_db/vcab_db/L_seq.fasta -dbtype prot")
os.system("makeblastdb -in ../seq_db/vcab_db/HV_seq.fasta -dbtype prot")
os.system("makeblastdb -in ../seq_db/vcab_db/LV_seq.fasta -dbtype prot")
os.system("makeblastdb -in ../seq_db/vcab_db/all_v_seq.fasta -dbtype prot")
os.system("makeblastdb -in ../seq_db/vcab_db/all_full_seq.fasta -dbtype prot")

# Run cal_angles.py & cal_interface_matrix.py
os.system("python cal_angles.py")
os.system("python cal_interface_matrix.py")
