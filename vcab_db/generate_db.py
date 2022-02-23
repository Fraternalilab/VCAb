import pandas as pd
import numpy as np
import requests
from Bio.Blast.Applications import NcbiblastpCommandline as ncbiblp

from Bio.PDB import StructureBuilder
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.SeqUtils import seq1

import os
import argparse
import time

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

######################## Part 1 Filter SAbDab/Get seq info ################################################
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

    return HL, H_df, L_df

# download the author-seq directly from PDBe & incorporate the seq_fragments to generate coordinate_seqs
def get_seq_from_pdbe(pdb,c):
    # c: chain id
    URL = 'https://www.ebi.ac.uk/pdbe/api/'
    n_pdb=pdb.lower()
    chainid=c
    mol_data=requests.get(URL + 'pdb/entry/molecules/' + n_pdb)
    structural_coverage=requests.get(URL+'pdb/entry/polymer_coverage/'+n_pdb+'/chain/'+chainid)

    if mol_data.status_code == 200 and structural_coverage.status_code == 200:
        # test if API is working
        mol_info=mol_data.json()[n_pdb]
        coverage_info=structural_coverage.json()[n_pdb]['molecules'][0]['chains'][0]['observed']
        coverage=[[i['start']['residue_number']-1,i['end']['residue_number']]for i in coverage_info]

        for i in range(len(mol_info)):
            # to find the molecular information(mol_name) of specific chain, because mol_info contains info of multiple entity
            chain=mol_info[i]['in_chains']
            if chainid in chain:
                fseq=mol_info[i]['pdb_sequence'] #seq[start-1:end]

                #check if sequence begins with (PCA):
                if '(PCA)' in fseq:
                    fseq=fseq.replace('(PCA)','')
                    coverage=[[i['start']['residue_number']-1,i['end']['residue_number']-1]for i in coverage_info]
                    # since '(PCA)' is at the start of the sequence, the ending point for the coverage needs to -1

                sseq=[fseq[i[0]:i[1]] for i in coverage]
                coor_seq=''.join(sseq)
                # return the full_seq, coordinate_seq
                return fseq,coor_seq


    elif mol_data.status_code == 200 and structural_coverage.status_code == 404:
        if len(c)>3: # recursion: break condition
            return 'mol: Error - %s, s_coverage: Error - %s' % (mol_data.status_code, structural_coverage.status_code),'mol: Error - %s, s_coverage: Error - %s' % (mol_data.status_code, structural_coverage.status_code)
        else:
            new_c=c*3 # due to the messy pdb annotation, some chain is labeled as "H" in sabdab, but actually the chain id is "HHH"
            return get_seq_from_pdbe(pdb,new_c)
    else:
        return 'mol: Error - %s, s_coverage: Error - %s' % (mol_data.status_code, structural_coverage.status_code),'mol: Error - %s, s_coverage: Error - %s' % (mol_data.status_code, structural_coverage.status_code)

# Add The author_seq(seq) & coordinate_seq  to the sabdab
def add_seq_to_sabdab (o_df):
    # add the following sequence: H_seq, H_coordinate_seq, L_full_seq, L_coordinate_seq
    df=o_df.copy()
    df=df.reset_index(drop=True)
    H_seq_lst=[]
    H_coor_seq_lst=[]
    L_seq_lst=[]
    L_coor_seq_lst=[]

    printProgressBar(0, len(df), prefix = 'Downloading sequence:', suffix = 'Complete', length = 50)
    for i in df.index:
        pdb=df.loc[i,'pdb']
        h=df.loc[i,'Hchain'][0]# Take only the first chain name
        l=df.loc[i,'Lchain'][0]
        H_seq,H_coor_seq=get_seq_from_pdbe(pdb,h)
        L_seq,L_coor_seq=get_seq_from_pdbe(pdb,l)

        H_seq_lst.append(H_seq)
        H_coor_seq_lst.append(H_coor_seq)
        L_seq_lst.append(L_seq)
        L_coor_seq_lst.append(L_coor_seq)
        printProgressBar(i + 1, len(df), prefix = 'Downloading sequence:', suffix = 'Complete', length = 50)

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
        file.write(full_name+'\n')
        file.write(seq+'\n')
    file.close()

######################## Part 2 Generate vcab except the VCB of the PDB file ########################
chain_type_names={'sp|P01857|IGHG1_HUMAN':"IgG1",
                 'sp|P01859|IGHG2_HUMAN':"IgG2",
                 'sp|P01860|IGHG3_HUMAN':"IgG3",
                 'sp|P01861|IGHG4_HUMAN':"IgG4",
                  'sp|P01876|IGHA1_HUMAN':"IgA1",
                  'sp|P01877|IGHA2_HUMAN':"IgA2",
                  'sp|P01871|IGHM_HUMAN':"IgM",
                  'sp|P01880|IGHD_HUMAN':"IgD",
                  'sp|P01854|IGHE_HUMAN':"IgE",
                'sp|A0M8Q6|IGLC7_HUMAN':"IgLC7",
                 'sp|P01834|IGKC_HUMAN':"IgKC",
                 'sp|P0CF74|IGLC6_HUMAN':"IgLC6",
                 'sp|P0CG04|IGLC1_HUMAN':"IgLC1",
                 'sp|P0DOY2|IGLC2_HUMAN':"IgLC2",
                 'sp|P0DOY3|IGLC3_HUMAN':"IgLC3"
                 }

# Structural Coverage Definition (from uniprot)
CHs_in_igs={}
CHs_in_igs['IgG1']=[(1,98),(99,110),(111,223),(224,330),'']
CHs_in_igs['IgG2']=[(1,98),(99,110),(111,219),(220,326),'']
CHs_in_igs['IgG3']=[(1,98),(99,160),(161,270),(271,376),'']
CHs_in_igs['IgG4']=[(1,98),(99,110),(111,220),(221,327),'']
CHs_in_igs['IgA1']=[(6,98),'',(125,220),(228,330),'']
CHs_in_igs['IgA2']=[(6,98),'',(112,207),(215,317),'']
CHs_in_igs['IgM']=[(1,105),'',(106,217),(218,323),(324,452)]
CHs_in_igs['IgD']=[(6,98),'',(175,263),(267,373),'']
CHs_in_igs['IgE']=[(6,103),'',(112,210),(214,318),(324,423)]
domains=pd.DataFrame(CHs_in_igs,index=['CH1','hinge','CH2','CH3','CH4'])

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
    df["chain_type"]=list(map(lambda x: chain_type_names[x],df["sseqid"])) # In the format of "IgG1"
    return df

def get_chain_type_VCB (iden_code,df):
    # df: the bl result of the full seqs, could be the bl results of H or L chain
    # note: for L chain: the returned type is the LSubtype
    # return 1. chain type of the seq; 2.the VC_Boundary of the sequence (author-submitted sequence)
    hit_info=df.loc[df["iden_code"]==iden_code]
    best_hit_info=hit_info.iloc[0,]

    chain_type=best_hit_info["chain_type"]
    VCB=best_hit_info["start_ab"]
    return chain_type,VCB

def aligned_to (start,end,d_name,ig):
    # return the aligned coverage for certain domain in specific Ig
    # start and end point of the aligned region in the standard isotype (uniprot)
    # d_name:domain name ('CH1',etc.)
    # ig:ig_name, 'igg1',etc

    d_region=domains[ig][d_name] #domain region
    if d_region == '':
        return d_name + ': not defined region'
    else:
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

def get_struc_cov_coor_VCB (iden_code,ig,df):
    # Note: df must be the bl result of the true seqs (coordinate seqs)
    # return: 1. align_info; 2.Structural Coverage; 3.VC_Boundary of the true sequence
    # iden_code: pdb_HL; ig: the isotype of the antibody

    true_hit_info=df.loc[df["iden_code"]==iden_code]
    best_true_hit_info=true_hit_info.iloc[0,]

    # Get the true_VCB info
    true_VCB=best_true_hit_info["start_ab"]

    if ig in ["IgG1","IgG2","IgG3","IgG4","IgA1","IgA2","IgM","IgD","IgE"]:
        # if the chain is heavy chain
        # Get the align_info
        s,e=best_true_hit_info["start_ref"],best_true_hit_info["end_ref"]
        s_t,e_t=best_true_hit_info["start_ab"],best_true_hit_info["end_ab"]#start and end point of the aligned region in the target antibody

        regions=list(domains.index)
        align_info=[(s_t,e_t),(s,e)]

        for r in regions:
            subresult=aligned_to (s,e,r,ig)
            if type(subresult)== list:
                align_info.append(subresult)

        # Get the struc_cov
        ig_d=[d for d in domains.index if domains [ig][d]!=''] # the domain list for the whole antibody
        pro_d=[a[0] for a in align_info[2:]] # the domain list for this antibody
        domain_info=",".join(pro_d)
        if true_VCB > 70:
            if ig_d==pro_d:
                struc_cov = 'full antibody'
            elif pro_d==['CH1'] or pro_d==['CH1','hinge']:
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

def generate_final_db (o_df,hfbl,lfbl,htbl,ltbl):
    #o_df: the df table from sabdab
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

        try:
            this_h_type,this_hf_vcb=get_chain_type_VCB (iden_code,hfbl)
            this_l_subtype,this_lf_vcb=get_chain_type_VCB (iden_code,lfbl)
        except:
            this_h_subtype,this_hf_vcb=("","")
            this_l_subtype,this_lf_vcb=("","")
            print (iden_code)

        # get the Ltype:
        if this_l_subtype=="IgKC":
            this_l_type="kappa"
        else:
            this_l_type="lambda"

        try:
            this_aln_info,this_domain_info,this_struc_cov,this_ht_vcb=get_struc_cov_coor_VCB (iden_code,this_h_type,htbl)
            __,__,__,this_lt_vcb=get_struc_cov_coor_VCB (iden_code,this_l_type,ltbl)
        except:
            this_aln_info,this_struc_cov,this_ht_vcb,this_lt_vcb,this_domain_info=("unidentified","unidentified","unidentified","unidentified","unidentified")

        Htype.append(this_h_type)
        Ltype.append(this_l_type)
        LSubtype.append(this_l_subtype)
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


######################### PART 2.1 CUT THE PDB FILES  ######################################################
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
        if (c+1<v_c) or (r.id[0]!=' '):
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
            h_vc=df.loc[i,'H_true_seq_VC_boundary']
            l_vc=df.loc[i,'L_true_seq_VC_boundary']

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


######################### PART 2.2 ADD PDB_VC_BOUNDARY  ######################################################
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

        try:
            structure=parser.get_structure(iden_code, f'{pdb_dir}{iden_code}_C.pdb')
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
            chains = {chain.id.upper():seq1(''.join(residue.resname for residue in chain)) for chain in structure.get_chains()}

            h_seq = chains[h]
            l_seq = chains[l]
        except:
            # file not found
            hchain_start = "pdb not found"
            lchain_start = "pdb not found"
            h_seq = ""
            l_seq = ""

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

######################### PART 2.3 ADD AUTHOR SEQS of V and C REGIONS  ######################################################
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

########################## Apply the functions #######################################
# The directory to run this is the directory of this python file
parser = argparse.ArgumentParser(description="Generate the VCAb database")
parser.add_argument("--in_file",help="the input tsv file containing paired H and L chains from antibody PDBs, such as the csv file downloaded from SAbDab")
args = parser.parse_args()

# 1.1
sabdab=pd.read_csv(args.in_file,sep='\t')

paired,H_df,L_df=sabdab_filter(sabdab)

paired_seq=add_seq_to_sabdab (paired) # add sequences to the Abs with paired chains
## For testing: show the acquired df now
#paired_seq.to_csv("paired.csv")
filtered_paired_seq,seq_err=filter_out_seq_error (paired_seq) # filter out abs with seqs can't downloaded from PDBe (ususally because of the wrong annotation of chain)
## For testing: show the acquired df now
#filtered_paired_seq.to_csv("paired.csv")
seq_err.to_csv("./unusual_cases/seq_err.csv")

cHL=collapse_by_coordinate_seqs (filtered_paired_seq)
## For testing: show the acquired df now
#cHL.to_csv("cHL_paired.csv")

# download the pdb files
generate_pdbid_list(cHL,"../pdb_struc/") # out_dir: ../pdb_struc
os.system("sh ../pdb_struc/pdb_download.sh -f ../pdb_struc/pdbid_lst.txt -o ../pdb_struc/full_pdb/ -p")

# generate the fasta files
convert_seq_from_df_to_fasta(cHL,'H_seq',"../seq_db/vcab_db")
convert_seq_from_df_to_fasta(cHL,'H_coordinate_seq',"../seq_db/vcab_db")

convert_seq_from_df_to_fasta(cHL,'L_seq',"../seq_db/vcab_db")
convert_seq_from_df_to_fasta(cHL,'L_coordinate_seq',"../seq_db/vcab_db")

# 1.2
hfqseqs="../seq_db/vcab_db/H_seq.fasta"
lfqseqs="../seq_db/vcab_db/L_seq.fasta"

htqseqs="../seq_db/vcab_db/H_coordinate_seq.fasta"
ltqseqs="../seq_db/vcab_db/L_coordinate_seq.fasta"

h_ref_db="../seq_db/ref_db/human_IGH_db/human_IGH.fasta"
l_ref_db="../seq_db/ref_db/human_light_chain_db/human_light_constant.fasta"

hbl=generate_bl_result (hfqseqs,h_ref_db,"h_seq_bl","./blast_result")
lbl=generate_bl_result (lfqseqs,l_ref_db,"l_seq_bl","./blast_result")

hcoorbl=generate_bl_result (htqseqs,h_ref_db,"h_coordinate_seq_bl","./blast_result")
lcoorbl=generate_bl_result (ltqseqs,l_ref_db,"l_coordinate_seq_bl","./blast_result")

total_vcab=generate_final_db (cHL,hbl,lbl,hcoorbl,lcoorbl)

# Extracted unusual cases:
vcab,unusual=find_unusual_cases (total_vcab)
f_vcab=vcab.loc[(vcab["Structural Coverage"]!="check V/C Annotation") & (vcab["Structural Coverage"]!="not classified")] # exclude unusual entries
sus_v_c=vcab.loc[vcab["Structural Coverage"]=="check V/C Annotation"] # These entries are probably Abs only containing V region.
struc_cov_unclassified=vcab.loc[vcab["Structural Coverage"]=="not classified"]
sus_v_c.to_csv("./unusual_cases/suspicious_v_c_annotation.csv")
unusual.to_csv("./unusual_cases/unusual.csv")


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

ff_vcab.to_csv("vcab.csv")

# 3. Generate files for the shiny app:
# 3.1. Generate POPSComp results
os.system("cd ../pops")
os.system("sh pops.sh ../pdb_struc/c_pdb/") # PDB structures with C region only are inputted for POPSComp analysis
os.system("cd -")

# 3.2. Generate BLAST databases
# generate fasta files
convert_seq_from_df_to_fasta(ff_vcab,'HV_seq',"../seq_db/vcab_db")
convert_seq_from_df_to_fasta(ff_vcab,'LV_seq',"../seq_db/vcab_db")
os.system("cat HV_seq.fasta LV_seq.fasta > all_v_seq.fasta") # Combining two fasta files together to generate the fasta file containing all sequences of V region
os.system("cat H_seq.fasta L_seq.fasta > all_full_seq.fasta") # Combining two fasta files together to generate the fasta file containing all sequences of V + C region (full sequence)
# make BLAST database. Note: BLAST should be installed on command line
os.system("makeblastdb -in ../seq_db/vcab_db/H_seq.fasta -dbtype prot")
os.system("makeblastdb -in ../seq_db/vcab_db/L_seq.fasta -dbtype prot")
os.system("makeblastdb -in ../seq_db/vcab_db/HV_seq.fasta -dbtype prot")
os.system("makeblastdb -in ../seq_db/vcab_db/LV_seq.fasta -dbtype prot")
os.system("makeblastdb -in ../seq_db/vcab_db/all_v_seq.fasta -dbtype prot")
os.system("makeblastdb -in ../seq_db/vcab_db/all_full_seq.fasta -dbtype prot")

# 3.3 Download all the PDB files, Cut pdb files (Done in previous steps)
