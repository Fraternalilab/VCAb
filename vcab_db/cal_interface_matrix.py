# To calculate the interface difference index and generate a mtrix holding this value
# Dongjun Guo, Apr.2022

import numpy as np
import pandas as pd
import Bio.PDB
from Bio import AlignIO
import os
import math
import time


# PART 1. read pops result_lstdef replace_back_T_F_chain (o_df):
def replace_back_T_F_chain (o_df):
    df=o_df.copy()
    for c,i in enumerate(df['Chain']):
        if i == 'TRUE':
            #print (c,i)
            df['Chain'][c]='T'
        if i == 'FALSE':
            #print (df['Chain'][c])
            df['Chain'][c]='F'

    return df

def read_pops_file(iden_code,pops_dir):
    # return popsFile: h_df, l_df, which contains residues with d_SASA >15
    h,l=iden_code.split('_')[1]

    o_df=pd.read_csv(f'{pops_dir}/{iden_code}_C_deltaSASA_rpopsResidue.txt',sep=' ')
    n_df=replace_back_T_F_chain (o_df)

    df=n_df.loc[o_df['D_SASA.A.2']>15] # only keep residues with D_SASA above the threshold

    h_df=df.loc[df['Chain']==h]
    l_df=df.loc[df['Chain']==l]

    h_df.reset_index(inplace=True)
    l_df.reset_index(inplace=True)
    return h_df,l_df

# Part 2. Run the pairwise profile aln
def run_pairwise_profile_aln(iden_code,chainType,ref_profile_dir,out_dir,df):
    # align the C_coor sequence to the reference profile
    ## Note: if the antibody is the Fab, ref_profile of CH1 seqs is used
    ## if the antibody is the whole antibody, ref_profile of the sequence covering all the domains is used
    ## in summary, the coverage of the ref_profile should corresponds to the structural coverage of the antibody
    # ref_profile: the name of the ref_profile
    # chainType can only be "H" or "L"

    col_name=""
    ref_profile=""

    struc_cov=df.loc[df["iden_code"]==iden_code,"Structural Coverage"].item()
    # assign different C_coor_seq and ref_profile to different chainType(H/L) and structural coverage (Fab/full antibody)
    if chainType.lower()=="h":
        col_name="HC_coordinate_seq"
        ref_profile=f"{ref_profile_dir}/unique_alleles_CH1_aln_profile"

        if struc_cov=='full antibody':
            ref_profile=f"{ref_profile_dir}/unique_alleles_aln_profile"

    if chainType.lower()=="l":
        col_name="LC_coordinate_seq"
        ref_profile=f"{ref_profile_dir}/unique_light_aln_profile"

    c_coor_seq=df.loc[df["iden_code"]==iden_code,col_name].item()

    seq_id=f"{iden_code}_{col_name}"
    with open(f"{out_dir}/{seq_id}.fasta","w") as f:
        f.write(f">{seq_id}\n")
        f.write(f"{c_coor_seq}\n")

    cmd_string=f"t_coffee {out_dir}/{seq_id}.fasta -profile={ref_profile} -outfile {out_dir}/{seq_id}_aln_to_ref.aln -in Mclustalw_pair"
    os.system(cmd_string)


# Part 3. Make sure the interface matrix generated are of the same size
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

def convert_all_domain_profile_to_CH1_profile (vcab_ad_aln_res_info,all_domain_profile,CH1_profile,seq_title):
    # Make sure the interface matrix generated from Fab and full_ig structure are of the same size
    # vcab_ad_aln_res_info is the ad_aln including all domains
    # seq_title is the type of the ref_seq. e.g. IGHE
    ad_aln=AlignIO.read(all_domain_profile,"clustal")
    ch1_aln=AlignIO.read(CH1_profile,"clustal")
    seq_in_ad_aln=str([ad_aln[i] for i in range(len(ad_aln)) if seq_title in ad_aln[i].id][0].seq)
    seq_in_ch1_aln=str([ch1_aln[i] for i in range(len(ch1_aln)) if seq_title in ch1_aln[i].id][0].seq)

    ad_aln_pos_to_ch1_aln={}
    # the dict contaning the conversion between the aln_pos (seq_in_ad_aln) and aln_pos (seq_in_ch1_aln)
    for k in range(len(seq_in_ad_aln)):
        # k is acturally the aln_pos in seq_in_ad_aln
        ad_aln_char=seq_in_ad_aln[k]
        if ad_aln_char !="-":
            ref_seq_pure_pos=map_aln_pos_to_pure_seq_pos (k,seq_in_ad_aln)
            ch1_aln_pos=map_pure_seq_pos_to_aln_pos (ref_seq_pure_pos,seq_in_ch1_aln)
            ad_aln_pos_to_ch1_aln[k]=ch1_aln_pos

    results={}
    n_vcab_res_info={k:v for k,v in vcab_ad_aln_res_info.items() if v!=[np.nan,np.nan]}
    for k in n_vcab_res_info.keys():
        if k in ad_aln_pos_to_ch1_aln.keys():
            # only consider the vcab residues which aligns to the ref_seq
            # ignore cases where the gap in ref_seq is generated when vcab_seq is aligned to ref_seq
            n_k=ad_aln_pos_to_ch1_aln[k]
            if n_k != None:
                # exclude the residues in the ad_aln but not in the
                results[n_k]=vcab_ad_aln_res_info[k]

    # fill the gaps information to the results dict
    for n_k in range(len(str(ch1_aln[0].seq))):
        if n_k not in results.keys():
            results[n_k]=[np.nan,np.nan]
    return dict(sorted(results.items()))


def map_profile_aln_to_residue_info (iden_code,chainType,ref_profile_dir,aln_out_dir,c_pdb_dir,pops_dir,df):
    # Collect the informations of residues/gaps of CH1/CL, which later would be the axis of the matrix
    # returned results is a dict, with items in this format:
    ## {aln_pos_in_pairwise_profile_aln:[res_obj,if_interface_residue]}
    ### if char at aln_pos is a gap, then res_obj=np.nan
    # aln is the pairwise_profile_aln of H_coor_seq/L_coor_seq to the corresponding_profile

    # Assign different values to different chainTypes (H/L)
    chainTypeNum=0
    aln_fn=f"{aln_out_dir}/{iden_code}_HC_coordinate_seq_aln_to_ref.aln"
    if chainType.lower()=="l":
        chainTypeNum=1
        aln_fn=f"{aln_out_dir}/{iden_code}_LC_coordinate_seq_aln_to_ref.aln"

    # Step 0. Generate pairwise profile aln file and read the file
    run_pairwise_profile_aln(iden_code,chainType,ref_profile_dir,aln_out_dir,df)
    aln=AlignIO.read(aln_fn,"clustal")

    # Step 1. Clean the aln: remove the gap-column generated in the ref_profile during the alignment
    ref_part=aln[1:] # The ref_profile part in the alignment
    ref_seq_num=len(ref_part) # The number of the sequence in the alignment of the ref_part (row_num)
    seq_len=len(ref_part[0].seq) # The seq_length(col_num) of in the aln==len(seq) in ref_part

    n_aln=aln[:,0:1] # initialize the new_aln with the first column of aln (which would be deleted later)
    for col in range(seq_len):
        ref_aln_col=ref_part[:,col]
        aln_col_with_id=aln[:,col:col+1]
        if ref_aln_col!="-"*ref_seq_num:
            # exclude the columns which contains all gaps in ref_part
            n_aln+=aln_col_with_id
    n_aln=n_aln[:,1:]

    # Step 2. Get the info stored in pdb file & pops results

    ## Acquire the chain_obj:
    parser=Bio.PDB.PDBParser()
    structure=parser.get_structure(iden_code,f"{c_pdb_dir}/{iden_code}_C.pdb")
    chainid=iden_code.split("_")[1][chainTypeNum]
    chain_obj=structure[0][chainid]

    ## Read the pops result (already filtered by the d_SASA cut-off)
    pops=read_pops_file(iden_code,pops_dir)[chainTypeNum]

    # Step 3. Collect the information & generate the result dict
    ## dict format: {aln_pos_in_pairwise_profile_aln:[res_obj,if_interface_residue]}
    result={}
    vcab_aln_seq=str(n_aln[0].seq)
    pure_pos_counter=0
    for k in range(len(vcab_aln_seq)):
        vcab_aln_char=vcab_aln_seq[k]
        if vcab_aln_char != "-":
            res_obj=list(chain_obj)[pure_pos_counter]
            pure_pos_counter+=1
            if_interface_res=len(pops.loc[(pops["ResidNe"]==res_obj.resname)&(pops["ResidNr"]==res_obj.id[1])])
            result[k]=[res_obj,if_interface_res]
        else:
            # at the gap position
            result[k]=[np.nan,np.nan]

    # Step 4. Make the dimension of the full_ig (H chain) the same as fab
    struc_cov=df.loc[df["iden_code"]==iden_code,'Structural Coverage'].item()
    if chainType.lower()=="h" and struc_cov=='full antibody':
        h_type=df.loc[df["iden_code"]==iden_code,"Htype"].item()
        h_type_dict={'IgA1':'IGHA1',
                     'IgA2':'IGHA2',
                     'IgD':'IGHD',
                     'IgG1':'IGHG1',
                     'IgG2':'IGHG2',
                     'IgG3':'IGHG3',
                     'IgG4':'IGHG4',
                     'IgM':'IGHM'}
        result=convert_all_domain_profile_to_CH1_profile (result,f"{ref_profile_dir}/unique_alleles_aln_profile",f"{ref_profile_dir}/unique_alleles_CH1_aln_profile",h_type_dict[h_type])

    return result

# Part 4. Calculate the distance matrix
def generate_distance_matrix(iden_code,ref_profile_dir,aln_out_dir,c_pdb_dir,pops_dir,df):
    h_res_info=map_profile_aln_to_residue_info (iden_code,"h",ref_profile_dir,aln_out_dir,c_pdb_dir,pops_dir,df)
    l_res_info=map_profile_aln_to_residue_info (iden_code,"l",ref_profile_dir,aln_out_dir,c_pdb_dir,pops_dir,df)

    result = np.zeros((len(h_res_info), len(l_res_info)), np.float)
    for row, h_res in h_res_info.items() :
        for col, l_res in l_res_info.items() :
            hres_obj,h_res_if_interface=h_res
            lres_obj,l_res_if_interface=l_res
            if (h_res_if_interface==1) and (l_res_if_interface==1):
                # would only give a value at this position when both residues are the interface residue.
                diff_vector  = hres_obj["CA"].coord - lres_obj["CA"].coord
                result[row, col] = np.linalg.norm(diff_vector)

    return result

# Part 5. Calculate interface distance index and collect the values into a matrix
def exclude_mtrx_not_calculated(o_df,fn):
    # fn: the file name of the file containing the iden_code with not calculated matrix(format: comma separated string) would be generated
    f=open(fn,"r")
    lst=f.read().split(",")
    f.close()
    df=o_df.loc[map(lambda x: x not in lst,o_df["iden_code"])]
    df=df.reset_index(drop=True)
    return df

def cal_distance_between_matrix (ab1,ab2,mtrx_dir):
    # both ab1 and ab2 are iden_code
    m1=np.loadtxt(f"{mtrx_dir}/{ab1}_interface_dist_mtrx.txt",dtype=float)
    m2=np.loadtxt(f"{mtrx_dir}/{ab2}_interface_dist_mtrx.txt",dtype=float)
    if m1.shape==m2.shape:
        item_diff=m1-m2
        item_diff_square=item_diff*item_diff
        item_diff_square_sum=item_diff_square.sum()
        dist=math.sqrt(item_diff_square_sum)
        return dist
    else:
        raise ValueError("Arrays must have the same size")

def generate_dm_of_interface_dm(df,mtrx_dir):
    result = np.zeros((len(df), len(df)), np.float)
    for r_counter,row in enumerate(df.index):
        ab1=df.loc[row,"iden_code"]
        for col in list(df.index)[r_counter+1:]:
            ab2=df.loc[col,"iden_code"]
            result[row,col]=cal_distance_between_matrix (ab1,ab2,mtrx_dir)
    labels={i:list(df.loc[i,["iden_code","Htype","Ltype"]].values) for i in df.index}
    return result,labels

######### APPLY THE FUNCTIONS ##############
vcab=pd.read_csv("./final_vcab.csv")
ref_profile_dir="../ch1_cl_interface_matrix/ref_profile/"
aln_out_dir="../ch1_cl_interface_matrix/aln_results/"
c_pdb_dir="../pdb_struc/c_pdb/"
pops_dir="../pops/result"

mtrx_out_dir="../ch1_cl_interface_matrix/matrix_results/"
os.system("sh ../ch1_cl_interface_matrix/check_folders.sh")

int_mtrx_not_calculated=[]
for i in vcab.index:
    iden_code=vcab.loc[i,"iden_code"]
    try:
        interface_mtrx=generate_distance_matrix(iden_code,ref_profile_dir,aln_out_dir,c_pdb_dir,pops_dir,vcab)
        np.savetxt(f"{mtrx_out_dir}/{iden_code}_interface_dist_mtrx.txt",interface_mtrx,fmt='%10.5f')
    except:
        int_mtrx_not_calculated.append(iden_code)

with open("../ch1_cl_interface_matrix/mtrx_not_calculated.txt", 'w') as f:
    f.write(",".join(int_mtrx_not_calculated))

# Calculate the distance matrix of dm (matrix of the interface distance index)
flt_vcab=exclude_mtrx_not_calculated(vcab,"../ch1_cl_interface_matrix/mtrx_not_calculated.txt")
# exclude the VCAb entries with no interface matrix (mainly because the POPSComp result is not available for the solution scattering method)
dm,ab_info_label=generate_dm_of_interface_dm(df,mtrx_out_dir)
np.savetxt(f"../ch1_cl_interface_matrix/dm_of_interface_dist_mtrx.txt",dm,fmt='%10.5f')

dm_df=pd.DataFrame(dm,columns=[ab_info_label[k][0] for k in ab_info_label.keys()])
dm_df.to_csv("../ch1_cl_interface_matrix/dm_of_interface_dist_mtrx.csv")

today=date.today()
update_date=today.strftime("%d/%m/%Y")
with open('VCAb_db_update_date.txt', 'w') as f:
    f.write(f"{update_date}")
