import numpy as np
import pandas as pd
import Bio.PDB
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.SeqUtils import seq1

import os
import math
from datetime import date
import json

aa_names= {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
          "-":"-"}

imgt_numbering=['1H',
 '1G',
 '1F',
 '1E',
 '1D',
 '1C',
 '1B',
 '1A',
 '1',
 '2',
 '3',
 '4',
 '5',
 '6',
 '7',
 '8',
 '9',
 '10',
 '11',
 '12',
 '13',
 '14',
 '15',
 '15A',
 '15B',
 '15C',
 '16',
 '17',
 '18',
 '19',
 '20',
 '21',
 '22',
 '23',
 '24',
 '25',
 '26',
 '27',
 '28',
 '29',
 '30',
 '31',
 '34',
 '35',
 '36',
 '37',
 '38',
 '39',
 '40',
 '41',
 '42',
 '43',
 '44',
 '45',
 '45A',
 '45B',
 '45C',
 '45D',
 '45E',
 '45F',
 '45G',
 '77',
 '78',
 '79',
 '80',
 '81',
 '82',
 '83',
 '84',
 '84A',
 '84B',
 '84C',
 '84D',
 '84E',
 '84F',
 '84G',
 '85G',
 '85F',
 '85E',
 '85D',
 '85C',
 '85B',
 '85A',
 '85',
 '86',
 '87',
 '88',
 '89',
 '90',
 '91',
 '92',
 '93',
 '94',
 '95',
 '96',
 '96A',
 '96B',
 '97',
 '98',
 '99',
 '100',
 '101',
 '102',
 '103',
 '104',
 '105',
 '106',
 '107',
 '108',
 '109',
 '110',
 '111',
 '112',
 '113',
 '114',
 '115',
 '116',
 '117',
 '118',
 '119',
 '120',
 '121',
 '122',
 '123',
 '124',
 '125',
 '126',
 '127',
 '128']

# 1. read pops result_lst
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

def read_pops_file(iden_code,pops_dir,if_c="_C"):
    """
    Return popsFile: h_df, l_df, which contains residues with d_SASA >15
    :args if_C: determines the POPS file readed.
                If "", read POPS for chain_pdb;
                If "_C", read POPS for c_PDB (the pdb containing C region only)
    """
    h,l=iden_code.split('_')[1]

    o_df=pd.read_csv(f'{pops_dir}/{iden_code}{if_c}_deltaSASA_rpopsResidue.txt',sep=' ')
    n_df=replace_back_T_F_chain (o_df)

    df=n_df.loc[o_df['D_SASA.A.2']>15] # only keep residues with D_SASA above the threshold

    h_df=df.loc[df['Chain']==h]
    l_df=df.loc[df['Chain']==l]

    h_df.reset_index(inplace=True)
    l_df.reset_index(inplace=True)
    return h_df,l_df

# 2. Convert positions between pure_seq & alned_seq (seq with gaps)
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

# 3. Convert IMGT numbering(author_seq) into pure_position(coor_seq)

def map_imgt_numbering_to_residue_info(iden_code,chainType,num_df,pdb_dir,pops_dir,if_C="_C",num_scheme=imgt_numbering,gap_included=False):
    """
    Returned a list in this format {imgt_numbering:[res_obj, if_interface_residue]}
    :args chainType: can only be "H" or "L"
    :args num_df: the df outputed by anarci_c containing the C_numbering results
    :args pdb_dir: the directory of the c_pdb files
    :args if_C: determines the structural & POPS file readed. If "", read chain_PDB; If "_C", read c_PDB (the pdb file contains C region only).
    :args num_scheme: None or the list containing all the numbering we want to included into the interface matrix
                          If the num_scheme ==None, then num_scheme== existed_numbering
    """
    pdbid=iden_code.split("_")[0]
    chainTypeNum=0
    if chainType.lower()=="l":
        chainTypeNum=1
    chain=iden_code.split("_")[1][chainTypeNum]
    id_code=f"{pdbid}_{chain}"

    #num_info=num_df.loc[num_df["Id"]==id_code]
    num_info=num_df.loc[map(lambda x: x[0:6]==id_code,num_df["Id"])]
    pure_num_info0=[(numbering,val.values[0]) for (numbering,val) in num_info.iloc[:,13:].items() if (val.values[0] not in["deleted","-"])]
    # pure_num_info0 removes positions generating gaps(both "deleted" and. "-") in the coor_seq
    pure_num_info={numbering:[pure_pos,res] for pure_pos,(numbering,res) in enumerate(pure_num_info0)}

    #return pure_num_info
    existed_numbering=list(pure_num_info.keys())
    num_seq_frag="".join([i[1] for i in list(pure_num_info.values())])

    # Get the structural_object of the chain
    parser=Bio.PDB.PDBParser()
    structure=parser.get_structure(iden_code,f"{pdb_dir}/{iden_code}{if_C}.pdb")

    chain_obj=structure[0][chain]

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

    # Read POPS File
    pops=read_pops_file(iden_code,pops_dir,if_C)[chainTypeNum]

    # 2. Convert IMGT_num_pos into pure_pos of coor_seq
    # IMGT_num --> pure_pos (imgt_seq) -->aln_pos (imgt_seq)=aln_pos(coor_seq) -->pure_pos (coor_seq)
    result={}
    #test=[]
    if num_scheme ==None:
        num_scheme=existed_numbering

    for i in num_scheme:
        result[i]=[np.nan,np.nan]

        if i in existed_numbering:
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
                result[i]=[np.nan,np.nan]

            else:
                # pure_position:
                coor_seq_pure_pos=map_aln_pos_to_pure_seq_pos(coor_seq_aln_pos,alned_coor_seq)

                res_obj=coor_seq_info[coor_seq_pure_pos]
                if_interface_res=len(pops.loc[(pops["ResidNe"]==res_obj.resname)&(pops["ResidNr"]==res_obj.id[1])])
                result[i]=[res_obj,if_interface_res]
    if gap_included:
        gap_info=[(numbering,val.values[0]) for (numbering,val) in num_info.iloc[:,13:].items() if (val.values[0]=="-")]
        #record the original numbering order (gap included, don't include "deleted")
        gap_included_num=[num for (num,res) in num_info.iloc[:,13:].items() if res.values[0]!="deleted"]

        for num,gap in gap_info:
            result[num]=["-",np.nan]
        result={k:result[k] for k in gap_included_num}
    return result
    #return test


# 4. Generate the distance matrix
def generate_distance_matrix (iden_code,pdb_dir,pops_dir,hcnum,lcnum):

    h_res_info=map_imgt_numbering_to_residue_info(iden_code,"h",hcnum,pdb_dir,pops_dir)
    l_res_info=map_imgt_numbering_to_residue_info(iden_code,"l",lcnum,pdb_dir,pops_dir)

    result = np.zeros((len(h_res_info), len(l_res_info)), float)
    for row, h_res in enumerate(h_res_info.values()) :
        for col, l_res in enumerate(l_res_info.values()) :
            hres_obj,h_res_if_interface=h_res
            lres_obj,l_res_if_interface=l_res
            if (h_res_if_interface==1) and (l_res_if_interface==1):
                # would only give a value at this position when both residues are the interface residue.
                diff_vector  = hres_obj["CA"].coord - lres_obj["CA"].coord
                result[row, col] = np.linalg.norm(diff_vector)

    return result

# 5. Calculate interface distance index and collect the values into a matrix
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
    if os.path.exists(f"{mtrx_dir}/{ab1}_interface_dist_mtrx.txt") is False or \
        os.path.exists(f"{mtrx_dir}/{ab2}_interface_dist_mtrx.txt") is False:
        return np.nan
    m1=np.loadtxt(f"{mtrx_dir}/{ab1}_interface_dist_mtrx.txt",dtype=float)
    m2=np.loadtxt(f"{mtrx_dir}/{ab2}_interface_dist_mtrx.txt",dtype=float)
    assert m1.shape==m2.shape,"Arrays must have the same size"

    item_diff=m1-m2
    item_diff_square=item_diff*item_diff
    item_diff_square_sum=item_diff_square.sum()
    dist=math.sqrt(item_diff_square_sum)
    return dist

def generate_dm_of_interface_dm(df,mtrx_dir):

    result = np.zeros((len(df), len(df)), float)
    for r_counter,row in enumerate(df.index):
        ab1=df.loc[row,"iden_code"]
        for col in list(df.index)[r_counter+1:]:
            # So that it will only calculate half of the dm matrix, because the matrix is symmetrical
            ab2=df.loc[col,"iden_code"]
            result[row,col]=cal_distance_between_matrix (ab1,ab2,mtrx_dir)
    labels={i:list(df.loc[i,["iden_code","Htype","Ltype"]].values) for i in df.index}
    return result,labels

def update_dm_of_interface_dm(df,old_mtrx_df_fn,mtrx_dir):
    # df: vcab, old_mtrx_df_fn: the csv file name of old dm_of_dm
    result = np.zeros((len(df), len(df)), float)
    old_mtrx_df=pd.read_csv(old_mtrx_df_fn).drop(columns=["Unnamed: 0"])

    old_mtrx=old_mtrx_df.to_numpy()

    # incorporate the old matrix as the part of new results
    result[:old_mtrx.shape[0],:old_mtrx.shape[1]]=old_mtrx

    # Find the newly added VCAb entries
    old_abs=list(old_mtrx_df.columns)
    new_abs=list(df["iden_code"].values)
    added_abs=[i for i in new_abs if i not in old_abs]

    new_abs_ordered=old_abs+added_abs # order the abs: old_abs first, then added_abs

    for row,ab1 in enumerate(new_abs_ordered):
        for col,ab2 in enumerate(new_abs_ordered):
            if row < len(old_abs) and col < len(old_abs):
                # skip the values in old matrix
                continue
            if col > row:
                result[row,col]=cal_distance_between_matrix (ab1,ab2,mtrx_dir)

    labels={i:[new_abs_ordered[i]]+list(df.loc[df["iden_code"]==new_abs_ordered[i],["Htype","Ltype"]].values[0]) for i in range(len(new_abs_ordered))}
    return result,labels

######### FUNCTIONS TO GENERATE RES_INFO & COMPLETE DISTANCE MATRIX ###########
def combine_v_c_numbering (v_res_info,c_res_info):
    v_res_df=pd.DataFrame(v_res_info,index=["residue_obj","if_interface"]).T
    c_res_df=pd.DataFrame(c_res_info,index=["residue_obj","if_interface"]).T

    v_res_df=v_res_df.reset_index().rename(columns={"index":"imgt_numbering"})
    c_res_df=c_res_df.reset_index().rename(columns={"index":"imgt_numbering"})

    res_df=pd.concat([v_res_df,c_res_df],keys=["V","C"]).reset_index().rename(columns={"level_0":"VorC"}).drop(columns=["level_1"])

    recorded=[]
    extra=[]
    # This is drop the V region end tail numbered in the C region
    for i in res_df.index:
        res=res_df.loc[i,"residue_obj"]
        if (res in recorded or type(res)==float) and type(res)!=str:
            extra.append(i)
        recorded.append(res)
    res_df=res_df.drop(index=extra).reset_index(drop=True)

    res_df["residue"]=[r if type(r)==str else r.resname for r in res_df["residue_obj"].values]
    res_df["res_code"]=[aa_names[r] if r in aa_names.keys() else "" for r in res_df["residue"].values]
    res_df["pdb_numbering"]=[pd.NA if type(r)==str else r.id[1] for r in res_df["residue_obj"].values]

    return res_df

def generate_complete_distance_matrix (iden_code,chain_pdb_dir,total_pops_dir,vhnum,vlnum,hcnum,lcnum,m_out_dir,n_out_dir):
    """
    Calculate the complete domain-packing contact matrix
    1. Include all the residues (non-gap)
    2. Include both V and C domains
    3. Calculate the paiwise distance for all, no matter if they are
    involved in interface or not
    4. Meanwhile, generate the dataframe, recording the info for each residue (gap included)
       (pdb numbering, imgt numbering, if it is involved in interface)
    :args m_out_dir: the out directory storing the total matrix
    :args n_out_dir: the our directory storing the residues info with imgt,pdb, VorC, if_interface info.

    """

    vh_res_info=map_imgt_numbering_to_residue_info(iden_code,"h",vhnum,chain_pdb_dir,total_pops_dir,if_C="",num_scheme=None,gap_included=True)
    vl_res_info=map_imgt_numbering_to_residue_info(iden_code,"l",vlnum,chain_pdb_dir,total_pops_dir,if_C="",num_scheme=None,gap_included=True)

    ch_res_info=map_imgt_numbering_to_residue_info(iden_code,"h",hcnum,chain_pdb_dir,total_pops_dir,if_C="",num_scheme=None,gap_included=True)
    cl_res_info=map_imgt_numbering_to_residue_info(iden_code,"l",lcnum,chain_pdb_dir,total_pops_dir,if_C="",num_scheme=None,gap_included=True)

    h_res_df=combine_v_c_numbering (vh_res_info,ch_res_info)
    l_res_df=combine_v_c_numbering (vl_res_info,cl_res_info)

    # take out the gaps to prepare for the calculation of complete interface matrix:
    h_res_df1=h_res_df.loc[map(lambda x:type(x)!=str,h_res_df["residue_obj"].values)].reset_index(drop=True)
    l_res_df1=l_res_df.loc[map(lambda x:type(x)!=str,l_res_df["residue_obj"].values)].reset_index(drop=True)

    result = np.zeros((len(h_res_df1), len(l_res_df1)), float)
    for row, h_res in enumerate(h_res_df1["residue_obj"].values) :
        for col, l_res in enumerate(l_res_df1["residue_obj"].values) :
            hres_obj=h_res
            lres_obj=l_res
            try:
                diff_vector  = hres_obj["CA"].coord - lres_obj["CA"].coord
            except:
                diff_vector = np.nan
            #    print ((hres_obj.id,lres_obj.id))
            result[row, col] = np.linalg.norm(diff_vector)

    np.savetxt(f"{m_out_dir}/{iden_code}.txt",result,fmt='%10.5f')

    h_res_df=h_res_df.drop(columns=["residue_obj"])
    l_res_df=l_res_df.drop(columns=["residue_obj"])
    h_res_df.to_csv(f"{n_out_dir}/{iden_code}_H_res_info.csv", index=False)
    l_res_df.to_csv(f"{n_out_dir}/{iden_code}_L_res_info.csv", index=False)

    return result,h_res_df,l_res_df

######### APPLY THE FUNCTIONS ##############
vcab=pd.read_csv("./result/final_vcab.csv").drop(columns=["Unnamed: 0"])

vhnum=pd.read_csv("./result/num_result/vnumbering_H.csv")
vlnum=pd.read_csv("./result/num_result/vnumbering_KL.csv")
if "Unnamed: 0" in vhnum.columns:
    vhnum=vhnum.drop(columns=["Unnamed: 0"])
if "Unnamed: 0" in vlnum.columns:
    vlnum=vlnum.drop(columns=["Unnamed: 0"])

hcnum=pd.read_csv("./result/num_result/cnumbering_H_C1.csv")
lcnum=pd.read_csv("./result/num_result/cnumbering_KL_C.csv")
if "Unnamed: 0" in hcnum.columns:
    hcnum=hcnum.drop(columns=["Unnamed: 0"])
if "Unnamed: 0" in lcnum.columns:
    lcnum=lcnum.drop(columns=["Unnamed: 0"])


pdb_dir="../pdb_struc/c_pdb/"
pops_dir="../pops/result"

chain_pdb_dir="../pdb_struc/chain_pdb/"
total_pops_dir="../pops/total_result/"

mtrx_out_dir="../ch1_cl_interface_matrix/matrix_results/"
complete_mtrx_dir="../ch1_cl_interface_matrix/complete_matrix"
res_info_out_dir="../res_info"
#os.system("sh ../ch1_cl_interface_matrix/check_folders.sh")

print ("Calculating interface matrix")
int_mtrx_not_calculated=[]
for i in vcab.index:
    iden_code=vcab.loc[i,"iden_code"]
    if os.path.exists(f"{mtrx_out_dir}/{iden_code}_interface_dist_mtrx.txt"):
        # skip the precalcuated interface matrix
        continue

    try:
        interface_mtrx=generate_distance_matrix(iden_code,pdb_dir,pops_dir,hcnum,lcnum)
        np.savetxt(f"{mtrx_out_dir}/{iden_code}_interface_dist_mtrx.txt",interface_mtrx,fmt='%10.5f')
    except:
        print (iden_code)
        int_mtrx_not_calculated.append(iden_code)

with open("../ch1_cl_interface_matrix/mtrx_not_calculated.txt", 'w') as f:
    f.write(",".join(int_mtrx_not_calculated))

print ("Calculating distance matrix of interface matrices")
# Calculate the distance matrix of dm (matrix of the interface distance index)
flt_vcab=exclude_mtrx_not_calculated(vcab,"../ch1_cl_interface_matrix/mtrx_not_calculated.txt")
# exclude the VCAb entries with no interface matrix (mainly because the POPSComp result is not available for the solution scattering method)

dm,ab_info_label=("","")
# if file exists, update the matrix; else, calculate the matrix from scratch:
if os.path.exists("../ch1_cl_interface_matrix/dm_of_interface_dist_mtrx.csv"):
    dm,ab_info_label=update_dm_of_interface_dm(vcab,"../ch1_cl_interface_matrix/dm_of_interface_dist_mtrx.csv",mtrx_out_dir)
else:
    dm,ab_info_label=generate_dm_of_interface_dm(vcab,mtrx_out_dir)
np.savetxt(f"../ch1_cl_interface_matrix/dm_of_interface_dist_mtrx.txt",dm,fmt='%10.5f')

print ("Storing matrix into csv")
with open('../ch1_cl_interface_matrix/dm_of_dm_label.json', 'w') as fp:
    json.dump(ab_info_label, fp)

dm_df=pd.DataFrame(dm,columns=[ab_info_label[k][0] for k in ab_info_label.keys()],index=[ab_info_label[k][0] for k in ab_info_label.keys()])
dm_df.to_csv("../ch1_cl_interface_matrix/dm_of_interface_dist_mtrx.csv")

# Generate res_info
print ("calculating res info")
total_int_mtrx_not_calculated=[]
for i in vcab.index:
    iden_code=vcab.loc[i,"iden_code"]
    if os.path.exists(f"{n_out_dir}/{iden_code}_H_res_info.csv"):
        # skip the precalcuated res_info
        continue

    try:
        __,__,__,=generate_complete_distance_matrix (iden_code,chain_pdb_dir,total_pops_dir,vhnum,vlnum,hcnum,lcnum,complete_mtrx_dir,res_info_out_dir)

    except:
        print (iden_code)
        total_int_mtrx_not_calculated.append(iden_code)
with open("../ch1_cl_interface_matrix/complete_matrix/mtrx_not_calculated.txt", 'w') as f:
    f.write(",".join(total_int_mtrx_not_calculated))

today=date.today()
update_date=today.strftime("%d/%m/%Y")
with open('VCAb_db_update_date.txt', 'w') as f:
    f.write(f"{update_date}")
