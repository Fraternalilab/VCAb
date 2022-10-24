import Bio.PDB
from Bio import SeqIO
from Bio import AlignIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

import numpy as np
import pandas as pd
from itertools import groupby

import os

from Bio.SeqUtils import seq1
from Bio.Seq import Seq

def torsion_angle(p):
    # Return the torsion angle of four points in 3D space
    # p: the list containing four points (x,y,z)
    # From: https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    b1 = b1/np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    angle=np.degrees(np.arctan2(y, x)) # [-pi,pi]
    if angle <0:
        angle+=360
    return angle

def get_frag (chain, start, end):
    # Return the fragment within seq_num range[start,end]
    # Chain is an Bio.PDB.Chain.Chain object
    # start & end defines the range for the fragment of interest
    # IMPORTANT:start & end are the PDB positions
    # Both residues at the start position and the end position are included
    frag=Bio.PDB.Chain.Chain("n")
    for res in chain:
        if res.id[1] in list(range(start,end+1)) and res.id[0]==' ':
            frag.add(res)
    return frag

# Converting positions between pure_seq and aln_seq:
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

# Extracting DSSP tables & add IMGT C_numbering info into the table
def extract_table_start (f_name):
    with open(f_name) as f:
        for i,line in enumerate(f):
            if '  #  RESIDUE AA' in line:
                return i
        f.close()

def extract_dssp_table_add_numbering (iden_code,dssp_dir,hcnum,lcnum):
    """
    Return the dssp table covering all the availble residues in the pdb file
    & add the numbering info into the table
    """
    # 1. Extract dssp table
    f_name=f"{dssp_dir}/{iden_code}.dssp"
    start=extract_table_start (f_name)
    f=open(f_name)
    result_lst=[]
    for i,line in enumerate(f.readlines()):
        if i > start:
            # The detailed explanation of dssp results:
            # https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html
            resnum_str=line[5:11] # residue number in pdb file, not the sequential number
            if resnum_str == "      ":
                # This always happens at the chain break, where the column for residue number is empty
                continue

            resnum_pdb=resnum_str.strip()
            chain_id=line[11] # the name of the chain id
            residue=line[13] # the name of the residue
            ss=line[16] # secondary structure

            result_lst.append([resnum_pdb,chain_id,residue,ss])

    result_df=pd.DataFrame(result_lst,columns=["resnum_pdb","chain_id","residue","2nd_struc"])
    pure_pdb_num=[]
    for i in list(result_df["resnum_pdb"].values):
        try:
            pure_pdb_num.append(int(i))
        except:
            pure_pdb_num.append(int(i[0:-1]))
    result_df["pure_resnum_pdb"]=pure_pdb_num

    nf_name=f_name.split("/")[-1]
    pure_f_name=nf_name.split(".dssp")[0]

    result_df["struc_name"]=pure_f_name

    # 2. Add the IMGT C_numbering result
    pdbid=iden_code.split("_")[0]
    h,l=iden_code.split("_")[1]

    def add_numbering_info(chain,num_df):
        dssp=result_df.loc[result_df["chain_id"]==chain].reset_index(drop=True)
        dssp["pure_pos"]=dssp.index
        dssp_seq="".join(dssp["residue"])

        id_code=f"{pdbid}_{chain}"
        num_info0=num_df.loc[num_df["Id"]==id_code]
        num_info=[(col,val.values[0]) for (col,val) in num_info0.iloc[:,13:].items() if (val.values[0] !="deleted")]

        pure_num_info=[(num,res) for (num,res) in num_info if res !="-"]
        num_seq_frag="".join([i[1] for i in pure_num_info])
        # 2.1 Perform the pairwise alignment
        matrix = matlist.blosum62
        pair_aln=pairwise2.align.globalds(num_seq_frag, dssp_seq, matrix,-10,-0.5,penalize_end_gaps=False)[0]
        # use blosum62 as the matrix, gap penalties follow the default value of EMBOSS needle
        # just take the first one as the aln result
        #return pair_aln
        alned_num_seq=pair_aln.seqA
        alned_dssp_seq=pair_aln.seqB

        # 2.2 Matching the pure_pos of dssp_seq with IMGT_num_pos:
        # pure_pos (dssp_seq) -->aln_pos(dssp_seq)=aln_pos(imgt_num_seq)-->pure_pos(imgt_num_seq)-->imgt_numbering
        imgt_num_lst=[]
        for i in dssp.index:
            pure_pos=dssp.loc[i,"pure_pos"]

            aln_pos=map_pure_seq_pos_to_aln_pos(pure_pos,alned_dssp_seq)
            if alned_num_seq[aln_pos]=="-":
                imgt_num_lst.append(None)
                continue
            imgt_seq_pure_pos=map_aln_pos_to_pure_seq_pos(aln_pos,alned_num_seq)

            imgt_num=pure_num_info[imgt_seq_pure_pos][0]
            imgt_num_lst.append(imgt_num)

        dssp["imgt_c_num"]=imgt_num_lst

        num=[]
        insertion=[]
        for i in dssp["imgt_c_num"].values:
            if type(i)!= str:
                num.append(0)
                insertion.append("")
            else:
                try:
                    num.append(int(i))
                    insertion.append("")
                except:
                    num.append(int(i[0:-1]))
                    insertion.append(i[-1])
        dssp["imgt_c_num_int"]=num
        dssp["insertion"]=insertion
        return dssp

    h_dssp=add_numbering_info(h,hcnum)
    l_dssp=add_numbering_info(l,lcnum)

    return h_dssp,l_dssp

def get_loops (chain,loop_pdb_pos,new_chain_id="n"):
    """
    Return the fragment as a chain object.
    The fragment can be both the loop of CH1/Cl and linker between V and C region
    :args chain: an Bio.PDB.Chain.Chain object
    :args loop_pdb_pos: the list containing all the pdb positions of the loop defined by IMGT
    """

    frag=Bio.PDB.Chain.Chain(new_chain_id)
    pure_pos_counter=0
    for res in chain:
        res_type,res_pdb0,res_pdb1=res.id
        res_pdb=f"{res_pdb0}{res_pdb1}".replace(" ","")

        for i in loop_pdb_pos:
            if res_pdb ==i:
                frag.add(res)

    return frag

################## 1. Calculate the CH1-CL interface angle ##################
def new_find_pdb_num_of_adjacent_E (dssp,num):
    """
    Find the pdb number of the adjacent point of the residue with IMGT numbering "num"
    : args num: the IMGT numbering
    1. For the extraction of the CH1/CL loop
     Use IMGT loop defination as the anchoring point
     IMGT loop defination (including both ends): [15.1,15.3],[45.1-45.7],[96.1-2]
     Thus, num should be 15,45,96
    """
    sub_df=pd.DataFrame([])
    if num in list(dssp["imgt_c_num_int"].values):
        sub_df=dssp.loc[(dssp["imgt_c_num_int"]==num)]
        #return loop_df
    else:
        sub_series_before=dssp.loc[(dssp["imgt_c_num_int"]<num)&(dssp["imgt_c_num_int"]!=0)].iloc[-1,:]
        sub_series_after=dssp.loc[(dssp["imgt_c_num_int"]>num)].iloc[0,:]
        sub_df=pd.DataFrame([sub_series_before,sub_series_after])
        sub_df=sub_df.loc[sub_df["2nd_struc"]!="E"]

    # return sub_df
    previous_E_index=dssp.loc[(dssp.index<list(sub_df.index)[0])&(dssp["2nd_struc"]=="E")].index[-1]
    next_E_index=dssp.loc[(dssp.index>list(sub_df.index)[-1])&(dssp["2nd_struc"]=="E")].index[0]

    loop_df=dssp.iloc[previous_E_index+1:next_E_index,:]
    loop_pdb=list(loop_df["resnum_pdb"].values)
    #return loop_df
    return loop_pdb


def new_get_CH1_CL_angle(hchain_obj,lchain_obj,h_dssp,l_dssp):
    # return the angle between CH1-CL, defined by the center of mass (COM)

    # Get pdb positions:

    ch1_loop=new_find_pdb_num_of_adjacent_E(h_dssp,15)+new_find_pdb_num_of_adjacent_E(h_dssp,45)+new_find_pdb_num_of_adjacent_E(h_dssp,96)
    cl_loop=new_find_pdb_num_of_adjacent_E(l_dssp,15)+new_find_pdb_num_of_adjacent_E(l_dssp,45)+new_find_pdb_num_of_adjacent_E(l_dssp,96)

    ch1_start=int(h_dssp.loc[map(lambda x: x != None,h_dssp["imgt_c_num"])].iloc[0,:]["resnum_pdb"])
    ch1_end=int(h_dssp.loc[map(lambda x: x != None,h_dssp["imgt_c_num"])].iloc[-1,:]["resnum_pdb"])

    cl_start=int(l_dssp.loc[map(lambda x: x != None,l_dssp["imgt_c_num"])].iloc[0,:]["resnum_pdb"])
    cl_end=int(l_dssp.loc[map(lambda x: x != None,l_dssp["imgt_c_num"])].iloc[-1,:]["resnum_pdb"])



    # Get the domain object
    ch1=get_frag (hchain_obj, ch1_start, ch1_end)
    cl=get_frag (lchain_obj, cl_start, cl_end)

    # Get the loops fragment:
    ch1_loop=get_loops (hchain_obj,ch1_loop)
    cl_loop=get_loops (lchain_obj,cl_loop)


    # Acquire center of mass (com) of these fragments
    ch1_com=ch1.center_of_mass()
    cl_com=cl.center_of_mass()
    ch1_loop_com=ch1_loop.center_of_mass()
    cl_loop_com=cl_loop.center_of_mass()
    #return (ch1_loop_com,cl_loop_com)
    #return (ch1_com,cl_com,ch1_loop_com,cl_loop_com)
    return torsion_angle([cl_loop_com,cl_com,ch1_com,ch1_loop_com])


################## 2. Calculate the elbow angle ##################
def get_linker_domain_region(dssp):
    vc_boundary=dssp.loc[map(lambda x: x == "1D",dssp["imgt_c_num"])].iloc[0,:]["pure_resnum_pdb"]
    c1_end=dssp.loc[map(lambda x: x != None,dssp["imgt_c_num"])].iloc[-1,:]["pure_resnum_pdb"]

    previous_E_index=dssp.loc[(dssp["pure_resnum_pdb"]<vc_boundary)&(dssp["2nd_struc"]=="E")].index[-1]
    next_E_index=dssp.loc[(dssp["pure_resnum_pdb"]>vc_boundary)&(dssp["2nd_struc"]=="E")].index[0]

    linker_df=dssp.iloc[previous_E_index+1:next_E_index,:]
    linker_pdb=list(linker_df["resnum_pdb"].values)
    return linker_pdb,vc_boundary,c1_end

def new_calc_elbow_angle(hchain_obj,lchain_obj,h_dssp,l_dssp):

    h_vc_linker_pos,h_vcb,ch1_end=get_linker_domain_region(h_dssp)
    l_vc_linker_pos,l_vcb,cl_end=get_linker_domain_region(l_dssp)

    vh=get_frag (hchain_obj, 0, h_vcb-1)
    ch1=get_frag (hchain_obj, h_vcb, ch1_end) # could lead to problem if it is a full ig

    vl=get_frag(lchain_obj,0,l_vcb-1)
    cl=get_frag (lchain_obj, l_vcb, cl_end)

    # Get the V & C fragment
    v_total=Bio.PDB.Chain.Chain("v")
    c_total=Bio.PDB.Chain.Chain("c")

    for res in vh:
        v_total.add(res)
    for res in vl:
        id_lst=list(res.id)
        res.id=(id_lst[0],id_lst[1],id_lst[2]+".l")
        # change the id of residue, because Bio.PDB would not add two residue of the same id into one chain
        v_total.add(res)

    for res in ch1:
        c_total.add(res)
    for res in cl:
        id_lst=list(res.id)
        res.id=(id_lst[0],id_lst[1],id_lst[2]+".l")
        c_total.add(res)

    # Get the V-C linker fragment:
    h_vc_linker=get_loops (hchain_obj,h_vc_linker_pos,"h")

    l_vc_linker_pos=[i+".l" for i in l_vc_linker_pos]
    #return lchain_obj,l_vc_linker_pos
    l_vc_linker=get_loops (lchain_obj,l_vc_linker_pos,"l")

    #return v_total,c_total,h_vc_linker,l_vc_linker
    # Acquire center of mass (com) of these fragments
    v_com=v_total.center_of_mass()
    c_com=c_total.center_of_mass()
    h_vc_linker_com=h_vc_linker.center_of_mass()
    l_vc_linker_com=l_vc_linker.center_of_mass()
    return torsion_angle([v_com,h_vc_linker_com,l_vc_linker_com,c_com])

################## 3. Calculate the CH1-CL interface angle and elbow angle together ##################

def calc_ch1_cl_and_elbow_angle (iden_code,dssp_dir,pdb_dir,hcnum,lcnum):
    # Acquire the structure:
    parser=Bio.PDB.PDBParser()
    structure=parser.get_structure(iden_code,f"{pdb_dir}/{iden_code}.pdb")

    # Acquire fragments in order to acquire COM to calculate the angles
    h,l=iden_code.split("_")[1]

    hchain_obj=structure[0][h]
    lchain_obj=structure[0][l]

    h_dssp,l_dssp=extract_dssp_table_add_numbering (iden_code,dssp_dir,hcnum,lcnum)
    #return hchain_obj,lchain_obj,h_dssp,l_dssp
    try:
        ch1_cl_angle=new_get_CH1_CL_angle(hchain_obj,lchain_obj,h_dssp,l_dssp)
    except:
        ch1_cl_angle=np.nan

    try:
    #if True:
        elbow_angle=new_calc_elbow_angle(hchain_obj,lchain_obj,h_dssp,l_dssp)
    except:
        elbow_angle=np.nan
    return ch1_cl_angle,elbow_angle

def combine_csvs(rf,cf):
    """
    For the update of csv files
    check for the files presented in both result folder and current folder
    1. Combine the common files in rf and cf into one, stored in rf.
    :args rf: result folder: stores the old data
    :args cf: current folder: stores the updated data(only the newly added part)
    returns nothing
    """
    common_files=[fn for fn in os.listdir(rf) if os.path.exists(f"{cf}/{fn}") and (".csv" in fn)]
    for f in common_files:
        df1=pd.read_csv(f"{rf}/{f}").drop(columns=["Unnamed: 0"])
        df2=pd.read_csv(f"{cf}/{f}").drop(columns=["Unnamed: 0"])
        df3=pd.concat([df1,df2]).reset_index(drop=True)
        if "Unnames: 0" in df3.columns:
            df3=df3.drop(columns=["Unnamed: 0"])

        df3.to_csv(f"{rf}/{f}")

if __name__=="__main__":
    vcabfn="./vcab.csv"
    hcnum_fn="./num_result/cnumbering_H_C1.csv"
    lcnum_fn="./num_result/cnumbering_KL_C.csv"

    dssp_dir="../dssp/dssp_results"
    pdb_dir="../pdb_struc/chain_pdb"

    vcab=pd.read_csv(vcabfn).drop(columns=["Unnamed: 0"])
    hcnum=pd.read_csv(hcnum_fn)
    lcnum=pd.read_csv(lcnum_fn)
    # print ("running dssp analysis")
    #os.system("sh ../dssp/run_dssp.sh")

    print("calculating angles")
    elbow_angles=[]
    ch1_cl_angles=[]
    for i in vcab.index:
        iden_code=vcab.loc[i,"iden_code"]
        try:
        #if True:
            c_angle,e_angle=calc_ch1_cl_and_elbow_angle (iden_code,dssp_dir,pdb_dir,hcnum,lcnum)
        except:
            c_angle,e_angle=(np.nan,np.nan)
        elbow_angles.append(e_angle)
        ch1_cl_angles.append(c_angle)

    vcab["elbow_angle"]=elbow_angles
    vcab["CH1-CL_interface_angle"]=ch1_cl_angles
    #vcab=vcab.rename(columns={'heavy_subclass':"heavy_vfamily",'light_subclass':"light_vfamily"})

    vcab.to_csv("./final_vcab.csv")
    combine_csvs("./result",".")
    combine_csvs("./result/num_result/","./num_result/")
    combine_csvs("./result/blast_result/","./blast_result/")
    combine_csvs("./result/blast_result/flt_bl_result/","./blast_result/flt_bl_result/")
