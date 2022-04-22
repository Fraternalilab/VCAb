# To calculate the elbow angles and the CH1-CL angle
# Dongjun Guo, Apr.2022

import Bio.PDB
from Bio import AlignIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio import SeqIO

import numpy as np
import pandas as pd
import json
from itertools import groupby

import os

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

def extract_table_start (f_name):
    with open(f_name) as f:
        for i,line in enumerate(f):
            if '  #  RESIDUE AA' in line:
                return i
        f.close()

def extract_dssp_table (f_name):
    # return the dssp table covering all the availble residues in the pdb file
    start=extract_table_start (f_name)
    f=open(f_name)
    result_lst=[]
    for i,line in enumerate(f.readlines()):
        if i > start:
            # The detailed explanation of dssp results:
            # https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html
            resnum_str=line[5:10] # residue number in pdb file, not the sequential number
            if resnum_str == "     ":
                # This always happens at the chain break, where the column for residue number is empty
                continue

            resnum_pdb=int(resnum_str)
            chain_id=line[11] # the name of the chain id
            residue=line[13] # the name of the residue
            ss=line[16] # secondary structure

            result_lst.append([resnum_pdb,chain_id,residue,ss])

    result_df=pd.DataFrame(result_lst,columns=["resnum_pdb","chain_id","residue","2nd_struc"])
    nf_name=f_name.split("/")[-1]
    pure_f_name=nf_name.split(".dssp")[0]

    result_df["struc_name"]=pure_f_name

    return result_df

def find_pdb_vcb_in_HandL(iden_code,df):
    ch1_start=int(df.loc[df["iden_code"]==iden_code,"pdb_H_VC_Boundary"].item())
    cl_start=int(df.loc[df["iden_code"]==iden_code,"pdb_L_VC_Boundary"].item())
    return ch1_start,cl_start

# PART 1. CH1-CL torsion angle calculation
# Strategy to extract the loops: find the adjacent beta-strand (represented by "E" in DSSP results) of the loop "anchoring" (ref) points
#### Anchor points (pure_pos) of the ref_seq has been precalculated and stored in .json file
#### Note: for the loop2 in CL (strand C-D), in order to escape the unusal beta-strand given by DSSP results (only two residues), two anchor points are selected for this loop. Thus, the loop boundary for this loop would be starting from the E previous to the first anchor point, to the E on the left
#### In this case, the anchoring point is the anchor region
def find_pdb_num_of_adjacent_E (iden_code,chainType,pure_pos,dssp_dir):
    # Find the adjacent beta-strand of the residue at the pure_pos of the sequence containing both V&C
    # pure_pos is a tuple,numbering starts from 0
    chainTypeNum=0
    if chainType.lower()=="l":
        chainTypeNum=1
    chainid=iden_code.split("_")[1][chainTypeNum]

    t_dssp=extract_dssp_table(f"{dssp_dir}/{iden_code}.dssp")
    chain_dssp=t_dssp.loc[t_dssp["chain_id"]==chainid] # extract the dssp results of the chain of interest
    chain_dssp=chain_dssp.reset_index(drop=True) # the index of the chain_dssp is the pure_pos of the sequence containing both V and C region


    previous_E_info=chain_dssp.loc[(chain_dssp.index<pure_pos[0])&(chain_dssp["2nd_struc"]=="E")].iloc[-1,:]
    previous_E_pdb_num=previous_E_info["resnum_pdb"].item()

    try:
        next_E_info=chain_dssp.loc[(chain_dssp.index>pure_pos[1])&(chain_dssp["2nd_struc"]=="E")].iloc[0,:]
        next_E_pdb_num=next_E_info["resnum_pdb"].item()
    except:
        # in case at the end terminal of fab structure, there is no beta-strand following the anchor point
        next_E_pdb_num=float("inf")

    return (previous_E_pdb_num+1,next_E_pdb_num-1)

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

def get_loops_pdb_boundary (iden_code,chainType,dssp_dir,df):
    # return the pure_seq_pos of the loops location
    # First do a pairwise alignment to its corresponding ref_seq
    # Then pure_pos(ref_seq) -> aln_pos(ref_seq) = aln_pos(vcab_seq) -> pure_pos(vcab_seq)
    # pure_pos(ref_seq) is pre-calculated and stored in global variable ref_loop_pure_pos (read from .json file)
    sub_df=df.loc[df["iden_code"]==iden_code]

    coor_seq=sub_df['H_coordinate_seq'].item() # Note: the query seq contains both V and C
    c_info=sub_df['Htype'].item()
    c_type=c_info.split("(")[0]
    if chainType.lower()=="l":
        coor_seq=sub_df['L_coordinate_seq'].item()
        c_type=sub_df['LSubtype'].item()

    ref_seq=ref_seq_db[c_type][2]

    # Perform the pairwise alignment
    matrix = matlist.blosum62
    pair_aln=pairwise2.align.globalds(coor_seq, ref_seq, matrix,-10,-0.5,penalize_end_gaps=False)[0]
    # use blosum62 as the matrix, gap penalties follow the default value of EMBOSS needle
    # just take the first one as the aln result
    #return pair_aln

    alned_vcab_seq=pair_aln.seqA
    alned_ref_seq=pair_aln.seqB

    # Convert the pure_pos of ref_loop_anchor_point into the pure_pos of vcab_seq_anchor_point,
    # Then find the adjacent "E"s around this pure_pos
    pdb_pos_lst=[]
    ref_loop_anchor_pure_pos=ref_seq_loop_anchor_site_pure_pos[c_type]

    for (pos_start,pos_end) in ref_loop_anchor_pure_pos:
        # Map the pure_seq_pos to the aln_pos of the query seq
        ref_aln_pos_start=map_pure_seq_pos_to_aln_pos (pos_start,alned_ref_seq)
        seq_aln_pos_start=ref_aln_pos_start # The aln_pos for VCAb seq and ref_seq are equal
        seq_pure_pos_start=map_aln_pos_to_pure_seq_pos (seq_aln_pos_start,alned_vcab_seq)

        ref_aln_pos_end=map_pure_seq_pos_to_aln_pos (pos_end,alned_ref_seq)
        seq_aln_pos_end=ref_aln_pos_end # The aln_pos for VCAb seq and ref_seq are equal
        seq_pure_pos_end=map_aln_pos_to_pure_seq_pos (seq_aln_pos_end,alned_vcab_seq)

        seq_pure_pos=(seq_pure_pos_start,seq_pure_pos_end)

        pdb_pos=find_pdb_num_of_adjacent_E (iden_code,chainType,seq_pure_pos,dssp_dir)
        pdb_pos_lst.append(pdb_pos)

    return pdb_pos_lst

def get_loops (chain,loop_pdb_boundary):
    # Return the loops fragment as a chain object
    # Chain is an Bio.PDB.Chain.Chain object
    # vcb is the pdb_vcb
    # loop_p_pos is a list containing 7 numbers indicating the loops position.
    ## i.e.[loop1_start,loop1_end,loop2_start,loop_2_end,loop_3_start, loop3_end, end_terminal]
    ## the positions above are PURE_SEQ POSITIONS, including the residues at that position
    frag=Bio.PDB.Chain.Chain("n")
    pure_pos_counter=0
    for res in chain:
        res_pdb=res.id[1]
        res_type=res.id[0]

        for i in loop_pdb_boundary[0:-1]: # exclude the final anchoring point for the C terminal
            if res_pdb >=i[0] and res_pdb <=i[1] and res_type==' ':
                frag.add(res)

    return frag

def get_CH1_CL_angle(iden_code,dssp_dir,pdb_dir,df):
    # return the angle between CH1-CL, defined by the center of mass (COM)
    # fn: the file name of the pdb structure
    # h,l: the chain id of H & L chains
    # ch1_r, cl_r, ch1_loop_r, cl_loop_r: the range of the ch1, cl, ch1_loop, cl_loop, in the format of lst/tuples

    # Acquire the structure:
    parser=Bio.PDB.PDBParser()
    structure=parser.get_structure(iden_code,f"{pdb_dir}/{iden_code}.pdb")

    # Acquire fragments in order to acquire COM to calculate the angles
    h,l=iden_code.split("_")[1]

    hchain_obj=structure[0][h]
    lchain_obj=structure[0][l]

    ch1_start,cl_start=find_pdb_vcb_in_HandL(iden_code,df)
    ch1_end,cl_end=3000,3000 # Note: this will cause problem if the antibody is a full ig

    ch1=get_frag (hchain_obj, ch1_start, ch1_end)
    cl=get_frag (lchain_obj, cl_start, cl_end)

    # get the loops fragment:
    # get the boundaries of the loops (in the format of pdb_numbers)
    ch1_loop_pdb_pos=get_loops_pdb_boundary (iden_code,"h",dssp_dir,df)
    ch1_loop=get_loops (hchain_obj,ch1_loop_pdb_pos)

    cl_loop_pdb_pos=get_loops_pdb_boundary (iden_code,"l",dssp_dir,df)
    cl_loop=get_loops (lchain_obj,cl_loop_pdb_pos)


    # Acquire center of mass (com) of these fragments
    ch1_com=ch1.center_of_mass()
    cl_com=cl.center_of_mass()
    ch1_loop_com=ch1_loop.center_of_mass()
    cl_loop_com=cl_loop.center_of_mass()
    #return (ch1_loop_com,cl_loop_com)
    #return (ch1_com,cl_com,ch1_loop_com,cl_loop_com)
    return torsion_angle([cl_loop_com,cl_com,ch1_com,ch1_loop_com])

# PART 2. Elbow angle calculation

def find_V_C_linker_pdb_pos (iden_code,dssp_dir,chainType,df):
    # extract the V_C_linker
    ## (defined as the linker between the last beta-strand in V region and the first beta-strand in C region)
    chainTypeNum=0
    if chainType.lower()=="l":
        chainTypeNum=1
    vc_boundary=find_pdb_vcb_in_HandL(iden_code,df)[chainTypeNum]
    chainid=iden_code.split("_")[1][chainTypeNum]

    t_dssp=extract_dssp_table(f"{dssp_dir}/{iden_code}.dssp")
    previous_E_info=t_dssp.loc[(t_dssp["chain_id"]==chainid)&(t_dssp["resnum_pdb"]<vc_boundary)&(t_dssp["2nd_struc"]=="E")].iloc[-1,:]
    previous_E_pdb_num=previous_E_info["resnum_pdb"].item()

    next_E_info=t_dssp.loc[(t_dssp["chain_id"]==chainid)&(t_dssp["resnum_pdb"]>vc_boundary)&(t_dssp["2nd_struc"]=="E")].iloc[0,:]
    next_E_pdb_num=next_E_info["resnum_pdb"].item()

    return (previous_E_pdb_num+1,next_E_pdb_num-1)

def calc_elbow_angle(iden_code,dssp_dir,pdb_dir,df):
    # Acquire the structure:
    parser=Bio.PDB.PDBParser()
    structure=parser.get_structure(iden_code,f"{pdb_dir}/{iden_code}.pdb")

    # Acquire fragments in order to acquire COM to calculate the angles
    h,l=iden_code.split("_")[1]

    hchain_obj=structure[0][h]
    lchain_obj=structure[0][l]

    h_vcb,l_vcb=find_pdb_vcb_in_HandL(iden_code,df)

    vh=get_frag (hchain_obj, 0, h_vcb-1)
    ch1=get_frag (hchain_obj, h_vcb, 3000) # could lead to problem if it is a full ig

    vl=get_frag(lchain_obj,0,l_vcb-1)
    cl=get_frag (lchain_obj, l_vcb, 3000)

    v_total=Bio.PDB.Chain.Chain("v")
    c_total=Bio.PDB.Chain.Chain("c")

    for res in vh:
        v_total.add(res)
    for res in vl:
        id_lst=list(res.id)
        res.id=(id_lst[0],id_lst[1]+3000,id_lst[2])
        # change the id of residue, because Bio.PDB would not add two residue of the same id into one chain
        v_total.add(res)

    for res in ch1:
        c_total.add(res)
    for res in cl:
        id_lst=list(res.id)
        res.id=(id_lst[0],id_lst[1]+3000,id_lst[2])
        c_total.add(res)

    # get the V-C linker fragment:
    # get the boundaries of the loops (in the format of pdb_numbers)
    h_vc_linker_pos=find_V_C_linker_pdb_pos (iden_code,dssp_dir,"h",df)
    h_vc_linker=get_frag (hchain_obj, h_vc_linker_pos[0], h_vc_linker_pos[1])

    l_vc_linker_pos=find_V_C_linker_pdb_pos (iden_code,dssp_dir,"l",df)
    l_vc_linker=get_frag (lchain_obj, l_vc_linker_pos[0]+3000, l_vc_linker_pos[1]+3000)

    #return v_total,c_total,h_vc_linker,l_vc_linker
    # Acquire center of mass (com) of these fragments
    v_com=v_total.center_of_mass()
    c_com=c_total.center_of_mass()
    h_vc_linker_com=h_vc_linker.center_of_mass()
    l_vc_linker_com=l_vc_linker.center_of_mass()
    return torsion_angle([v_com,h_vc_linker_com,l_vc_linker_com,c_com])

if __name__=="__main__":
    vcab=pd.read_csv("./new_vcab.csv").drop(columns=["Unnamed: 0"])
    os.system("sh ../dssp/run_dssp.sh")

    ref_seq_db=json.load(open("./ref_Seq_info.json","r"))
    ref_seq_loop_anchor_site_pure_pos=json.load(open("./ref_seq_loops_ref_points.json","r"))
    dssp_dir="../dssp/dssp_results"
    pdb_dir="../pdb_struc/chain_pdb"

    elbow_angles=[]
    ch1_cl_angles=[]
    for i in vcab.index:
        iden_code=vcab.loc[i,"iden_code"]
        if vcab.loc[i,'Structural Coverage']=="full antibody":
            e_angle=np.nan
            c_angle=np.nan
        else:
            try:
                e_angle=calc_elbow_angle(iden_code,dssp_dir,pdb_dir,vcab)
            except:
                e_angle=np.nan # the elbow angle can't be calculated because there is a problem to find the beta-strand position
            try:
                c_angle=get_CH1_CL_angle(iden_code,dssp_dir,pdb_dir,vcab)
            except:
                c_angle=np.nan # the elbow angle can't be calculated because there is a problem to find the beta-strand position
        elbow_angles.append(e_angle)
        ch1_cl_angles.append(c_angle)

    vcab["elbow_angle"]=elbow_angles
    vcab["CH1-CL_interface_angle"]=ch1_cl_angles

    vcab.to_csv("final_vcab.csv")
