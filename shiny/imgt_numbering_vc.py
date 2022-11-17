import pandas as pd
from anarci_vc import number,run_anarci,chain_type_to_class

def number_VorC_seq (title, sequence,VorC,hmmerpath=""):
    """
    VorC can only be "V" or "C"
    """
    scheme= "imgt" if VorC=="V" else "imgt_c"
    database= "ALL" if VorC=="V" else "C_ONLY"
    seq_tuple=[(title,sequence)]
    usr_in,num_info,top_hit,five_hits=run_anarci(seq_tuple, ncpu=1,scheme=scheme, database=database, allow=set(["H","K","L","H_C1","K_CC","L_CC"]),allowed_species=None,hmmerpath=hmmerpath)
    num_info=num_info[0]
    if num_info==None:
        return None
    else:
        num_result=num_info[0][0]
        num_df=pd.DataFrame(num_result,columns=["imgt_numbering","res_code"])
        num_df["number"]=[i[0] for i in num_df["imgt_numbering"].values]
        num_df["insertion"]=[i[1] for i in num_df["imgt_numbering"].values]
        num_df["imgt_numbering"]=["".join([str(i[0]),i[1]]).strip() for i in num_df["imgt_numbering"].values]
        #num_df["VorC"]=VorC
        return num_df

def number_usr_input_seq (title, sequence, region,hmmerpath=""):
    """
    region can be "v_region" or "full_seq"
    """

    if region=="v_region":
        vnum=number_VorC_seq (title, sequence,"V",hmmerpath)
        if vnum is not None:
            vnum["VorC"]="V"
        return vnum
    elif region=="full_seq":
        vnum=number_VorC_seq (title, sequence,"V",hmmerpath)
        cnum=number_VorC_seq (title, sequence,"C",hmmerpath)

        if (vnum is not None) and (cnum is not None):
            # Check if the begining of C numbering overlapped with the end of V numbering
            check=["1H","1G","1F","1E"]
            c_check=cnum.loc[map(lambda x: x in check,cnum["imgt_numbering"]),"res_code"]
            c_check_index=c_check.index
            c_check=[i for i in list(c_check.values) if i !="-"]

            vnum_copy=vnum.loc[vnum["res_code"]!="-"]
            v_check=list(vnum_copy.iloc[(-len(c_check)):,:]["res_code"].values)
            if v_check==c_check:
                cnum=cnum.drop(index=c_check_index)


        res_df=pd.concat([vnum,cnum],keys=["V","C"]).reset_index().rename(columns={"level_0":"VorC"}).drop(columns=["level_1"])


        return res_df
