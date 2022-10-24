# Dongjun Guo, Joseph Ng. Aug.2022

import pandas as pd
import json
import os
import urllib.request
from datetime import date

from Bio import SeqIO
from Bio.SeqUtils import seq1
from Bio.Seq import Seq
import Bio.PDB


from anarci_vc import number,run_anarci,chain_type_to_class

############ 1.1 Collect all the sequence/newly added sequence since last update from PDB ###########

def collect_all_seq_info (fn,o_fn=None):
    """
    1. Convert the fasta file (downloaded from PDB: https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz)
    into a dictionary in this format: {pdbid:{chainid:seq}}
    2. Only collect the newly added PDB sequences since last update (if any)

    :args fn:file name of the new fasta file
    :args o_fn: file name of the old fasta file
    """

    # Check if there are the old_pdb_seqs_existed:
    if o_fn:
        o_records=SeqIO.parse(o_fn,"fasta")
        o_records_id=[r.id for r in o_records]

        n_records=list(SeqIO.parse(fn,"fasta"))
        n_records_id=[r.id for r in n_records]

        added_id=set(n_records_id)-set(o_records_id)

        records=[r for r in n_records if r.id in added_id]

    else:
        records=SeqIO.parse(fn,"fasta")


    # Convert the Seq_records into this format:{pdbid:{chainid:seq}}
    result={}
    for record in records:

        title_lst=record.description.split(" ")
        seq=str(record.seq).replace("X","")

        mol_type=title_lst[1]
        if mol_type!="mol:protein":
            # Only include protein sequences, exclude DNA/other molecules
            continue

        pdbc=title_lst[0]
        pdb=pdbc.split("_")[0]
        chainid=pdbc.split("_")[1]

        if pdb not in result.keys():
            result[pdb]={chainid:seq}
        else:
            result[pdb][chainid]=seq

    return result

############ 1.2 Get antibody sequences ###########
# 1.2.1 Get the antibody seqs with both V and C regions
def extract_numbering_from_csv(fn1,fn2):
    """
    fn1,fn2 must belongs to the same region (V or C)
    """
    df1=pd.read_csv(fn1)
    df2=pd.read_csv(fn2)
    result={}

    def extract_info_from_one_csv(df):
        for i in df.index:
            pdb,chain=df.loc[i,"Id"].split("_")
            locus=chain_type_to_class[df.loc[i,"chain_type"]]

            numbering_info=df.loc[i,:].iloc[13:] # numbering information, series
            numbering=[(number,residue) for number,residue in numbering_info.items() if residue !="deleted"]

            if pdb not in result.keys():
                result[pdb]={}

            result[pdb][chain]=(''.join([r[1] for r in numbering]),
                                [r[0] for r in numbering],
                                locus
                               )

    extract_info_from_one_csv(df1)
    extract_info_from_one_csv(df2)
    return result

def restructure_numbering_table (o_df,horl):
    """
    return an new df with new column "pdb", "H"/"L"(the chain id column)
    args: o_df: original df
    args: horl: can only be "H" or "L"
    """
    df=o_df.copy()
    df["pdb"]=[i.split("_")[0] for i in df["Id"]]
    df[horl.upper()]=[i.split("_")[1] for i in df["Id"]]
    return df

def filter_for_vc (o_hvn,o_hcn,o_lvn,o_lcn):
    """
    filter for paired antibodies containing both V and C region

    """
    hvn=o_hvn.iloc[:,0:13]
    hcn=o_hcn.iloc[:,0:8]
    lvn=o_lvn.iloc[:,0:13]
    lcn=o_lcn.iloc[:,0:8]

    hvc=pd.merge(hvn,hcn,on="Id",how="inner",suffixes=("_vh", "_ch")) # H chains containing both V&C
    lvc=pd.merge(lvn,lcn,on="Id",how="inner",suffixes=("_vl", "_cl")) # H chains containing both V&C

    n_hvc=restructure_numbering_table (hvc,"H")
    n_lvc=restructure_numbering_table (lvc,"L")

    return n_hvc,n_lvc

# 1.2.2 Re-organize the stored C-numbering results
# For some kind of reason, when anarci storing results into csv, the column 15A&15B,45A&45B is switched, now we want to switch them back
def df_column_switch(o_df, col_lst):
    # col_lst is in this format[(column1,column2)]
    # Then column 1 and column 2 (together with the column values) would be switched
    df=o_df.copy()
    i = list(df.columns)
    for (column1,column2) in col_lst:
        if column1 in df.columns and column2 in df.columns:
            a, b = i.index(column1), i.index(column2)
            i[b], i[a] = i[a], i[b]
    df = df[i]
    return df

############ 1.3 Get paired antibody sequences ###########

## 1.3.1 Download mmcif files so that it can calculate the distance between conserved Cys and determine the paired H&L chains
def generate_pdbid_list (hvc,lvc,out_dir="."):
    """
    Only get the pdbid of antibodies containing both V and C regions

    """
    pdbid=list(set(hvc["pdb"].values)&set(lvc["pdb"].values)) # The intersection of pdbids in both hvc and lvc
    with open(f'{out_dir}/ab_vc_pdbid_lst.txt','w') as output:
        output.write(','.join(pdbid))
    return pdbid

## 1.3.2 get paired sequences
def getSeqFromMMCIF(pdbid,struc_dir):
    """
    get coordinate sequence from each chain in the given cif_file

    :args pdbid: 4-letter pdbid
    :args struc_dir: directory to store the .cif files
    :returns : a dictionary with key = chain ID,
               value = tuple: (one-letter aa code coordinate sequence for the chain,
               list of pdb numbering per residue)
    """

    out = {}

    parser = Bio.PDB.MMCIFParser()
    structure = parser.get_structure(pdbid, f"{struc_dir}/{pdbid}.cif")
    for model in structure:
        for chain in model:
            seq=seq1(''.join(residue.resname for residue in chain)).replace("X","")

            pdb_numbering = []
            for residue in chain.get_residues():
                res_id = residue.get_full_id()[3]
                if res_id[0] == ' ':
                    pdb_numbering.append(''.join([str(i) for i in res_id[1:]]).strip())
            out[chain.get_id()] = (Seq('').join(list(seq)), pdb_numbering)
    return out

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

def new_pairHandL(pdbid, seq_dict,struc_dir):
    """
    pair H and L chains together by inspecting distance between CYS nearest to L104 and H104
    (IMGT numbering) cysteines residue. If <= 22 Angstrom, consider as paired
    (ref SAbDab, Dunbar et al NAR)
    :args pdbid:4-letter pdbid
    :args seq_dict: a dictionary containing only 1 PDB, returned from numberWithANARCI, in this format: {chainid:(seqs,IMGT-numbering,locus(H/L))}
    :args struc_dir: directory to store the cif files
    :returns: a tuple of two elements:
              - a dict object as output from numberWithANARCI
              - list of tuple, each a pair of chain IDs (Hchain, Lchain)
                indicating H-L pairing
    """

    seqs = getSeqFromMMCIF(pdbid,struc_dir)
    imgt_seqs = seq_dict

    parser = Bio.PDB.MMCIFParser()
    structure = f'{struc_dir}/{pdbid}.cif'
    structure = parser.get_structure(pdbid, structure)

    # find the position of the cysteines
    cysteines = {}
    for chain, item in imgt_seqs.items():
        gapped_seq, imgt_numbering, locus = item
        cys = [i for i, n in enumerate(imgt_numbering) if n == '104'][0] # This position should be cys in theory. The string position corresponds to IMGT-104
        all_cys_imgt_pos=[c for c,i in enumerate(str(gapped_seq)) if i=="C"]

        real_cys=cys
        if all_cys_imgt_pos!=[]:
            cys_rp_diff=[abs(i-cys) for i in all_cys_imgt_pos] # The distance between the real cys and IMGT-104
            min_cys_rp_diff=min(cys_rp_diff) # The min distance between the therotical CYS and real CYS
            real_cys=all_cys_imgt_pos[cys_rp_diff.index(min_cys_rp_diff)] # real_cys str_position

        # map aln_pos (imgt_pos with gaps) to pure_seq_pos
        real_cys_pure=map_aln_pos_to_pure_seq_pos (real_cys,gapped_seq)
        #print (chain,seqs[chain][1],real_cys,real_cys_pure)

        cys_pdb_pos = seqs[chain][1][real_cys_pure] # the pdb_numbering position corresponding to H/L 104 [IMGT]
        cysteines[chain] = {'CYS': cys_pdb_pos, 'locus': locus}

    # fetch the Bio.PDB.Residue objects for each selected cysteine
    cysteines_obj = {} # list of residue object in the structure
    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            if chain_id in cysteines.keys():
                for residue in chain:
                    res_id = ''.join([str(i) for i in residue.get_full_id()[3][1:]]).strip()
                    if res_id == cysteines[chain_id]['CYS']:
                        cysteines_obj[chain_id] = residue
    #return cysteines_obj
    # calculate pairwise Ca distance
    pair = []
    cys_chains=list(cysteines_obj.keys())
    for c,chain1 in enumerate(cys_chains):
        for chain2 in cys_chains[c+1:]:
            if chain1 != chain2:
                dist = cysteines_obj[chain2]['CA'] - cysteines_obj[chain1]['CA']
                if dist <= 22:
                    if cysteines[chain1]['locus'] == 'H':
                        hchain, lchain = (chain1, chain2)
                    else:
                        hchain, lchain = (chain2, chain1)
                    pair.append( (hchain, lchain) )

    coor_seq={}
    if pair!=[]:
        all_chains=list(sum(pair,()))
        coor_seq={chain:seqs[chain] for chain in all_chains}

    return coor_seq,list(set(pair))


def get_all_pHL_seqs(pdb_lst,seq_dict,vnum_dict,hvc,lvc,struc_dir):
    """
    :args pdb_lst: list of pdbids of antibodies containing both V&C, H&L
    :args seq_dict: original seq_dict with seqs not numbered
    :args vnum_dict: the reformatted v_numbering_result from run_anarci
    :args hvc,lvc: the csv output containing hmm information of the H/L seqs containing both V and C
    :args struc_dir: the directory of the structure files
    """

    result={}

    errors=[]

    for pdb in pdb_lst:
        numbered=vnum_dict[pdb]
        try:
        #if True:
            coor_seqs,pairs=new_pairHandL(pdb, numbered,struc_dir)

            if pairs==[]: # Only includes paired H&L seqs
                continue

            for pair in pairs:
                h=pair[0]
                l=pair[1]

                vh_num_seq=numbered[h][0]
                vh_imgt_num=",".join(numbered[h][1])
                vl_num_seq=numbered[l][0]
                vl_imgt_num=",".join(numbered[l][1])

                id_code=f"{pdb}_{h}{l}"
                result[id_code]=[pdb,h,l,
                            hvc.loc[hvc["Id"]==f"{pdb}_{h}","identity_species"].item(),
                            hvc.loc[hvc["Id"]==f"{pdb}_{h}","v_gene"].item(),
                            lvc.loc[lvc["Id"]==f"{pdb}_{l}","identity_species"].item(),
                            lvc.loc[lvc["Id"]==f"{pdb}_{l}","v_gene"].item(),
                            str(seq_dict[pdb][h]),
                            str(seq_dict[pdb][l]),
                            str(coor_seqs[h][0]),
                            str(coor_seqs[l][0]),
                            f"{vh_num_seq}:{vh_imgt_num}",
                            f"{vl_num_seq}:{vl_imgt_num}",
                            ",".join(coor_seqs[h][1]),
                            ",".join(coor_seqs[l][1])
                            ]

                assert len(result[id_code])==15, f"column length is shorter:{result[id_code]}"

        except:
            """
            the error happens because it can not return values for things like hvc.loc[hvc["Id"]==f"{pdb}_{h}","identity_species"].item()
            This happens because one pdb have multiple H-L pairs,
            for some chains, they are in the dataframe hvc,
            for some chains, they are not.
            But because they all have the same pdbid, so they are all in the pdb_lst of ab_vc_pdb.
            These entries can simply be excluded.
            """
            errors.append(pdb)


    col_names=["pdb","Hchain","Lchain","vh_species_hmm","heavy_vfamily","vl_species_hmm","light_vfamily",
                "H_seq","L_seq","H_coordinate_seq","L_coordinate_seq",
                "VHseq_IMGT_numbering","VLseq_IMGT_numbering","H_PDB_numbering","L_PDB_numbering"
                ]
    result_df=pd.DataFrame(result).T
    result_df.columns=col_names
    return result_df,errors


############################ Apply the functions ############################
if __name__=="__main__":

    print ("start")
    today=date.today()
    update_date=today.strftime("%d%m%Y")
    downloaded_pdb_seqs=f"pdb_seqres_{update_date}.txt" # The file name (with directory) of the downloaded pdb seqs
    protein_pdb_seqs=f"pdb_protein_seqs_updated_{update_date}.json" # The output file name (with directory) of the filtered protein pdb seqs(contain only the newly added one)

    if os.path.exists("./num_result")==False:
        os.mkdir("num_result")
    if os.path.exists("./blast_result")==False:
        os.mkdir("blast_result")
    if os.path.exists("./result")==False:
        os.mkdir("./result/")
        os.system("mv *.csv result/")
        os.system("cp -r num_result result/")
        os.system("cp -r blast_result result/")
    if os.path.exists("../seq_db/vcab_db/fasta")==False:
        os.mkdir("../seq_db/vcab_db/fasta")
        os.system("mv ../seq_db/vcab_db/*.fasta ../seq_db/vcab_db/fasta")

    # Download PDB seqs:
    urllib.request.urlretrieve('https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz',f"pdb_seqres_{update_date}.txt.gz")
    os.system(f"gzip -d pdb_seqres_{update_date}.txt.gz")

    # Check if the old pdb sequence file existed
    pdb_seq_fns=[i for i in os.listdir(".") if "pdb_seqres_" in i]
    old_pdb_seq_fn=[i for i in pdb_seq_fns if i !=f"pdb_seqres_{update_date}.txt"]
    old_pdb_seq_fn=old_pdb_seq_fn[0] if old_pdb_seq_fn!=[] else None


    struc_dir="../pdb_struc"

    # The prefix of the output file name (with directory) of the numbering csv file
    vnum_out=f"num_result/vnumbering"
    o_cnum_out=f"num_result/o_cnumbering"
    cnum_out=f"num_result/cnumbering"


    pdb_seq_dict=collect_all_seq_info (downloaded_pdb_seqs,old_pdb_seq_fn)
    with open(protein_pdb_seqs, "w") as outfile:
        json.dump(pdb_seq_dict, outfile)

    print ("numbering V&C")
    # Transfer the pdb_seq_dict into the list of tuples:[(seqid,seq)], where seqid is "pdbid_chainid"
    pdb_seq_tuple_lst=[(pdbid+"_"+chainid,str(seq)) for pdbid, info in pdb_seq_dict.items() for chainid,seq in info.items()]
    v_result=run_anarci( pdb_seq_tuple_lst, ncpu=1,scheme="imgt", database="ALL", allow=set(["H","K","L"]),assign_germline=True,allowed_species=None,output=True,csv=True,outfile=vnum_out)
    c_result=run_anarci( pdb_seq_tuple_lst, ncpu=1,scheme="imgt_c", database="C_ONLY", allow=set(["H","K","L","H_C1","K_CC","L_CC"]),allowed_species=None,output=True,csv=True,outfile=o_cnum_out)


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

    # Get the pdbid/cif files of the ab structures containinf V&C H&L
    print ("Downloading PDB structures...")
    # For test:
    #hvc=pd.read_csv("hvc.csv").drop(columns=["Unnamed: 0"])
    #lvc=pd.read_csv("lvc.csv").drop(columns=["Unnamed: 0"])
    #pdb_seq_dict=json.load(open("pdb_protein_seqs20092022.json","r"))

    ab_vc_pdb=generate_pdbid_list (hvc,lvc,out_dir=struc_dir)
    os.system(f"sh {struc_dir}/pdb_download.sh -f {struc_dir}/ab_vc_pdbid_lst.txt -o {struc_dir}/full_pdb/ -c")

    # Reformat the v_result
    n_vresult=extract_numbering_from_csv(f"{vnum_out}_H.csv",f"{vnum_out}_KL.csv")
    pHL,__=get_all_pHL_seqs(ab_vc_pdb,pdb_seq_dict,n_vresult,hvc,lvc,f"{struc_dir}/full_pdb")
    pHL.to_csv("paired_ab.csv")
