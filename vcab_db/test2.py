import pandas as pd

def extract_hit_bl_result_for_shiny (bl_df,horltype,new_bl_name,df):
    # extract only the bl_result of the hit chain type (when the hit chain type is known)
    # bl_df: the blast result dataframe to be extracted
    # horltype can only be "Htype" or "Ltype"
    new_bl_lst=[]
    for i in df.index:
        iden_code=df.loc[i,"iden_code"]
        ctype=df.loc[i,horltype]
        allele_info=ctype.split("(")[1]
        if "," in allele_info:
            allele=allele_info.split(",")[0]
        else:
            allele=allele_info.split(":")[0]

        sub_bl_df=pd.DataFrame(bl_df.loc[(bl_df["iden_code"]==iden_code)&(bl_df["matched_alleles"]==allele)].iloc[0,:]).T
        new_bl_lst.append(sub_bl_df)
    new_bl=pd.concat(new_bl_lst)

    new_bl.to_csv(f"{new_bl_name}.csv")

hbl = pd.read_csv('hbl.csv')
ff_vcab = pd.read_csv('new_vcab.csv')
extract_hit_bl_result_for_shiny (hbl,"Htype","./blast_result/best_h_seq_bl_result",ff_vcab)
