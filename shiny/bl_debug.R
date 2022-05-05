library (rBLAST)
library (taxonomizr)

# Directories of blast db:
# ref db:
igh_bl <- "../seq_db/ref_db/human_IGH_db/human_IGH.fasta"
# ref_IGH_seq from uniprot
light_bl <- "../seq_db/ref_db/human_light_chain_db/human_light_constant.fasta"
# ref_L_seq from uniprot
all_ref_bl <- "../seq_db/ref_db/all_ref_db/all_ref.fasta"
# ref_db containing all the reference sequences

# Blast db:
VCAb_fseq_bl <- "../seq_db/vcab_db/all_full_seq.fasta"
# all the full sequence from both H and L chains
VCAb_vseq_bl <- "../seq_db/vcab_db/all_v_seq.fasta"
# all the V region sequence from both H and L chains

VCAbH_bl <- "../seq_db/vcab_db/H_seq.fasta"
# all H_full_seq in VCAb
VCAbL_bl <- "../seq_db/vcab_db/L_seq.fasta" 
# all L_full_seq in VCAb
HV_bl <- "../seq_db/vcab_db/HV_seq.fasta"
# all HV seq in VCAb
LV_bl <- "../seq_db/vcab_db/LV_seq.fasta"
# all LV seq in VCAb

generate_blast_result <- function(o_seq,db_dir,suffix=""){
  str_seq <- as.character(o_seq)
  if (str_seq==""){
    return (NULL)
  }
  else {
      db <- blast(db=db_dir,type="blastp")
      aa_seq <- AAStringSet(str_seq) # Convert the string format into String Set
      seq_pred <- predict(db,aa_seq) # seq_pred is the dataframe containing all the blast results.
      
      ## Don't need this step: the blast result is automatically ranked by "Bits", a measurement of how well query&subject seqs are aligned together.
      # Rank seq_pred
      #seq_pred <- seq_pred[order(seq_pred$E,seq_pred$Perc.Ident,decreasing=TRUE),]
      blast_df <- seq_pred[,2:12]
      
      # generate the column of "iden_code"
      split_id <- function(i,sep){
        strsplit(i,sep)[[1]][1]
      }
      colnames(blast_df) <- paste(colnames(blast_df),suffix,sep="")
      blast_df$iden_code <- unlist(lapply(blast_df$SubjectID,split_id,"-"))
      return (blast_df)}
    
    

  
} # return the dataframe of the complete result of the blast
input_seq <- "EVQLVESGAEVKKPGASVKVSCKVSGYTLTELSMHWVRQAPGKGLEWMGGFDPEDGETMYAQKFQGRVTMTEDTSTDTAYMESSLRSEDTAVYYCATSTAVAGTPDLFDYYYGMDVWGQGTTVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPK"
blast_df <- generate_blast_result(input_seq,VCAbH_bl) # return the dataframe of the complete result of the blast
