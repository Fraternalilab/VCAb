library (shiny)
library (DT)
library (rBLAST)
library (taxonomizr)
library (NGLVieweR)
library (tibble)
library (shinyhelper)
library (ggplot2)

####################### DIRECTORIES FOR ALL THE USED FILES #######################
vcab_dir="../vcab_db/final_vcab.csv"
pops_parent_dir <- "../pops/result/"
f_pdb_dir <- "../pdb_struc/full_pdb/"
pdb_parent_dir <- "../pdb_struc/chain_pdb/"


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
# Example: generate files required for the establishment of blast db:
# makeblastdb("~/Desktop/antibody/human_IGH_db/human_IGH.fasta",dbtype="prot")

# Files required to plot the seq_cov plot
dom_info_dir="../seq_db/ref_db/all_alleles/H_chains/h_alleles_domain_info.csv"
total_dom_info=read.csv(dom_info_dir)
total_dom_info=total_dom_info[,!(names(total_dom_info) %in% c("X"))]

h_author_seq_bl_dir="../vcab_db/blast_result/best_h_seq_bl_result.csv"
h_coor_seq_bl_dir="../vcab_db/blast_result/best_h_coordinate_seq_bl_result.csv"
t_h_author_bl=read.csv(h_author_seq_bl_dir)
t_h_coor_bl=read.csv(h_coor_seq_bl_dir)

# Files required to rank the antibody according to interface similarity:
## dm stands for distance matrix: the matrix holding the interface difference index value
dm_dir="../ch1_cl_interface_matrix/dm_of_interface_dist_mtrx.csv"
dm_df=read.csv(dm_dir)

# Read the vcab database (the csv file)
o_vcab=read.csv(file=vcab_dir) # original vcab table

# Read the unusual case csv
all_unusual_cases = read.csv(file = '../vcab_db/unusual_cases/all_unusual_cases.csv')
all_unusual_cases = all_unusual_cases[, 3:ncol(all_unusual_cases)]

# Round some values to two digit in order to make the table looks better:
vcab=o_vcab
colnames(vcab)[which(names(vcab) == "elbow_angle")] <- "elbow_angle0"
vcab$elbow_angle <- unlist(lapply(vcab$elbow_angle0,function(x){round(x,2)}))
vcab$CH1_CL.interface_angle <- unlist(lapply(vcab$CH1.CL_interface_angle,function(x){round(x,2)}))
vcab=vcab[,!(names(vcab) %in% c("X","CH1.CL_interface_angle","elbow_angle0"))] 

#Note: the positions of seqs:33:36
##vcab=vcab[!(vcab$pdb %in% c("2rcj","7bm5")),]
update_date_fn="../vcab_db/VCAb_db_update_date.txt"
update_date=readChar(update_date_fn, file.info(update_date_fn)$size)

release_fn="../release"
release=readChar(release_fn, file.info(release_fn)$size)

####################### Functions #######################
# Generate blast table: return the best blast result (order: highest iden, alignment_length)
generate_blast_result <- function(o_seq,db_dir,suffix=""){
  str_seq <- as.character(o_seq)
  if (str_seq==""){
    return (NULL)
  }
  else {
    withProgress(message="BLASTing...",value=0,{
      db <- blast(db=db_dir,type="blastp")
      aa_seq <- AAStringSet(str_seq) # Convert the string format into String Set
      seq_pred <- predict(db,aa_seq) # seq_pred is the dataframe containing all the blast results.
      
      ## In web, the blast result is automatically ranked by "Bits", a measurement of how well query&subject seqs are aligned together.
      # Rank seq_pred
      seq_pred <- seq_pred[order(seq_pred$Bits,decreasing=TRUE),]
      blast_df <- seq_pred[,2:12]
      rownames(blast_df) <- NULL
      
      # generate the column of "iden_code"
      split_id <- function(i,sep){
        strsplit(i,sep)[[1]][1]
      }
      colnames(blast_df) <- paste(colnames(blast_df),suffix,sep="")
      blast_df$iden_code <- unlist(lapply(blast_df$SubjectID,split_id,"-"))
      return (blast_df)
    })
    
  }
} # return the dataframe of the complete result of the blast

blast_paired_chains <- function (ab_title,hseq,lseq,region){
  # Generate the blast result only (doesn't extract VCAb information) of the paired H&L chain
  
  # region is the region selected by the user, in order to choose the specific bl_db
  # ab_title would be the value of QueryID, when multiple seqs in fasta file is inputted
  # ab_title would be empty in "manual" mode. i.e. only one H/L pair is inputted.
  
  H_bl <- ifelse(region=="v_region", HV_bl, VCAbH_bl)
  L_bl <- ifelse(region=="v_region", LV_bl, VCAbL_bl)
  
  H_df <- generate_blast_result(hseq,H_bl,".H")
  L_df <- generate_blast_result(lseq,L_bl,".L")
  # Merge the result of VH blast with VL blast:
  # if both not empty, merge; else, return the one that is not empty
  t_df <- if (is.null(H_df)) L_df else if (is.null(L_df)) H_df else merge(H_df,L_df,by="iden_code",suffixes=c(".H",".L"))
  
  # Add the avg_ident to the merged table
  t_df <- within(t_df,avg_ident <- if (is.null(H_df)) Perc.Ident.L else if (is.null(L_df)) Perc.Ident.H else (Perc.Ident.H+Perc.Ident.L)/2)
  t_df <- t_df[with(t_df,order(avg_ident,decreasing=TRUE)),] # sort the table based on avg_ident
  
  rownames(t_df) <- NULL # reset the row index of df
  
  if (ab_title==""){
    t_df <- t_df[,!(names(t_df) %in% c("avg_ident"))] # The column "avg_ident" would not be displayed
    return (t_df)
  }
  else{
    # seq_title would be displayed in the bl_df:
    t_df$QueryID <- ab_title
    # Rearrange the columns of t_df and drop the column avg_ident
    # Note: using merge might break the row order defined by avg_ident
    t_df <- t_df[,c("QueryID","iden_code",
                    "SubjectID.H","Perc.Ident.H","Alignment.Length.H","Mismatches.H","Gap.Openings.H","Q.start.H","Q.end.H","S.start.H","S.end.H","E.H","Bits.H",
                    "SubjectID.L","Perc.Ident.L","Alignment.Length.L","Mismatches.L","Gap.Openings.L","Q.start.L","Q.end.L","S.start.L","S.end.L","E.L","Bits.L")]
    
    return (t_df)
  }
  
}

check_uploaded_file <- function(f_path,seq_max=200){
  tryCatch(
    expr={
      set <- readAAStringSet(f_path)
      # Check the length of the fasta file:
      if (length(set)>seq_max){
        # warning window to user
        showModal(modalDialog(title="File input error", "The fasta file can contain 200 sequences max."))
        return (NULL)
      }
      else{
        return (f_path)
      }
    },
    error = function (e){
      # The file os not a valid fasta file
      showModal(modalDialog(title="File input error", "Please input a valid fasta file"))
      return (NULL)
    },
    warning = function (w){
      showModal(modalDialog(title="File input error", "Please input a valid fasta file"))
      return (NULL)
    }
  )
  
  
}

# Allow the user to upload their own fasta files
uploaded_file_blast_unpaired <- function(f_path,region){
  # f_path: uploaded file path; # db_dir: the directory of the database used for blast
  # region: the region of interest selected by the user(input$sele_bl_ab), in order to select which database would be blasted against.
  
  bl_db <- ifelse(region=="v_region", VCAb_vseq_bl, VCAb_fseq_bl)
  set <- readAAStringSet(f_path) # convert the file into a string set
  
  # read each seq one by one, then blast
  result <- list()
  for (i in 1:length(set)){
    ab_title <- names(set)[i]
    ab_seq <- set[[i]] # return an AAString object
    bl_result <- generate_blast_result(ab_seq,bl_db)
    if (is.null(bl_result)){
      # jump to the next for loop if the bl_result is empty
      next
    }
    best_match <- bl_result[1:3,] # display the top three blast result
    best_match$QueryID <- ab_title
    result <- rbind(result,best_match)
  }
  result <- result[,c("QueryID","SubjectID","iden_code","Perc.Ident","Alignment.Length","Mismatches","Gap.Openings","Q.start","Q.end","S.start","S.end","E","Bits")]
  return (result)
}

uploaded_file_blast_paired <- function (f_path,region){
  # if the number of fasta sequence is within 0-200, return the blast result of vh&vl seq, and the list of unpaired_seqs (if any).
  # f_path: should be the file path of fasta files containing V seqs
  set <- readAAStringSet(f_path)
  
  # Note: the format of the title for the sequence should be exactly like this: AbName-HorL
  seq_title <- names(set)
  title_split <- strsplit(names(set),"-")
  ab_name <- unlist(lapply(title_split,function(i){
    i[1]
  }))
  chain_type <- unlist(lapply(title_split,function(i){
    # i[2]: chain name
    if (i[2] %in% c("H","L")) i[2] else "unknown"
  }))
  ab_seq <-unlist(lapply(seq_title,function(i){
    as.character(set[[i]])
  }))
  
  # Build the df containing seq title, ab_name, chain_type, ab_seq
  seq_info_df <- data.frame(seq_title,ab_name,chain_type,ab_seq)
  
  # Store the paired and unpaired seqs as dataframes
  #paired_seqs <- seq_info_df[NULL,]
  unpaired_seqs <- NULL
  bl_result <- NULL # store all the 
  
  # Inspect chains with the same ab_name
  for (ab in unique(ab_name)){
    o_subset <- seq_info_df %>% dplyr::filter(ab_name==ab)
    subset <- o_subset[with(o_subset,order(chain_type)),] # sort the subset according to chain_type (H first, then L)
    HLtypes <- sort(subset$chain_type)
    if (identical(HLtypes,c("H","L"))) {
      # paired seqs: seqs with the same ab_name, and seqs contain only one H and one L.
      #paired_seqs <- rbind(paired_seqs,subset)
      
      # get the seq for H & L chain, respectively.
      H_seq <- subset$ab_seq[subset$chain_type=="H"]
      L_seq <- subset$ab_seq[subset$chain_type=="L"]
      
      
      # blast for both HV & LV to get the match
      bl_df <- blast_paired_chains(ab,H_seq,L_seq,region)
      
      bl_result <- rbind(bl_result,bl_df[1:3,]) # only attach the top 3 entries to bl_result
      
    }
    else{
      unpaired_seqs <- rbind(unpaired_seqs,subset)
    }
    
  }
  
  # return bl_result of paired & df of unpaired_seqs
  return (list("paired_bl"=bl_result,"unpaired"=unpaired_seqs))
}

# Extract the entry with the same iden_code in the vcab table
extract_entry <- function(title,df){
  name <- strsplit(title,"-")[[1]][1]
  return (df %>% dplyr::filter(iden_code==name))
}

# Filter the table by features (attributes)
filter_the_rows <- function (iso_txt,Ltype_txt,struc_cov,exp_method,res_cut,if_antigen,df=vcab){
  # select the rows to be displayed
  # return a set of TRUE/FALSE value
  
  # all the inputs of this function are the input: e.g. iso_txt should be in this format: input$iso_txt
  chain_type_filter <- function(col_name,user_input){
    if (user_input=="All") TRUE else unlist(lapply(df[[col_name]],function(x){strsplit(x,"\\(")[[1]][1]==user_input}))
  }
  select_filter <- function(check_col_name,select_input_value){
    # check_col_name: the name of the column in the table to be checked / the column which serves as the filter standard
    # o_select_input_value: the value inputted by the user, e.g. input$iso_txt
    if(select_input_value=="All") TRUE else df[[check_col_name]] %in% c(select_input_value) # To select rows, to include the "all" option
  }
  
  antigen_filter <- function(){
    if (if_antigen=="Any") TRUE else if (if_antigen=="No") df[["antigen_chain"]]=="" else df[["antigen_chain"]]!=""
  }
  first_filter <- chain_type_filter("Htype",iso_txt)&chain_type_filter("Ltype",Ltype_txt)&select_filter("Structural.Coverage",struc_cov)&select_filter("method",exp_method)&antigen_filter()
  final_filter <- if (is.na(res_cut)) first_filter else first_filter&(df$resolution <= res_cut) # Keep entries with resolution smaller than the threshold, if the res_cut is inputted.
  
  return (final_filter)
}

# From https://community.rstudio.com/t/how-to-use-shiny-action-button-in-datatable-through-shiny-module/39998
shinyInput <- function(FUN, len, id, ns, ...) {
  inputs <- character(len)
  for (i in seq_len(len)) {
    inputs[i] <- as.character(FUN(paste0(id, i), ...))
  }
  inputs
}

# Add "show" button in every entry of the table
addShow <- function(df,ns){
  Actions = shinyInput(actionButton,nrow(df),'button_',
                       label = "Show",
                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4"
                       #onclick = paste0('Shiny.onInputChange(\"' , ns("select_button"), '\", this.id)')
  )
  return (add_column(df, Structure=Actions, .after="iden_code"))
}

# Generate the total_info table containing both ab_info and blast result
generate_total_info <- function(bl_df,ns){
  # Get the ab_info of all the entries in the blast table
  ab_names <- bl_df$iden_code
  ab_info_lst <- lapply(ab_names,extract_entry,df=vcab)
  ab_info_df <- as.data.frame(do.call(rbind,ab_info_lst))
  
  final_ab_info_df <- addShow(ab_info_df,ns)
  
  
  # drop the columns of iden_code in the bl_df
  bl_df <- bl_df[!(names(bl_df) %in% c("iden_code"))]
  
  # Merge the blast_df with ab_info_df: 
  # use cbind() instead of merge to preserve the original blast order, since merge() would be order the df by the "by" col
  #bl_ab_df <- merge(bl_df,ab_info_df,by="iden_code")
  bl_ab_df <- cbind(bl_df,final_ab_info_df)
  return(bl_ab_df)
}

generate_pops_info <- function(pdb_c){
  # type_pdb_c: the name of the antibody in this format: "igg1_k_7c2l_HL"
  # pdb_c: the name of the antibody in this format: "7c2l_HL"
  
  pops_dir <- paste(pops_parent_dir,pdb_c,"_C_deltaSASA_rpopsResidue.txt",sep='')
  pops_file <- read.csv(file=pops_dir, sep=' ')
  new_pops <- pops_file%>% dplyr::filter(D_SASA.A.2>15) #filtered pops file
  new_pops <- new_pops[,c("Chain","ResidNe","ResidNr","D_SASA.A.2")] # only show these four columns
  
  
  Hchain <- substr(strsplit(pdb_c,"_")[[1]][2],1,1)
  Lchain <- substr(strsplit(pdb_c,"_")[[1]][2],2,2)
  
  hpops <- new_pops[new_pops$Chain==Hchain,]
  lpops <- new_pops[new_pops$Chain==Lchain,]
  
  rownames(hpops) <- NULL
  rownames(lpops) <- NULL
  
  return (list("hpops"=hpops,"lpops"=lpops))
}

get_coverage_pos_plot <- function(iden_code,horltype,ref_dom,t_author_bl,t_coor_bl,df=vcab){
  # horltype should be only "Htype" or "Ltype"
  ctype=df[df$iden_code==iden_code,horltype]
  allele_info=strsplit(ctype,"\\(")[[1]][2]
  allele=strsplit(allele_info,",")[[1]][1]
  
  # (1) annotation of domain starts/end points of the reference allele
  dom_positions = ref_dom[ref_dom$q==allele,]
  rownames(dom_positions) <- NULL
  
  ref_start=head(dom_positions,n=1)$a
  ref_end=tail(dom_positions,n=1)$b
  
  # (2) annotation of the coverage of reference / author sequence / coordinate sequence
  author_hit_info=head(t_author_bl[(t_author_bl$iden_code==iden_code)&(t_author_bl$matched_alleles==allele),],n=1)
  coor_hit_info=head(t_coor_bl[(t_coor_bl$iden_code==iden_code)&(t_coor_bl$matched_alleles==allele),],n=1)
  
  author_start=author_hit_info$start_ref
  author_end=author_hit_info$end_ref
  coor_start=coor_hit_info$start_ref
  coor_end=coor_hit_info$end_ref
  
  coverage_pos <- data.frame(
    a = c(ref_start, author_start, coor_start), b = c(ref_end, author_end, coor_end), # start and end
    q = c(allele, 'author', 'coords') # indicate sequence type
  )
  
  # order the sequence type annotation for the plot
  dom_positions$q <- factor(dom_positions$q, 
                            levels = c(allele, 'author', 'coords'),
                            labels = c(paste0('reference\n',allele), 
                                       'author-submitted\nsequence (HC)',
                                       'atomic\ncoordinate (HC)'))
  coverage_pos$q <- factor(coverage_pos$q, 
                           levels = c(allele, 'author', 'coords'),
                           labels = c(paste0('reference\n',allele), 
                                      'author-submitted\nsequence (HC)',
                                      'atomic\ncoordinate (HC)'))
  
  #return (list("dom_pos"=dom_positions,"cov_pos"=coverage_pos))
  # determine the positions of domain labels (text on the reference allele bar)
  dom_positions$label_pos <- apply(dom_positions[, c("a", "b")], MARGIN = 1, mean)
  
  #_________________________________________________
  # the ggplot2 functions
  
  
  prepare_plot <- function(tb)
  {
    # prepare the ggplot canvas
    ggplot(tb, aes(xmin = a, xmax = b, y = q)) + 
      scale_y_discrete(drop = FALSE, name = "") +
      scale_x_continuous(labels = seq(1, ref_end, by = 100), name = "AA position",
                         breaks = seq(1, ref_end, by = 100)) + theme_bw()
  }
  
  draw_coverage <- function(tb, start_column_name, end_column_name, y_column_name)
  {
    # draw the structural coverage (ie the lighter colour bars)
    # start_column_name and end_column_name corresponds to start/end points
    # of each sequence type (ref / coords/ author)
    # y_column_name is the column that indicate ref/coords/author
    geom_rect(data = tb, 
              aes_string(xmin = start_column_name, xmax = end_column_name),
              ymin = as.numeric(tb[, y_column_name]) - 0.1,
              ymax = as.numeric(tb[, y_column_name]) + 0.1,
              fill = "grey80")
  }
  
  draw_dom_position <- function(tb, y_column_name)
  {
    # draw the domain boundaries with dark coloured rectangles
    # y_column_name refers to ref/author/coords (all should be 'ref'! - 
    # just so that ggplot2 knows on which horizontal line to put the rect)
    geom_rect(ymin = as.numeric(tb[, y_column_name]) - 0.1,
              ymax = as.numeric(tb[, y_column_name]) + 0.1) 
  }
  
  # actual plotting
  g <- prepare_plot(dom_positions)
  g <- g + draw_coverage(coverage_pos, start_column_name = "a",
                         end_column_name = "b", y_column_name = "q")
  g <- g + draw_dom_position(dom_positions, "q")
  g <- g + geom_text(aes(label = dom, x = label_pos), colour = "white")
  g <- g + theme(axis.text = element_text(size = 12))
  # 'g' is the ggplot2 object to be rendered/printed:
  return (g)
  
}

get_similar_interface <- function(iden_code,dm_df){
  sub_df_row=dm_df[dm_df$X==iden_code,!(names(dm_df) %in% c("X"))]
  rownames(sub_df_row) <- NULL
  t_sub_df_row=t(sub_df_row) # first extract the rows, then transpose the df into column
  rownames(t_sub_df_row) <- lapply(rownames(t_sub_df_row),function(x){substring(x,2)})
  colnames(t_sub_df_row) <- c("row_extraction")
  
  sub_df_col=data.frame(col_extraction=dm_df[,paste0("X",iden_code)])
  rownames(sub_df_col) <- dm_df$X
  
  sub_df <- cbind(t_sub_df_row,sub_df_col) # combine the distance listed in the row (iden_code) and col(Xiden_code)
  sub_df$o_interface_difference_index <- apply(sub_df,1,max) # The max value of the row_extraction and col_extraction is the interface_diff_index
  sub_df$interface_difference_index <- unlist(lapply(sub_df$o_interface_difference_index,function(x){round(x,2)}))
  similar_interface_df <- head(sub_df[order(sub_df$interface_difference_index),],11)
  similar_interface_df$iden_code <- rownames(similar_interface_df)
  
  similar_interface_df<- similar_interface_df[,c("iden_code","interface_difference_index")]
  rownames(similar_interface_df) <- NULL
  names(similar_interface_df)[2]<-"interface.difference.index"
  return (similar_interface_df)
}


####################### UI #######################
ui <- fluidPage(
  titlePanel(title=div(img(
    src="./VCAb_logo.png",
    width = 270, height = 80
  ),
  "V and C region bearing Antibody Database"),
  windowTitle = "VCAb antibody database"),
  tags$em(paste0("Last updated: ",update_date)),
  navbarPage("",
             tabPanel("Search",
                      fluidRow(
                        column(7,
                               wellPanel(
                                 # show the user input
                                 tabsetPanel(id = "tabs",
                                             tabPanel("PDB",
                                                      textInput("pdb_txt","Enter the pdb ID","7c2l")
                                                      
                                             ),
                                             tabPanel("Features",
                                                      selectInput("iso_txt","Isotype:",choices=c("All","IgA1","IgA2","IgD","IgE","IgG1","IgG2","IgG3","IgG4","IgM")) %>%
                                                        helper(type="inline",title="Isotype", 
                                                               content=c("There are nine isotypes in human, classified by the sequence of C region on H chain.",
                                                                         "Each isotype has different function.")),
                                                      
                                                      selectInput("Ltype_txt","Light chain type:",choices=c("All","kappa","lambda")) %>%
                                                        helper(type="inline",title="Light chain type", 
                                                               content=c("There are two light chain types in human, classified by the sequence of C region on L chain.")),
                                                      selectInput("struc_cov","Structural Coverage:",choices=c("All",sort(unique(vcab$Structural.Coverage)))) %>%
                                                        helper(type="inline", title="Structural Coverage",
                                                               content=c("In VCAb, the structural coverage is classified as Fab and full antibody.",
                                                                         "Full antibody covers both Fab and Fc region")),
                                                      selectInput("if_antigen","If has antigen:",choices=c("Any","Yes","No")) %>%
                                                        helper(type="inline",title="If has antigen",
                                                               content=c("If the pdb file of this entry containing the antigen chain",
                                                                         "Any: include the antibody in the results no matter if it has antigen or not",
                                                                         "Yes: only include the antibody if the pdb file contains the antibody chain",
                                                                         "No: only include the antibody if the pdb file doesn't contain any antibody chain")),
                                                      
                                                      selectInput("exp_method","Experimental Method:",choices=c("All",sort(unique(vcab$method))),multiple=FALSE,selected="All") %>%
                                                        helper(type="inline",title="Experimental Method",
                                                               content=c("The experimental method used to acquire the structure.")),
                                                      numericInput("res_cut","Resolution Threshold:",NULL,min=1,max=5) %>%
                                                        helper(type="inline",title="Resolution Threshold",
                                                               content=c("This is used to acquire structures with resolution below the threshold.",
                                                                         "The threshold can be set to value from 1 to 5.")) 
                                                      # Just empty the input to allow the user to select ab without the limit of resolution.
                                                      
                                             ),
                                             
                                             tabPanel("Sequence",
                                                      tabsetPanel(id="seq_tabs",
                                                                  tabPanel("Search individual sequence",
                                                                           textAreaInput("seq_txt","Enter the amino acid sequence of the chain",width="600px",rows=5,resize="both",
                                                                                         value="EVQLVESGAEVKKPGASVKVSCKVSGYTLTELSMHWVRQAPGKGLEWMGGFDPEDGETMYAQKFQGRVTMTEDTSTDTAYMESSLRSEDTAVYYCATSTAVAGTPDLFDYYYGMDVWGQGTTVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPK"),
                                                                           
                                                                           # Radio buttons instead of checkBox are used because it only allow the user to select one option.
                                                                           radioButtons(inputId="seq_type", label="The type of this chain:", 
                                                                                        choices=c("Heavy Chain" = "Hseq",
                                                                                                  "Light Chain" = "Lseq",
                                                                                                  "Don't know" = "unknown_seq")),
                                                                           checkboxInput(inputId="two_chains",label="Add the other paired H/L chain",FALSE) %>%
                                                                             helper(type="inline",title="Add the other paired H/L chain",
                                                                                    content=c("Sequences of the paired H & L chains can be acquired via single cell sequencing, select this check box to enable the input of the other paired H or L chain.",
                                                                                              "",
                                                                                              "When the paired chains are inputted, the similarity of both H chain and L chain would be taken into consideration, in order to find VCAb entries with similar sequence.",
                                                                                              "Please note: When the other paired chain is added, the chain types selected for these two sequences must be one H chain and one L chain.")),
                                                                           
                                                                           conditionalPanel(
                                                                             condition="input.two_chains==1",
                                                                             textAreaInput("seq_txt_2","Enter the amino acid sequence of the chain",width="600px",rows=5,resize="both",
                                                                                           value="EIVMTQSPLSSPVTLGQPASISCRSSQSLVHSDGNTYLSWLQQRPGQPPRLLIYKISNRFSGVPDRFSGSGAGTDFTLKISRVEAEDVGVYYCTQATQFPYTFGQGTKVDIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"),
                                                                             
                                                                             # Radio buttons instead of checkBox are used because it only allow the user to select one option.
                                                                             radioButtons(inputId="seq_type_2", label="The type of this chain:", 
                                                                                          choices=c("Heavy Chain" = "Hseq",
                                                                                                    "Light Chain" = "Lseq"))
                                                                           )
                                                                           
                                                                           
                                                                  ),
                                                                  tabPanel("Search in batch",
                                                                           fileInput("upload","Upload a fasta file (200 seqs max)"),
                                                                           radioButtons(inputId="up_paired", label=" Sequences inside the uploaded file are: ",
                                                                                        choices=c(
                                                                                          "paired H and L chains" = "paired",
                                                                                          "not-paired chains" = "unpaired"
                                                                                        )) %>%
                                                                             helper(type="inline",title="What is paired H & L chain?",
                                                                                    content=c("Sequences of the paired H & L chains can be acquired via single cell sequencing.",
                                                                                              "",
                                                                                              "In order for your sequences to be picked up as \"paired\", the title of paired H & L sequences in the uploaded file should be in this format: ",
                                                                                              "AntibodyName-HorL, where AntibodyName is the name of the antibody, such as \"7c2l_HL\"; HorL can only be letter \"H\" or \"L\", in order to specify if the sequence belongs to H or L chain.",
                                                                                              "AntibodyName and HorL are connected by a hyphen.",
                                                                                              "",
                                                                                              "Under the paired mode, the similarity of both H chain and L chain would be taken into consideration, in order to find VCAb entries with similar sequence.",
                                                                                              "If the title of the fasta sequence is not in the format specified above, the sequence would not be picked up as paired, and the unpaired sequences would be tabulated below the result table."))
                                                                           
                                                                  )),
                                                      
                                                      # horizontal line
                                                      tags$hr(), 
                                                      radioButtons(inputId="sele_bl_db", label="Select the region of your interest:",
                                                                   choices=c("V region"="v_region",
                                                                             "Full sequence (V & C)"="full_seq")) %>%
                                                        helper(type="inline",title="Selection of the database to BLAST against",
                                                               content = c("This selection would determine the database to be BLAST against.",
                                                                           "",
                                                                           "If the \"V region\" is selected, the sequence would be BLAST against the database containing only sequences of V region, meaning the search would be based on the V region similarity, without the consideration of C region.",
                                                                           "If \"Full sequence (V & C)\" is selected, the search would be based on the sequence similarity of both V and C region. "))
                                             ),
                                             tabPanel("CH1-CL Interface",
                                                      tags$em("Get VCAb antibody entries with similar residue contacts at the CH1-CL interface."),
                                                      br(),br(),
                                                      selectizeInput("pdb_interface","Enter the iden_code",choices=unique(vcab$iden_code),
                                                                     options=list(maxOptions =5,
                                                                                  placeholder = 'Please type iden_code here',
                                                                                  onInitialize = I('function() { this.setValue(""); }'))) %>%
                                                        #selectInput("pdb_interface","Enter the iden_code",choices=unique(vcab$iden_code),selectize = TRUE) %>%
                                                        helper(type="inline",title="iden_code", 
                                                               content=c("In VCAb, each entry has a unique iden_code, in the format of \"PDBID_HL\", where PDBID is the four-character PDB ID of the antibody, HL are the chain ID of the heavy/light chain.",
                                                                         "If you don't know the heavy/light chain ID, just type the PDBID in the searching box. VCAb iden_code with this PDBID will be automatically listed in the option list as you are typing, then you can click on the corresponding VCAb iden_code to select it.",
                                                                         "CH1-CL interface similarity is ranked by a metric which we termed 'interface difference index'. This is based on considering residue contacts between the CH1 and CL domains for each structure, and comparing these contacts between every pair of structures. The smaller the value, the more similar CH1-CL interface it has, with respect to the query iden_code. For detailed explanation, please go to the VCAb Documentation.",
                                                                         "NOTE: Currently the interface similarity search function only support antibodies within VCAb."
                                                               ))
                                                      #textInput("pdb_interface","Enter the iden_code","7c2l_HL")
                                                      
                                             )
                                             
                                 )
                                 
                               )
                        ),
                        
                        column(5,
                               wellPanel(
                                 tabsetPanel(
                                   tabPanel("Structural Viewer",
                                            # show the structure viewer
                                            textOutput("struc_selected_message") %>%
                                            helper(type="inline",title="Structure viewer",
                                                               content=c("Here users can interactively inspect (zoom, hover etc.) the structure selected in the Antibody Information (bottom left) panel. Explanations:",
									 "1. Amino acid numbering follows the numbering scheme in the visualised .pdb file.",
									 "2. Hovering over the structure you will see pop up bubble with information of the residue indicated by your mouse. This is the format of this print-out message:",
									 "atom : [ 'residue three-letter code' ] 'residue number' : 'chain identifier' . 'atom name' ()",
									 "(items indicated in quotation marks will be changed as the mouse moves.)"
                                                               )),
					    NGLVieweROutput("structure"),
                                            div(img(
                                              src="./struc_viewer_legend.png",
                                              width = 560#, height = 80
                                            )),
                                            checkboxInput(inputId = "if_full_view", label = "View all the chains with the same PDB ID"),
                                            conditionalPanel(
                                              condition="input.if_full_view==1",
                                              checkboxGroupInput("full_view_option","Display options:",
                                                                 c("Color antigen chain in this pdb"="color_antigen",
                                                                   "Show other ligand(s) in this pdb"="show_ligand",
                                                                   "Zoom in to the heavy-light chain pair of the selected VCAb entry"="zoom_in_to_HL"))
                                              
                                              
                                              
                                              #)
                                            ),
                                            downloadButton("download_struc",label="Download the displayed structure (.pdb) file")#,
                                            
                                   ),
                                   tabPanel("Sequence Coverage",
                                            plotOutput("seq_cov_plot")
                                   )
                                 )
                                 
                                 
                               )
                        )
                      ),
                      fluidRow(
                        column(7,
                               wellPanel(
                                 # show the antibody information table
                                 actionButton("search","Search"),
                                 textOutput("chain_type_message"),
                                 selectInput("file_col","Select the additional column(s) you want to display",
                                             choices=colnames(vcab)[!colnames(vcab) %in% c('X','iden_code','Htype','Ltype','Structural.Coverage')], 
                                             multiple=TRUE)%>%
                                   helper(type="inline",title="Display additional columns",
                                          content=c("For simplicity, some columns of the antibody table below are hidden",
                                                    "By selecting the column names listed here, you can acquire the information in these columns",
                                                    "Multiple columns can be selected at the same time")),
                                 fluidRow(
                                   column(4,
                                          downloadButton("download_subset",label="Download search results")
                                   ),
                                   column(8,
                                          checkboxInput(inputId = "if_download_pdb",label="Select this to download pdb files of VCAb entries (.zip)")%>%
                                            helper(type="inline",title="Download search results",
                                                   content=c("The download file would be a zip file containing the antibody table listed here, and pdb files if you check the checkbox",
                                                             "pdb files downloaded contaning one Heavy-Light chain pair indicated in the iden_code would be downloaded",
                                                             "It might take some time if the table containing too many VCAb entries and you want to download the pdb files for each one of them"
                                                   ))
                                   )
                                   
                                 ),
                                 
                                 br(),
                                 #checkboxInput(inputId="filter_result",label="Filter the result by features",FALSE),
                                 conditionalPanel(#condition="input.filter_result==1",
                                   condition="input.tabs==\"Sequence\" ",
                                   # The following options are the same as the "Features" tab:
                                   hr(),
                                   
                                   column(5,
                                          strong ("Filter the results by features:"),
                                          br(),br(),
                                          selectInput("flt_iso_txt","Isotype:",choices=c("All","IgA1","IgA2","IgD","IgE","IgG1","IgG2","IgG3","IgG4","IgM")) %>%
                                            helper(type="inline",title="Isotype", 
                                                   content=c("Any one of: IgA1, IgA2, IgD, IgE, IgG1, IgG2, IgG3, IgG4, IgM, or 'All' (i.e. any isotype). There are nine isotypes in human, classified by the sequence of C region on H chain.",
                                                             "Each isotype has different function.")),
                                          
                                          selectInput("flt_Ltype_txt","Light chain type:",choices=c("All","kappa","lambda")) %>%
                                            helper(type="inline",title="Light chain type", 
                                                   content=c("Any one of: kappa, lambda, or All (i.e. either kappa or lambda). There are two light chain types in human, classified by the sequence of C region on L chain.")),
                                          selectInput("flt_struc_cov","Structural Coverage:",choices=c("All",sort(unique(vcab$Structural.Coverage)))) %>%
                                            helper(type="inline", title="Structural Coverage",
                                                   content=c("In VCAb, the structural coverage is classified as Fab and full antibody.",
                                                             "Full antibody covers both Fab and Fc region"))
                                   ),
                                   column(5, offset=2,
                                          
                                          selectInput("flt_if_antigen","If has antigen:",choices=c("Any","Yes","No")) %>%
                                            helper(type="inline",title="If has antigen",
                                                   content=c("If the pdb file of this entry containing the antigen chain",
                                                             "Any: include the antibody no matter if the pdb file contains the antibody chain or not",
                                                             "Yes: only include the antibody if the pdb file contains the antibody chain",
                                                             "No: only include the antibody if the pdb file doesn't contain any antibody chain")),
                                          selectInput("flt_exp_method","Experimental Method:",choices=c("All",sort(unique(vcab$method))),multiple=FALSE,selected="All") %>%
                                            helper(type="inline",title="Experimental Method",
                                                   content=c("The experimental method used to acquire the structure.")),
                                          numericInput("flt_res_cut","Resolution Threshold:",NULL,min=1,max=5) %>%
                                            helper(type="inline",title="Resolution Threshold",
                                                   content=c("This is used to acquire structures with resolution below the threshold.",
                                                             "The threshold can be set to value from 1 to 5.")) 
                                          # Just empty the input to allow the user to select ab without the limit of resolution.
                                   )#,
                                   #actionButton("filter","Filter the Results")
                                   
                                 ),
                                 
                                 hr(),
                                 DT::dataTableOutput("ab_info_table"),
                                 textOutput("unpair_message"),
                                 DT::dataTableOutput("unpaired_table")
                                 
                               )
                        ),
                        column(5,
                               wellPanel(
                                 tabsetPanel(id="residue_list_panel",
                                             tabPanel("CH1-CL interface residues",
                                                      # show the filtered(DSASA <= 15) POPSComp table: show H_pops & L_pops separately
                                                      textOutput("pops_message") %>%
                                                      helper(type="inline",title="CH1-CL interface residues",
                                                               content=c("Interface analysis was performed using POPSComp (Cavallo, Kleinjung and Fraternali NAR (2003), doi: 10.1093/nar/gkg601. (github: https://github.com/Fraternalilab/POPScomp).",
									 "Amino acid numbering follows the numbering scheme in the displayed/analysed .pdb file."
                                                               )),
                                                      br(),
                                                      
                                                      actionButton("clear_sele_res","Clear selected residues"),
                                                      checkboxInput(inputId = "zoom_in_sele_res",label="Zoom in to selected residues",FALSE),
                                                      tabsetPanel(
                                                        tabPanel("H chain residues",
                                                                 DT::dataTableOutput("h_pops")),
                                                        tabPanel("L chain residues",
                                                                 DT::dataTableOutput("l_pops"))
                                                      )
                                                      
                                             ),
                                             tabPanel("Disulfide Bond",
                                                      br(),
                                                      # show the filtered(DSASA <= 15) POPSComp table: show H_pops & L_pops separately
                                                      actionButton("clear_sele_disulfide","Clear selected residues"),
                                                      checkboxInput(inputId = "zoom_in_sele_disulfide",label="Zoom in to selected Cysteines",FALSE),
                                                      DT::dataTableOutput("disulfide_info")
                                             )
                                             
                                 )
                               )
                        )
                      )
             ),
	     tabPanel("About",
		      h5(paste0("Version: ", release)),
		      br(),
                      h5("Documentation can be found in the following link:"),
                      tags$a(href="https://github.com/Fraternalilab/VCAb/wiki", "VCAb github wiki"),
                      h5(),
		      h5("If you have used VCAb in your work, please cite: "),
		      br(),
		      h5("Guo Dongjun, Joseph Chi-Fung Ng, Deborah K Dunn-Walters, Franca Fraternali. VCAb: An accurate and queryable database of isotype annotation for human antibody structures. Under review, 2022"),
		      br(),
   		      div(img(src="./Fig2_VCAb_db.png",
                              width = 1080#, height = 80
			      )),
		      h5("VCAb (V and C region bearing antibody) database is established with the purpose to clarify the annotation of isotype and structural coverage of human antibody structures, and provide an accessible and easily consultable resource. For each antibody entry, users can search for its sequence, isotype, structure and details of the CH1-CL interface. The structure and the CH1-CL interface residues of the antibody can be visualized and inspected in the web server. Users can search the VCAb by entering the PDB identifiers, attributes (e.g. isotype, structural coverage, experimental methods, etc.), single sequence or sequences in batches. Researchers interested in antibody annotations and structures would benefit from the VCAb database, especially due to the curated information it provides on isotype, light chain type and the CH1-CL interface residues. "),
		      br()
	     ),
             tabPanel("Statistics",
                      sidebarLayout(
                        sidebarPanel(
                          checkboxGroupInput(inputId = "stat", label="The statistics of VCAb based on:",
                                             choices=c("Isotypes"="Htype",
                                                       "Light chain types"="Ltype",
                                                       "Structural Coverage"="Structural.Coverage")),
                          uiOutput("total_ab_enties")
                        ),
                        mainPanel(
                          plotOutput("statistics_plt")
                        )
                      )
                      
                      
             ),
             tabPanel("Download",
                      ## Allow the user to download the entire database
                      downloadButton("download",label="Download all entries in VCAb"),
		      br(),
		      downloadButton("download_unusual",label="Download all 'unusual' antibody structures removed in VCAb"),
                      br(),br(),br(),br(),br(),br()
                      
             )
  )
  
)

####################### Server #######################
server <- function(input,output,session){
  
  observe_helpers(withMathJax = TRUE) # allow the display of the helper message when the user click on the question mark next to some text.
  
  ### THE PANEL of antibody information table #################################################
  ns <- session$ns
  
  
  ## Assign different contents to the panel of ab_info_table when different tabs are selected in the input panel
  tabs_value <- reactive({input$tabs}) # the title of the tabs
  seq_sub_tab_value <- reactive({input$seq_tabs})
  ab_info <- reactiveValues() # The reactive value to hold the table of antibody information
  two_chains <- reactive({input$two_chains})
  pdbs <- vcab$pdb
  
  observeEvent(input$search,{
    # Initialize everything:
    initialize_everything()
    # Do the actual job
    if (tabs_value() == "PDB"){
      if(input$pdb_txt %in% pdbs){
        o_df <- vcab %>% dplyr::filter(pdb==input$pdb_txt)
        final_df <- addShow(o_df,ns)
        ab_info$ab_info_df <- final_df
        ab_info$chain_type_message <- NULL
      }
      else{
        ab_info$ab_info_df <- NULL
        ab_info$chain_type_message <- "This pdb is not in the database."
      }
    }
    else if (tabs_value() == "Features"){
      #o_df <- vcab[select_filter("Htype",iso_txt)&select_filter("Ltype",Ltype_txt)&select_filter("Structural.Coverage",struc_cov)&select_filter("method",exp_method),]
      #o_df <- if (is.na(res_cut)) o_df else o_df[o_df$resolution <= res_cut,] # Keep entries with resolution smaller than the threshold, if the res_cut is inputted.
      o_df <- vcab[filter_the_rows(input$iso_txt,input$Ltype_txt,input$struc_cov,input$exp_method,input$res_cut,input$if_antigen),]
      rownames(o_df) <- NULL # reset the index
      
      final_df <- addShow(o_df,ns)
      ab_info$ab_info_df <- final_df
      ab_info$chain_type_message <- NULL
    }
    else if (tabs_value() == "Sequence"){
      
      ### Identify if the user upload the fasta file ###
      if (is.null(input$upload)){
        ## IF the user didn't upload the file, 
        ## Display chain type message, the table containing ab_info & blast_result
        
        # If the user doesn't know if the sequence is H/L chain, identify it. Otherwise, display the type(isotype, L chain type) of the chain:
        chain_type_mess <- function(seq_txt,seq_type){
          # Identify chain type by blasting against the reference sequence
          ref_db_dir <- ifelse(seq_type == "Hseq",igh_bl,ifelse(seq_type == "Lseq",light_bl,all_ref_bl))
          bl_result <- generate_blast_result(seq_txt,ref_db_dir)
          if (is.null(bl_result)){
            return (NULL)
          }
          else{
            type_title <- bl_result[1,1]
            type <- tail(strsplit(type_title,"\\|")[[1]],n=1)
            
            if(seq_type=="unknown_seq"){
              Hnames <- c("IGHG1_HUMAN","IGHG2_HUMAN","IGHG3_HUMAN","IGHG4_HUMAN","IGHM_HUMAN","IGHA1_HUMAN","IGHA2_HUMAN","IGHD_HUMAN","IGHE_HUMAN")
              if (type %in% Hnames){
                return ("This chain is likely to be a H chain. Please select the 'Heavy chain' option under the sequence box for further identification.")
              }
              else{
                return ("This chain is likely to be a L chain. Please select the 'Light chain' option under the sequence box for further identification.")
              }
            }
            else{
              return (paste("The type of this chain matches to ", type,
                            ".","\n","The top ten entries in VCAb which are the most similar to the input sequence are shown below."))
            }
          }
          
        }
        
        ## Check if the user choose to input two chains ##
        if (two_chains()==1){
          # If there are two chains:
          # Check if the chain types inputted by the user contain one H and one L
          if (setequal(c(input$seq_type,input$seq_type_2),c("Hseq","Lseq"))){
            # The chain types inputted by the user is correct, do the job:
            
            # Identify which one is H chain, which one is L chain:
            h_seq <- input$seq_txt
            l_seq <- input$seq_txt_2
            if (input$seq_type=="Lseq"){
              h_seq <- input$seq_txt_2
              l_seq <- input$seq_txt
            }
            
            t_df <- blast_paired_chains("",h_seq,l_seq, input$sele_bl_db) # the total blast table for VH&VL sequence
            bl_df <- t_df
            
            bl_ab_df <- generate_total_info(bl_df,ns)
            new_bl_ab_df <- bl_ab_df[!(names(bl_ab_df) %in% c("avg_ident"))] # drop the column of "avg_ident"
            
            ab_info$ab_info_df <- new_bl_ab_df
            
            
          }
          else{
            # remind the user to input correct chain types
            showModal(modalDialog(title="Wrong input for \"chain type\" ", 
                                  "When sequences for both chains of the antibody is inputted, the chain types selected must be exactly one H and one L chain."))
          }
          
        }
        else{
          # If there is only one chain inputted by the user
          
          # Find out which db should be used for blast:
          blast_v_db <- ifelse(input$seq_type == "Hseq",HV_bl,ifelse(input$seq_type == "Lseq",LV_bl,VCAb_vseq_bl)) # the v_seq db
          #note: when the user select "unknown type" for V region seqs, it would blast against all the V seqs
          
          blast_f_db <- ifelse(input$seq_type == "Hseq",VCAbH_bl,ifelse(input$seq_type == "Lseq",VCAbL_bl,all_ref_bl)) # the full_seq db
          blast_db_dir <- ifelse(input$sele_bl_db=="v_region", blast_v_db, blast_f_db)
          
          # BLAST:
          blast_df <- generate_blast_result(input$seq_txt,blast_db_dir) # return the dataframe of the complete result of the blast
          if (is.null(blast_df)==FALSE){
            VCAb_blast_df <- blast_df[1:10,] # Only show the top ten hits
            
            # To observe if the user choose the correct chain type
            if(input$seq_type != "unknown_seq"){
              top_blast_ident <- VCAb_blast_df[1,"Perc.Ident"]
              if (top_blast_ident<50){
                showModal(modalDialog(title="Are you sure you select the correct chain type?", "The Perc. Ident for the top hits are too low, 
                                      you might want to use the 'Don't know' option to further confirm the chain type of this sequence."))
              }
              }
            
            # Record the order of the blast result
            VCAb_blast_df$blast_order <- 1:nrow(VCAb_blast_df)
            
            bl_ab_df <- generate_total_info(VCAb_blast_df,ns) # Get the total_info table containing both bl&ab info
            
            new_bl_ab_df <- bl_ab_df[order(bl_ab_df$blast_order),]
            rownames(new_bl_ab_df) <- NULL
            new_bl_ab_df <- new_bl_ab_df[!(names(new_bl_ab_df) %in% c("blast_order"))]
            
            ab_info$ab_info_df <- new_bl_ab_df
            #ab_info$ab_info_df <- VCAb_blast_df
            ab_info$chain_type_message <- chain_type_mess(input$seq_txt,input$seq_type)
            }
          
          
          
        }
        
        
      }
      else{
        # IF the user upload the file, 
        # only display the table containing ab_info & blast_result, without the chain type message.
        
        uploaded <- input$upload
        uploaded_path <- check_uploaded_file(uploaded$datapath) # check if the uploaded file is a fasta file containing 200 seqs max.
        
        # Check if the uploaded file is a fasta file containing 200 seqs max.
        if (is.null(uploaded_path)){
          # When error happens, i.e. it's not a fasta file containing max.200 seqs
          ab_info$ab_info_df <- NULL
          ab_info$unpaired <- NULL
        }
        else{
          # When the uploaded file is a valid fasta file containing 200 seqs max.:
          ## Check if the user choose to input paired chains ##
          if (input$up_paired=="paired"){
            # if the uploaded file containing paired H & L chains
            pair_info <- uploaded_file_blast_paired(uploaded_path,input$sele_bl_db)
            unpaired <- pair_info$unpaired # allow the user to download this table later
            paired_bl <- pair_info$paired_bl
            paired_bl_ab <- generate_total_info(paired_bl,ns)
            
            ab_info$ab_info_df <- paired_bl_ab
            ab_info$unpaired <- unpaired
            
          }
          else{
            # if the sequence in uploaded file is unpaired
            bl_result <- uploaded_file_blast_unpaired(uploaded_path,input$sele_bl_db) # blast result
            bl_ab_df <- generate_total_info(bl_result,ns) # Get the total_info table containing both bl&ab info
            
            ab_info$ab_info_df <- bl_ab_df
          }
          
          
        }
        
      }
      
    }
    
    else if (tabs_value()=="CH1-CL Interface"){
      #IF the user wants to search according to the interface similarity
      iden_code <- input$pdb_interface
      if (iden_code==""){
        ab_info$ab_info_df <- NULL
        ab_info$chain_type_message <- "Please enter a valid iden_code"
      }
      else{
        interface_info <- get_similar_interface(iden_code,dm_df)
        total_table <- generate_total_info(interface_info,ns)
        ab_info$ab_info_df <- total_table
      }
      
      
      
    }
    else{
      
      ab_info$ab_info_df <- NULL
      ab_info$chain_type_message <- NULL
    }
    if (is.null(ab_info$ab_info_df)==FALSE){
      for (i in 1:nrow(ab_info$ab_info_df)){
        ab_info$ab_info_df[i,"Structure"] <- sprintf(
          ' %s <a href="http://pdbe.org/%s" target="_blank" onmousedown="event.preventDefault(); event.stopPropagation(); return false;"; >PDBe</a> <a href="http://rcsb.org/structure/%s" target="_blank" onmousedown="event.preventDefault(); event.stopPropagation(); return false;"; >RCSB</a> <a href="http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/structureviewer/?pdb=%s" target="_blank" onmousedown="event.preventDefault(); event.stopPropagation(); return false;"; >SAbDab</a>',
          ab_info$ab_info_df[i,"Structure"], substr(ab_info$ab_info_df[i,"iden_code"], 1, 4), substr(ab_info$ab_info_df[i,"iden_code"], 1, 4) , substr(ab_info$ab_info_df[i,"iden_code"], 1, 4)
        )
      }
    }
    
  })
  
  # Make the ab_info_df downloadable for users
  
  output$download_subset <- downloadHandler(
    filename = function(){
      o_time=Sys.time()
      o1 <- gsub(" ","",o_time)
      o2 <- gsub(":","",o1)
      o3 <- gsub("-","",o2)
      final_time <- o3
      paste0("vcab_user_search_",final_time,".zip")
    },
    content=function(f_name){
      temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
      dir.create(temp_directory)
      
      # The result table:
      search_result <- ab_info$ab_info_df[!(names(ab_info$ab_info_df) %in% c("Structure"))] # exclude the "Structure" Column
      colnames(search_result) <- stringr::str_replace_all(colnames(search_result),"<br>","\\.")
      table_name <- "vcab_user_search_result_table.csv"
      write.csv(search_result,file.path(temp_directory,table_name))
      
      # PDB files:
      if (input$if_download_pdb==1){
        pdb_dir=dir.create(file.path(temp_directory,"chain_pdb"))
        for (i in ab_info$ab_info_df["iden_code"]){
          pdb_name<-paste0(i,".pdb")
          #file.copy(paste0(pdb_parent_dir,i,".pdb"),file.path(temp_directory,pdb_name))
          file.copy(paste0(pdb_parent_dir,i,".pdb"),file.path(temp_directory,"chain_pdb",pdb_name))
        }
      }
      
      zip::zip(zipfile=f_name,files=dir(temp_directory),root=temp_directory)
    },
    contentType="application/zip"
  )
  
  # Replace the "." with <br> to allow the col header change to a new line
  observeEvent(input$search,{
    if (!(is.null(ab_info$ab_info_df))){
      colnames(ab_info$ab_info_df) <- stringr::str_replace_all(colnames(ab_info$ab_info_df),"\\.","<br>")
    }
  })
  
  
  # Find which columns will be displayed
  bl_col <- c("QueryID","SubjectID","Perc<br>Ident","E")
  two_bl_col <- c("QueryID","Perc<br>Ident<br>H","E<br>H","Perc<br>Ident<br>L","E<br>L")
  
  file_col_idx <- reactive({
    if (tabs_value()=="PDB" | tabs_value()=="Features"){
      return (names(ab_info$ab_info_df) %in% c('iden_code','Structure','Htype','Ltype','Structural<br>Coverage',stringr::str_replace_all(input$file_col,"\\.","<br>")))
    }
    else if ((seq_sub_tab_value()=="Search individual sequence" & two_chains()==1)| (seq_sub_tab_value()=="Search in batch" & input$up_paired=="paired")){
      return (names(ab_info$ab_info_df) %in% c(two_bl_col,'iden_code','Structure','Htype','Ltype','Structural<br>Coverage',stringr::str_replace_all(input$file_col,"\\.","<br>")))
      #return (names(ab_info$ab_info_df))
    }
    else if (tabs_value()=="CH1-CL Interface"){
      return (names(ab_info$ab_info_df) %in% c("interface<br>difference<br>index",'iden_code','Structure','Htype','Ltype','Structural<br>Coverage',stringr::str_replace_all(input$file_col,"\\.","<br>")))
    }
    else {
      return (names(ab_info$ab_info_df) %in% c(bl_col,'iden_code','Structure','Htype','Ltype','Structural<br>Coverage',stringr::str_replace_all(input$file_col,"\\.","<br>")))
    }
  })
  
  # Filter the rows based on the user input
  file_row_idx <- function(df){
    reactive({
      if (tabs_value()=="Sequence"){
        return (filter_the_rows(input$flt_iso_txt,input$flt_Ltype_txt,input$flt_struc_cov,input$flt_exp_method,input$flt_res_cut,input$flt_if_antigen,df))
      }
      else{
        return (TRUE)
      }
    })
  }
  
  
  
  
  
  
  
  # Show the information
  output$chain_type_message <- renderText({ab_info$chain_type_message})
  output$ab_info_table <- DT::renderDataTable({
    DT::datatable(ab_info$ab_info_df[file_row_idx(ab_info$ab_info_df)(),file_col_idx()],selection="single",escape=F, 
                  options=list(scrollX = TRUE))
    #DT::datatable(ab_info$ab_info_df,selection="single",escape=F,options=list(scrollX=TRUE)) # for test
  })
  
  # Display the message if there are unpaired seq in the fasta file
  unpair_mess <- reactive({if (is.null(ab_info$unpaired)) "" else "Unpaired sequences identified in the fasta file:"})
  
  output$unpair_message <- renderText({unpair_mess()})
  output$unpaired_table <- DT::renderDataTable({
    DT::datatable(ab_info$unpaired,selection="none",options=list(scrollX = TRUE))
  })
  
  
  
  ### THE PANEL of POPSComp table ##################################################################################################
  
  struc_selected <- eventReactive(input$ab_info_table_rows_selected, {
    # return the "id" of the structure selected by the user in the following format: H(type)L(type)_pdb_HL
    # IMPORTANT: this format might be changed if I decided to use a different format for the PDB file
    selectedRow <- as.numeric(input$ab_info_table_rows_selected)
    Ab_info_df <- ab_info$ab_info_df
    pdb_c <- Ab_info_df[selectedRow,'iden_code']
    return (pdb_c)
  })
  
  
  # Generate the df of popscomp table & the message before the table. Hold them in reactive values, so that these values can be re-initialized easily.
  pops_info <- reactiveValues()
  observe({
    pdb_c <- struc_selected()
    # In case POPSComp results is not available. e.g. the structures solved by solution scattering
    tryCatch(
      expr={
        total_pops <- generate_pops_info(pdb_c)
        pops_info$h_df <- total_pops$hpops
        pops_info$l_df <- total_pops$lpops
        pops_info$mess <- paste("All the residues of",struc_selected(),"identified to be involved in the CH1-CL interface is shown below. \n
                                Please click the rows to label the interface residues of your interest.")
      },
      error = function (e){
        pops_info$h_df <- NULL
        pops_info$l_df <- NULL
        pops_info$mess <- paste("The POPSComp analysis of ",struc_selected()," is not available, 
                                probably because the PDB structure doesn't include the information about the side chains.")
        return (NULL)
      },
      warning = function (w){
        pops_info$h_df <- NULL
        pops_info$l_df <- NULL
        pops_info$mess <- paste("The POPSComp analysis of ",struc_selected()," is not available, 
                                probably because the PDB structure doesn't include the information about the side chains.")
        return (NULL)
      }
    )
    
    
      })
  
  output$pops_message <- renderText({pops_info$mess})
  output$h_pops <- DT::renderDataTable({
    DT::datatable(pops_info$h_df,
                  options = list(orderClasses = TRUE,
                                 pageLength = 10,
                                 scrollX = TRUE,
                                 selection='multiple'))
  })
  output$l_pops <- DT::renderDataTable({
    DT::datatable(pops_info$l_df,
                  options = list(orderClasses = TRUE,
                                 pageLength = 10,
                                 scrollX = TRUE,
                                 selection='multiple'))
  })
  
  # Allow the user to reset the selected interface residues
  h_pops_proxy=dataTableProxy('h_pops')
  l_pops_proxy=dataTableProxy("l_pops")
  observeEvent(input$clear_sele_res,{
    h_pops_proxy %>% selectRows(NULL)
    l_pops_proxy %>% selectRows(NULL)
  })
  
  ###### Disulfide Information ######
  disulfide <- reactiveValues()
  observe({
    pdb_c <- struc_selected()
    disulfide_str <- vcab[vcab$iden_code==pdb_c,"disulfide_bond"]
    
    d_rows=strsplit(disulfide_str,", ")[[1]]
    pairs=c()
    dist=c()
    for (r in d_rows){
      cols=strsplit(r,":")[[1]]
      pairs=c(pairs,cols[1])
      dist=c(dist,cols[2])
    }
    disulfide$df <-data.frame(
      disulfide_bond=pairs,
      distance=dist
    )
    
  })
  output$disulfide_info <- DT::renderDataTable({
    DT::datatable(disulfide$df,
                  options = list(orderClasses = TRUE,
                                 pageLength = 10,
                                 scrollX = TRUE,
                                 selection='multiple')
    )
  })
  
  # Allow the user to reset the selected interface residues
  disulfide_proxy=dataTableProxy('disulfide_info')
  observeEvent(input$clear_sele_disulfide,{
    disulfide_proxy %>% selectRows(NULL)
  })
  
  
  ### THE PANEL of 3D Structural viewer ##################################################################################################
  ###### TABPANEL:Draw the seq_pos figure ######
  output$seq_cov_plot<- renderPlot({
    pdb_c <- struc_selected()
    get_coverage_pos_plot (pdb_c,"Htype",total_dom_info,t_h_author_bl,t_h_coor_bl,df=vcab)
  })
  
  ###### TABPANEL: Display the antibody structure ######
  # Get the directory of the pdb file:
  pdb_dir_val <- reactiveValues()
  which_pdb_dir <- reactive({
    pdb_c <- struc_selected()
    pdbid <- strsplit(pdb_c,"_")[[1]][1]
    
    if (input$if_full_view==1){
      return (paste0(f_pdb_dir,pdbid,".pdb"))
    }
    else{
      return (paste(pdb_parent_dir,pdb_c,".pdb",sep=''))
    }
  })
  observeEvent(input$ab_info_table_rows_selected,{
    type_pdb_c <- struc_selected()
    pdbid <- strsplit(type_pdb_c,"_")[[1]][1]
    #pdb_dir_val$pdb <- pdbid
    pdb_dir_val$iden_code <- type_pdb_c
    
    #pdb_dir_val$dir <- which_pdb_dir()
    pdb_dir_val$mess <- paste("The structure of",struc_selected(),"is shown below:")
  })
  
  output$struc_selected_message <- renderText({pdb_dir_val$mess})
  
  # Select the residues to be labeled
  select_residues_to_be_labeled <- function(rows_selected, pops_table){
    # return a string, containig all the residues to be labeled, in this format "134:H or 135:H or 136:L or 137:L" (res_position:chain_name)
    # row_selected: the selected rows of pops table. e.g. input$h_pops_rows_selected, is a string
    # pops_table: e.g. pops_info$h_df
    selected_pops_row_idx <- as.numeric(rows_selected)
    selected_pops <- pops_table[selected_pops_row_idx,]
    
    position_str_vec <- paste (selected_pops$ResidNr, ":", selected_pops$Chain, sep='')
    position_str <- paste(position_str_vec,collapse=' or ')
    
    return (position_str)
  }
  
  select_disulfide_residues_to_be_labeled <- function(disulfide_rows,disulfide_table){
    # Select the residues involved in disulfide bond
    sele_disulfide_row_idx <- as.numeric(disulfide_rows)
    sele_disulfide <- disulfide_table[sele_disulfide_row_idx,"disulfide_bond"]
    disulfide_pos_vec=c()
    for (pair in sele_disulfide){
      residues=strsplit(pair,"-")[[1]]
      for (r in residues){
        remove_res_name=substring(r,4)
        str1=gsub("\\(",":",remove_res_name)
        res_num_chain=gsub("\\)","",str1)
        disulfide_pos_vec <- c(disulfide_pos_vec,res_num_chain)
      }
    }
    
    disulfide_pos_str <-paste(disulfide_pos_vec,collapse=' or ')
    return (disulfide_pos_str)
  }
  h_res_select <- reactive({
    select_residues_to_be_labeled(input$h_pops_rows_selected,pops_info$h_df)
  })
  
  l_res_select <- reactive({
    select_residues_to_be_labeled(input$l_pops_rows_selected,pops_info$l_df)
  })
  
  disulfide_select <- reactive({
    select_disulfide_residues_to_be_labeled(input$disulfide_info_rows_selected,disulfide$df)
  })
  
  output$test_out <-renderText({
    disulfide_select()
  })
  
  # This repetitive code seems unavoidable:
  res_select <- reactive({
    h_selected_pops_row_idx <- as.numeric(input$h_pops_rows_selected)
    h_pops_table <- pops_info$h_df
    l_selected_pops_row_idx <- as.numeric(input$l_pops_rows_selected)
    l_pops_table <- pops_info$l_df
    
    h_selected_pops <- h_pops_table[h_selected_pops_row_idx,]
    l_selected_pops <- l_pops_table[l_selected_pops_row_idx,]
    selected_pops <- rbind(h_selected_pops,l_selected_pops)
    
    position_str_vec <- paste (selected_pops$ResidNr, ":", selected_pops$Chain, sep='')
    position_str <- paste(position_str_vec,collapse=' or ')
    #paste(select_residues_to_be_labeled(input$h_pops_rows_selected,pops_info$h_df)," or ",select_residues_to_be_labeled(input$l_pops_rows_selected,pops_info$l_df))
    
    # Include the selected disulfide bonds
    sele_disulfide_row_idx <- as.numeric(input$disulfide_info_rows_selected)
    disulfide_table<- disulfide$df
    sele_disulfide <- disulfide_table[sele_disulfide_row_idx,"disulfide_bond"]
    disulfide_pos_vec=c()
    for (pair in sele_disulfide){
      residues=strsplit(pair,"-")[[1]]
      for (r in residues){
        remove_res_name=substring(r,4)
        str1=gsub("\\(",":",remove_res_name)
        res_num_chain=gsub("\\)","",str1)
        disulfide_pos_vec <- c(disulfide_pos_vec,res_num_chain)
      }
    }
    
    disulfide_pos_str <-paste(disulfide_pos_vec,collapse=' or ')
    total_pos_str <- paste0(position_str," or ",disulfide_pos_str)
    #return (total_pos_str)
    return (position_str)
  })
  
  # Allow HC & LC to be colored with different colors, select the residues based on VCBoundary:
  HC_color_select <- reactive({
    pdb_c <- struc_selected()
    hl <- strsplit(pdb_c,"_")[[1]][2]
    hChain <- substr(hl,1,1)
    
    # Found the boundary between V and C:
    h_vc <- o_vcab[o_vcab$iden_code==pdb_c,"pdb_H_VC_Boundary"]
    return(paste(h_vc,"-1000:",hChain,sep=""))
  })
  
  LC_color_select <- reactive({
    pdb_c <- struc_selected()
    hl <- strsplit(pdb_c,"_")[[1]][2]
    lChain <- substr(hl,2,2)
    
    # Found the boundary between V and C:
    l_vc <- o_vcab[o_vcab$iden_code==pdb_c,"pdb_L_VC_Boundary"]
    return(paste(l_vc,"-1000:",lChain,sep=""))
  })
  
  HV_color_select <- reactive({
    pdb_c <- struc_selected()
    hl <- strsplit(pdb_c,"_")[[1]][2]
    hChain <- substr(hl,1,1)
    
    # Found the boundary between V and C:
    h_vc <- o_vcab[o_vcab$iden_code==pdb_c,"pdb_H_VC_Boundary"]
    return(paste("0-",h_vc,":",hChain,sep=""))
  })
  
  LV_color_select <- reactive({
    pdb_c <- struc_selected()
    hl <- strsplit(pdb_c,"_")[[1]][2]
    lChain <- substr(hl,2,2)
    
    # Found the boundary between V and C:
    l_vc <- o_vcab[o_vcab$iden_code==pdb_c,"pdb_L_VC_Boundary"]
    return(paste("0-",l_vc,":",lChain,sep=""))
  })
  
  antigen_color_select <- reactive({
    pdb_c <- struc_selected()
    antigen_info <- vcab[vcab$iden_code==pdb_c,"antigen_chain"]
    if (antigen_info==""){
      return ("")
    }
    else{
      HL <- strsplit(pdb_c,"")[[1]]
      antigen_chain_vec=strsplit(antigen_info," \\| ")[[1]]
      antigen_chain_vec <- antigen_chain_vec[!antigen_chain_vec %in% HL] #The antigen annotation fetched from SAbDab is not always correct. In some cases, the H/L chain ID is in this column
      str1=unlist(lapply(antigen_chain_vec,function(x){paste0(":",x)}))
      
      return (paste(str1,collapse=" or "))
    }
  })
  
  if_full_view_value <- reactive({input$if_full_view})
  full_view_option_value<- reactive({input$full_view_option})
  
  HL_chain_selection <- reactive({
    pdb_c<-struc_selected()
    HL <- strsplit(pdb_c,"_")[[1]][2]
    hl_vec <- strsplit(HL,"")[[1]]
    str1=unlist(lapply(hl_vec,function(x){paste0(":",x)}))
    
    return (paste(str1,collapse=" or "))
  })
  
  zoom_in_sele_res_val <- reactive({input$zoom_in_sele_res})
  zoom_in_sele_disulfide_val <- reactive({input$zoom_in_sele_disulfide})
  
  multi_condition <- function(con1,con2){return (con1&con2)}
  
  output$structure <- renderNGLVieweR({
    # Colors:
    # V region: grey
    # C region: CH1 red, CL blue
    # labeled residues: yellow
    NGLVieweR(which_pdb_dir()) %>%
      #NGLVieweR(pdb_dir_val$dir) %>%
      addRepresentation("cartoon", param = list(name = "cartoon", colorScheme =
                                                  "element",colorValue="gray")) %>% #set the whole struc. grey
      addRepresentation("cartoon", param = list(name = "cartoon", colorScheme =
                                                  "element",colorValue="red",sele=HC_color_select())) %>% #set the HC region red
      addRepresentation("cartoon", param = list(name = "cartoon", colorScheme =
                                                  "element",colorValue="blue",sele=LC_color_select())) %>% #set the LC region blue
      addRepresentation("cartoon", param = list(name = "cartoon", colorScheme =
                                                  "element",colorValue="salmon",sele=HV_color_select())) %>% #set the HV region red
      addRepresentation("cartoon", param = list(name = "cartoon", colorScheme =
                                                  "element",colorValue="cornflowerblue",sele=LV_color_select())) %>% #set the LV region blue
      #addRepresentation("ball+stick", param = list(name = "ball+stick", colorScheme =
      #                                            "element",colorValue="cornflowerblue",sele="hetero")) %>% #show the heteroatoms
      # Conditional pipe:
      
      `if`(h_res_select()!=":", addRepresentation(., "ball+stick", param = list(colorScheme = "element",colorValue = "yellow",sele = h_res_select())),.) %>%
      `if`(l_res_select()!=":", addRepresentation(., "ball+stick", param = list(colorScheme = "element",colorValue = "green",sele = l_res_select())),.) %>%
      `if`(disulfide_select()!="", addRepresentation(., "ball+stick", param = list(colorScheme = "element",colorValue = "orange",sele = disulfide_select())),.) %>%
      
      `if`(res_select()!=":", addRepresentation(.,"label",param = list(sele = res_select(),labelType = "format",labelFormat = "%(resname)s %(resno)s", labelGrouping = "residue",color = "white",fontFamiliy = "sans-serif",xOffset = 1,yOffset = 0,zOffset = 0,fixedSize = TRUE,radiusType = 1,radiusSize = 1.5,showBackground = FALSE)),.) %>%
      `if`(disulfide_select()!="", addRepresentation(.,"label",param = list(sele = disulfide_select(),labelType = "format",labelFormat = "%(resname)s %(resno)s", labelGrouping = "residue",color = "white",fontFamiliy = "sans-serif",xOffset = 1,yOffset = 0,zOffset = 0,fixedSize = TRUE,radiusType = 1,radiusSize = 1.5,showBackground = FALSE)),.) %>%
      
      `if`(multi_condition("color_antigen" %in% full_view_option_value(),antigen_color_select()!=""),addRepresentation(., "cartoon", param = list(name = "cartoon", colorScheme ="element",colorValue="black",sele=antigen_color_select())),.) %>%
      `if`("show_ligand" %in% full_view_option_value(), addRepresentation(., "ball+stick", param = list(name = "ball+stick", colorScheme ="element",colorValue="cyan",sele="ligand")),.) %>%
      
      `if`("zoom_in_to_HL" %in% full_view_option_value(), zoomMove(.,center=HL_chain_selection(),zoom=HL_chain_selection()),.) %>%
      `if`(multi_condition(zoom_in_sele_res_val()==1,res_select()!=":"), zoomMove(.,center=res_select(),zoom=res_select()),.) %>%
      `if`(multi_condition(zoom_in_sele_disulfide_val()==1,disulfide_select()!=""), zoomMove(.,center=disulfide_select(),zoom=disulfide_select()),.) %>%
      
      
      
      
      
      
      stageParameters(backgroundColor = "grey") %>%
      #setSize('20','20') %>%
      setQuality("high") %>%
      setFocus(0) 
    
  })
  
  # Make the PDB file downloadable for users
  output$download_struc <- downloadHandler(
    filename = function(){
      o_time=Sys.time()
      o1 <- gsub(" ","",o_time)
      o2 <- gsub(":","",o1)
      o3 <- gsub("-","",o2)
      final_time <- o3
      paste(pdb_dir_val$iden_code,"_",final_time,".pdb",sep="")
    },
    content=function(file){
      #file.copy(pdb_dir_val$dir,file)
      file.copy(which_pdb_dir(),file)
    }
  )
  
  # Make the pdb link available to user:
  output$pdb_link <-renderUI({
    #tags$a(href=paste("http://pdbe.org/",pdb,sep=""),"PDBe",target="_blank")
    if (is.null(pdb_dir_val$iden_code)){
      link="https://www.rcsb.org/"
    }
    else{
      pdb <- toupper(substr(pdb_dir_val$iden_code,1,4))
      link=paste("https://www.rcsb.org/structure/",pdb,sep="")
    }
    em("Go to",a(href=link,"RCSB PDB entry",target="_blank"), "of the antibody above")
  })
  
  
  ### Initialize everything when new tabs are clicked ############################################################################################
  # When different tabs are clicked, initialize everything in the other three panels
  initialize_everything <- function(){
    # The ab_info panel:
    ab_info$ab_info_df <- NULL
    ab_info$chain_type_message <- NULL
    ab_info$unpaired <- NULL
    
    # Don't update this for now: (don't delete this code!)
    # update input$file_col
    #updateSelectInput(session,"file_col",label="Select the additional column(s) you want to display",
    #                  choices=colnames(vcab)[!colnames(vcab) %in% c('pdb','Hchain','Lchain','iden_code','Htype','Ltype','Structural.Coverage')], 
    #                  #multiple=TRUE,
    #                  selected="")
    
    # The popsComp panel:
    pops_info$mess <- NULL
    pops_info$h_df <- NULL
    pops_info$l_df <- NULL
    
    disulfide$df <- NULL
    
    
    # The 3D structural viewer panel:
    pdb_dir_val$iden_code <- NULL
    #pdb_dir_val$dir <- ""
    pdb_dir_val$mess <- NULL
    updateCheckboxInput(session,"if_full_view",value=0)
    updateCheckboxInput(session,"zoom_in_sele_res",value=0)
    updateCheckboxInput(session,"zoom_in_sele_disulfide",value=0)
    updateCheckboxGroupInput(session,"full_view_option","Display options:",
                             c("Color antigen chain in this pdb"="color_antigen",
                               "Show other ligand(s) in this pdb"="show_ligand",
                               "Zoom in to the heavy-light chain pair of the selected VCAb entry"="zoom_in_to_HL"),selected=NULL)
  }
  
  # When new tab is selected, initialize everything
  observeEvent(input$tabs,{
    initialize_everything()
    
  })
  
  # when new seq tab is selected, initialize everything
  observeEvent(input$seq_tabs,{
    initialize_everything()
    
  })
  
  # when new tab in residue_list_panel is selected, initialize checkbox only
  observeEvent(input$residue_list_panel,{
    updateCheckboxInput(session,"zoom_in_sele_res",value=0)
    updateCheckboxInput(session,"zoom_in_sele_disulfide",value=0)
    
  })
  
  
  ###### Statistics ############################################################################################
  output$total_ab_enties <- renderUI({
    em(paste("In total, there are",as.character(nrow(vcab)),"antibody entries in VCAb."))
  })
  output$statistics_plt <- renderPlot({
    in_value <- input$stat 
    
    if (length(in_value)==1){
      ggplot(vcab,aes_string(in_value))+
        geom_bar() +
        scale_y_log10() + 
        geom_text(stat="count",aes(label=..count..), vjust=-1)
    }
    else if (length(in_value)==2){
      in_value <- sort(in_value)
      ggplot(vcab,aes_string(in_value[1], fill = in_value[2]))+
        geom_bar(stat="count", position = position_dodge2(width = 0.5)) +
        scale_y_log10() + 
        geom_text(stat="count",aes(label=..count..), vjust=-1, position = position_dodge2(width = 0.8)) 
    }
    else if (length(in_value)==3){
      in_value <- sort(in_value)
      ggplot(vcab,aes_string(in_value[1], fill = in_value[2]))+
        geom_bar(stat="count", position = position_dodge2(width = 0.5)) +
        scale_y_log10() +
        facet_wrap(~ get(in_value[3])) + 
        geom_text(stat="count",aes(label=..count..), vjust=-1, position = position_dodge2(width = 0.8))
    }
    
    
  })
  
  ###### Download the VCAb database ############################################################################################
  vcab_download=o_vcab
  output$download <- downloadHandler(
    filename = "VCAb.csv",
    content=function(file){
      write.csv(vcab_download,file)
    }
  )
  all_unusual = all_unusual_cases
  output$download_unusual <- downloadHandler(
    filename = 'all_unusual_cases.csv', 
    content = function(file){
      write.csv(all_unusual, file)
    }
  )
  
}


shinyApp(ui=ui,server=server)
