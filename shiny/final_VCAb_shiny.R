library (shiny)
library (DT)
library (rBLAST)
library (taxonomizr)
library (NGLVieweR)
library (tibble)

####################### DIRECTORIES FOR ALL THE USED FILES #######################
vcab_dir="new_vcab.csv"
pops_parent_dir <- "./pops/result/"
pdb_parent_dir <- "./pdb_struc/chain_pdb/"

# Directories of blast db:
# ref db:
igh_bl <- "./seq_db/ref_db/human_IGH_db/human_IGH.fasta"
# ref_IGH_seq from uniprot
light_bl <- "./seq_db/ref_db/human_light_chain_db/human_light_constant.fasta"
# ref_L_seq from uniprot
all_ref_bl <- "./seq_db/ref_db/all_ref_db/all_ref.fasta"
# ref_db containing all the reference sequences

# Blast db:
#VCAb_fseq_bl <- "~/Desktop/antibody/VCAb_full_seq_db/all_full_seq_HL.fasta"
# all the full sequence from both H and L chains
VCAbH_bl <- "./seq_db/vcab_db/H_full_db/H_full_seq.fasta"
# all H_full_seq in VCAb
VCAbL_bl <- "./seq_db/vcab_db/L_full_db/L_full_seq.fasta" 
# all L_full_seq in VCAb
HV_bl <- "./seq_db/vcab_db/HV_full_db/HV_full_seq.fasta"
# all HV seq in VCAb
LV_bl <- "./seq_db/vcab_db/LV_full_db/LV_full_seq.fasta"
# all LV seq in VCAb
# Example: generate files required for the establishment of blast db:
# makeblastdb("~/Desktop/antibody/human_IGH_db/human_IGH.fasta",dbtype="prot")


# Read the vcab database (the csv file)
o_vcab=read.csv(file=vcab_dir) # original vcab table
vcab=o_vcab #[,c(32,2:4,18,19,41,37:39,6:16,33:36)]
#Note: the positions of seqs:33:36
vcab=vcab[!(vcab$pdb %in% c("2rcj","7bm5")),]



# Generate blast table: return the best blast result (order: highest iden, alignment_length)
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
    return (blast_df)
  }
} # return the dataframe of the complete result of the blast

blast_vh_vl <- function (ab_title,vhseq,vlseq){
  # ab_title would be the value of QueryID
  H_df <- generate_blast_result(vhseq,HV_bl,".H")
  L_df <- generate_blast_result(vlseq,LV_bl,".L")
  # Merge the result of VH blast with VL blast:
  # if both not empty, merge; else, return the one that is not empty
  t_df <- if (is.null(H_df)) L_df else if (is.null(L_df)) H_df else merge(H_df,L_df,by="iden_code",suffixes=c(".H",".L"))
  
  t_df <- within(t_df,avg_ident <- if (is.null(H_df)) Perc.Ident.L else if (is.null(L_df)) Perc.Ident.H else (Perc.Ident.H+Perc.Ident.L)/2)
  t_df <- t_df[with(t_df,order(avg_ident,decreasing=TRUE)),]
  
  rownames(t_df) <- NULL # reset the row index of df
  
  if (ab_title==""){
    t_df <- t_df[,!(names(t_df) %in% c("avg_ident"))]
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

# Extract the entry with the same iden_code in the best-fit blast result
extract_entry <- function(title,df){
  name <- strsplit(title,"-")[[1]][1]
  return (df %>% dplyr::filter(iden_code==name))
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
                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                       onclick = paste0('Shiny.onInputChange(\"' , ns("select_button"), '\", this.id)'))
  #return (cbind(Structure=Actions,df))
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
read_fasta_file_and_blast <- function(f_path){
  # f_path: uploaded file path; # db_dir: the directory of the database used for blast
  set <- readAAStringSet(f_path) # convert the file into a string set
  # read each seq one by one, then blast
  result <- list()
  for (i in 1:length(set)){
    ab_title <- names(set)[i]
    ab_seq <- set[[i]] # return an AAString object
    bl_result <- generate_blast_result(ab_seq,VCAb_fseq_bl) 
    best_match <- bl_result[1:3,] # display the top three blast result
    best_match$QueryID <- ab_title
    result <- rbind(result,best_match)
  }
  result <- result[,c("QueryID","SubjectID","iden_code","Perc.Ident","Alignment.Length","Mismatches","Gap.Openings","Q.start","Q.end","S.start","S.end","E","Bits")]
  return (result)
}


find_paired_V_seqs_and_blast <- function (f_path){
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
      bl_df <- blast_vh_vl(ab,H_seq,L_seq) # Probably can be downloaded by user?
      
      bl_result <- rbind(bl_result,bl_df[1:3,]) # only attach the top 3 entries to bl_result
      
    }
    else{
      unpaired_seqs <- rbind(unpaired_seqs,subset)
    }
    
  }
  
  # return bl_result of paired & df of unpaired_seqs
  return (list("paired_bl"=bl_result,"unpaired"=unpaired_seqs))
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

ui <- fluidPage(
  titlePanel("VCAb Antibody Structural Database"),
  
  fluidRow(
    column(7,
           wellPanel(
             # show the user input
             tabsetPanel(id = "tabs",
                         tabPanel("PDB",
                                  textInput("pdb_txt","Enter the pdb ID","7c2l")
                                  
                         ),
                         tabPanel("Search by Attributes",
                                  selectInput("iso_txt","Isotype:",choices=c("All",sort(unique(vcab$Htype)))),
                                  selectInput("Ltype_txt","Light chain type:",choices=c("All",sort(unique(vcab$Ltype)))),
                                  selectInput("struc_cov","Structural Coverage:",choices=c("All",sort(unique(vcab$Structural.Coverage)))),
                                  selectInput("exp_method","Experimental Method:",choices=c("All",sort(unique(vcab$method))),multiple=FALSE,selected="All"),
                                  numericInput("res_cut","Resolution Threshold:",NULL,min=1,max=5) # how to allow the user to select all the ab?
                                  
                         ),
                         
                         tabPanel("Sequence",
                                  
                                  textAreaInput("seq_txt","Enter the amino acid sequence of the chain",width="600px",rows=5,resize="both"),
                                  
                                  # Radio buttons instead of checkBox are used because it only allow the user to select one option.
                                  radioButtons(inputId="seq_type", label="The type of this chain", 
                                               choices=c("Heavy Chain" = "Hseq",
                                                         "Light Chain" = "Lseq",
                                                         "Don't know" = "unknown_seq")),
                                  checkboxInput(inputId="two_chains",label="Add another chain",FALSE),
                                  conditionalPanel(
                                    condition="input.two_chains==1",
                                    textAreaInput("seq_txt_2","Enter the amino acid sequence of the chain",width="600px",rows=5,resize="both"),
                                    
                                    # Radio buttons instead of checkBox are used because it only allow the user to select one option.
                                    radioButtons(inputId="seq_type_2", label="The type of this chain", 
                                                 choices=c("Heavy Chain" = "Hseq",
                                                           "Light Chain" = "Lseq"))
                                  ),
                                  fileInput("upload","Or upload a fasta file (200 seqs max)"),
                                  radioButtons(inputId="sele_bl_db", label="Select the region of your interest",
                                               choices=c("V region"="v_region",
                                                         "Full sequence"="full_seq"))
                         ),
                         tabPanel("V_seq+isotype",
                                  textAreaInput("HV_seq_txt","Enter the sequence of the V region of H chain",width="600px",rows=3,resize="both"),
                                  textAreaInput("LV_seq_txt","Enter the sequence of the V region of L chain",width="600px",rows=3,resize="both"),
                                  fileInput("upload_v","Or upload a fasta file containing both VH & VL seqs (200 seqs max)"),
                                  selectInput("iso_select2","Select the name of the isotype you want to search",
                                              choices=c("All Isotypes",sort(unique(vcab$Htype)))
                                  )
                         ),
                         tabPanel("All Entries",
                                  ## Allow the user to download the entire database
                                  downloadButton("download",label="Download")
                                  
                         )
                         
             )
             
           )
    ),
    
    column(5,
           wellPanel(
             textOutput("struc_selected_message"),
             # show the structure viewer
             NGLVieweROutput("structure")
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
                         choices=colnames(vcab)[!colnames(vcab) %in% c('pdb','Hchain','Lchain','iden_code','Htype','Ltype','Structural.Coverage')], 
                         multiple=TRUE),
             
             downloadButton("download_subset",label="Download the searching result"),
             DT::dataTableOutput("ab_info_table"),
             textOutput("unpair_message"),
             DT::dataTableOutput("unpaired_table")
             
           )
    ),
    column(5,
           wellPanel(
             textOutput("pops_message"),
             # show the filtered(DSASA <= 15) POPSComp table: show H_pops & L_pops separately
             tabsetPanel(
               tabPanel("H chain residues",
                        DT::dataTableOutput("h_pops")),
               tabPanel("L chain residues",
                        DT::dataTableOutput("l_pops"))
             )
           )
    )
  ),
  # Add the link to SAbDab
  fluidRow(column(12,
                  "The link to ",
                  tags$a(href="http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/",
                         "SAbDab",
                         target="_blank")))
  
)

server <- function(input,output,session){
  
  ### THE PANEL of antibody information table #################################################
  ns <- session$ns
  
  
  ## Assign different contents to the panel of ab_info_table when different tabs are selected in the input panel
  tabs_value <- reactive({input$tabs}) # the title of the tabs
  ab_info <- reactiveValues() # The reactive value to hold the table of antibody information
  pdbs <- vcab$pdb
  
  observeEvent(input$search,{
    # Initialize everything:
    #ab_info$ab_info_df <- NULL
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
        ab_info$chain_type_message <- "This pdb is not in the database."
      }
    }
    else if (tabs_value() == "Search by Attributes"){
      select_filter <- function(check_col_name,select_input_value,df=vcab){
        if(select_input_value=="All") TRUE else df[[check_col_name]] %in% c(select_input_value)
      }
      
      o_df <- vcab[select_filter("Htype",input$iso_txt)&select_filter("Ltype",input$Ltype_txt)&select_filter("Structural.Coverage",input$struc_cov)&select_filter("method",input$exp_method),]
      o_df <- if (is.na(input$res_cut)) o_df else o_df[o_df$resolution <= input$res_cut,]
      rownames(o_df) <- NULL
      
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
        chain_type_mess <- function(){
          # Identify chain type
          ref_db_dir <- ifelse(input$seq_type == "Hseq",igh_bl,ifelse(input$seq_type == "Lseq",light_bl,all_ref_bl))
          type_title <- generate_blast_result(input$seq_txt,ref_db_dir)[1,1]
          type <- tail(strsplit(type_title,"\\|")[[1]],n=1)
          #type <- chain_type_name()
          if(input$seq_type=="unknown_seq"){
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
        
        blast_db_dir <- ifelse(input$seq_type == "Hseq",VCAbH_bl,ifelse(input$seq_type == "Lseq",VCAbL_bl,all_ref_bl))
        blast_df <- generate_blast_result(input$seq_txt,blast_db_dir) # return the dataframe of the complete result of the blast
        
        
        VCAb_blast_df <- blast_df[1:10,] # Only show the top ten hits
        VCAb_blast_df <- VCAb_blast_df[,c("iden_code","Perc.Ident","E")]
        
        # To observe if the user choose the correct chain type
        observeEvent(input$search,{
          if(input$seq_type != "unknown_seq"){
            top_blast_ident <- VCAb_blast_df[1,"Perc.Ident"]
            if (top_blast_ident<50){
              showModal(modalDialog(title="Are you sure you select the correct chain type?", "The Perc. Ident for the top hits are too low, 
                                    you might want to use the 'Don't know' option to further confirm the chain type of this sequence."))
            }
            }
          })
        
        # Record the order of the blast result
        VCAb_blast_df$blast_order <- 1:nrow(VCAb_blast_df)
        
        bl_ab_df <- generate_total_info(VCAb_blast_df,ns) # Get the total_info table containing both bl&ab info
        ##### IMPORTANT: Avoid extract one ab_entry for multiple times
        
        new_bl_ab_df <- bl_ab_df[order(bl_ab_df$blast_order),]
        rownames(new_bl_ab_df) <- NULL
        new_bl_ab_df <- new_bl_ab_df[!(names(new_bl_ab_df) %in% c("blast_order"))]
        
        ab_info$ab_info_df <- new_bl_ab_df
        #ab_info$ab_info_df <- VCAb_blast_df
        ab_info$chain_type_message <- chain_type_mess()
        }
      else{
        ## IF the user upload the file, 
        ## only display the table containing ab_info & blast_result, without the chain type message.
        uploaded <- input$upload
        uploaded_path <- check_uploaded_file(uploaded$datapath) # check if the uploaded file is a fasta file containing 200 seqs max.
        if (is.null(uploaded_path)){
          ab_info$ab_info_df <- NULL
        }
        else{
          bl_result <- read_fasta_file_and_blast(uploaded_path) # blast result
          bl_ab_df <- generate_total_info(bl_result,ns) # Get the total_info table containing both bl&ab info
          ab_info$ab_info_df <- bl_ab_df
        }
      }
      
    }
    else if (tabs_value() == "V_seq+isotype"){
      ## Generate the merged table containing both blast & Ab info, and this table is filtered with the isotype inputted by user
      
      if (is.null(input$upload_v)){
        # User only input one pair of VH & VL seqs
        t_df <- blast_vh_vl("",input$HV_seq_txt,input$LV_seq_txt) # the total blast table for VH&VL sequence
        
        bl_ab_df <- generate_total_info(t_df,ns)
        new_bl_ab_df <- bl_ab_df[!(names(bl_ab_df) %in% c("avg_ident"))] # drop the column of "avg_ident"
        
        # Filter the bl_ab_df with the isotype inputted by user (if all option, don't filter)
        new_bl_ab_df <- if(input$iso_select2=="All Isotypes") new_bl_ab_df else new_bl_ab_df %>% dplyr::filter(Htype==input$iso_select2)
        
        #final_bl_ab_df <- addShow(new_bl_ab_df)
        
        ab_info$ab_info_df <- new_bl_ab_df
        
        
      }
      else{
        # User input a fasta file containing VH & VL seqs
        up_v <- input$upload_v
        up_v_path <- check_uploaded_file(up_v$datapath) # check if the uploaded file is a fasta file containing 200 seqs max.
        if (is.null(up_v_path)){
          ab_info$ab_info_df <- NULL
          ab_info$unpaired <- NULL
          #print ("hhhhhh")
        }
        else{
          pair_info <- find_paired_V_seqs_and_blast(up_v_path)
          unpaired <- pair_info$unpaired # allow the user to download this table later
          paired_bl <- pair_info$paired_bl
          paired_bl_ab <- generate_total_info(paired_bl,ns)
          
          # Filter the table by the isotype selected by the user:
          paired_bl_ab <- if(input$iso_select2=="All Isotypes") paired_bl_ab else paired_bl_ab %>% dplyr::filter(Htype==input$iso_select2)
          
          #final_paired_bl_ab <- addShow(paired_bl_ab)
          ab_info$ab_info_df <- paired_bl_ab
          ab_info$unpaired <- unpaired
        }
      }
      
      ab_info$chain_type_message <- NULL
      
      
    }
    else{
      # Download all the VCAb entries
      ab_info$ab_info_df <- NULL
      ab_info$chain_type_message <- NULL
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
      paste("vcab_user_search_",final_time,".csv",sep="")
    },
    content=function(file){
      search_result <- ab_info$ab_info_df[!(names(ab_info$ab_info_df) %in% c("Structure"))] # exclude the "Structure" Column
      colnames(search_result) <- stringr::str_replace_all(colnames(search_result),"<br>","\\.")
      write.csv(search_result,file)
    }
  )
  
  # Replace the "." with <br> to allow the col header change to a new line
  observeEvent(input$search,{
    if (!(is.null(ab_info$ab_info_df))){
      colnames(ab_info$ab_info_df) <- stringr::str_replace_all(colnames(ab_info$ab_info_df),"\\.","<br>")
    }
  })
  
  
  # Find which columns will be displayed
  bl_col <- c("QueryID","SubjectID","Perc<br>Ident","E")
  v_bl_col <- c("QueryID","Perc<br>Ident<br>H","E<br>H","Perc<br>Ident<br>L","E<br>L")
  
  file_col_idx <- reactive({
    if (tabs_value()=="PDB" | tabs_value()=="Search by Attributes"){
      return (names(ab_info$ab_info_df) %in% c('iden_code','Structure','Htype','Ltype','Structural<br>Coverage',input$file_col))
    }
    else if (tabs_value()=="Sequence"){
      return (names(ab_info$ab_info_df) %in% c(bl_col,'iden_code','Structure','Htype','Ltype','Structural<br>Coverage',input$file_col))
      #return (names(ab_info$ab_info_df))
    }
    else {
      return (names(ab_info$ab_info_df) %in% c(v_bl_col,'iden_code','Structure','Htype','Ltype','Structural<br>Coverage',input$file_col))
    }
  })
  
  
  # Show the information
  output$chain_type_message <- renderText({ab_info$chain_type_message})
  output$ab_info_table <- DT::renderDataTable({
    DT::datatable(ab_info$ab_info_df[,file_col_idx()],selection="single",escape=F, 
                  options=list(scrollX = TRUE))
    #DT::datatable(ab_info$ab_info_df,selection="none",escape=F,options=list(scrollX=TRUE)) # for test
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
  
  
  ### THE PANEL of 3D Structural viewer ##################################################################################################
  
  # Directory of the pdb file:
  pdb_dir_val <- reactiveValues()
  observeEvent(input$ab_info_table_rows_selected,{
    type_pdb_c <- struc_selected()
    pdb_dir_val$dir <- paste(pdb_parent_dir,type_pdb_c,".pdb",sep='')
    pdb_dir_val$mess <- paste("The structure of",struc_selected(),"is shown below:")
  })
  
  
  #output$struc_selected_message <- renderText({paste("The structure of",struc_selected_mess(),"is shown below:")})
  output$struc_selected_message <- renderText({pdb_dir_val$mess})
  
  # Select the residues to be labeled
  h_res_select <- reactive({
    selected_pops_row_idx <- as.numeric(input$h_pops_rows_selected)
    pops_table <- pops_info$h_df
    
    selected_pops <- pops_table[selected_pops_row_idx,]
    
    position_str_vec <- paste (selected_pops$ResidNr, ":", selected_pops$Chain, sep='')
    position_str <- paste(position_str_vec,collapse=' or ')
    return (position_str)
  })
  
  l_res_select <- reactive({
    selected_pops_row_idx <- as.numeric(input$l_pops_rows_selected)
    pops_table <- pops_info$l_df
    
    selected_pops <- pops_table[selected_pops_row_idx,]
    
    position_str_vec <- paste (selected_pops$ResidNr, ":", selected_pops$Chain, sep='')
    position_str <- paste(position_str_vec,collapse=' or ')
    return (position_str)
  })
  
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
  })
  
  # Allow HC & LC to be colored with different colors:
  H_color_select <- reactive({
    pdb_c <- struc_selected()
    hl <- strsplit(pdb_c,"_")[[1]][2]
    hChain <- substr(hl,1,1)
    
    # Found the boundary between V and C:
    h_vc <- o_vcab[o_vcab$iden_code==pdb_c,"pdb_H_VC_Boundary"]
    return(paste(h_vc,"-1000:",hChain,sep=""))
  })
  
  L_color_select <- reactive({
    pdb_c <- struc_selected()
    hl <- strsplit(pdb_c,"_")[[1]][2]
    lChain <- substr(hl,2,2)
    
    # Found the boundary between V and C:
    l_vc <- o_vcab[o_vcab$iden_code==pdb_c,"pdb_L_VC_Boundary"]
    return(paste(l_vc,"-1000:",lChain,sep=""))
  })
  
  output$structure <- renderNGLVieweR({
    # Colors:
    # V region: grey
    # C region: CH1 red, CL blue
    # labeled residues: yellow
    NGLVieweR(pdb_dir_val$dir) %>%
      addRepresentation("cartoon", param = list(name = "cartoon", colorScheme =
                                                  "element",colorValue="gray")) %>% #set the whole struc. grey
      addRepresentation("cartoon", param = list(name = "cartoon", colorScheme =
                                                  "element",colorValue="red",sele=H_color_select())) %>% #set the HC region red
      addRepresentation("cartoon", param = list(name = "cartoon", colorScheme =
                                                  "element",colorValue="blue",sele=L_color_select())) %>% #set the LC region blue
      # Conditional pipe:
      `if`(h_res_select()!=":", addRepresentation(., "ball+stick", param = list(colorScheme = "element",colorValue = "yellow",sele = h_res_select())),.) %>%
      `if`(l_res_select()!=":", addRepresentation(., "ball+stick", param = list(colorScheme = "element",colorValue = "green",sele = l_res_select())),.) %>%
      `if`(res_select()!=":", addRepresentation(.,"label",param = list(sele = res_select(),labelType = "format",labelFormat = "[%(resname)s]%(resno)s", labelGrouping = "residue",color = "black",fontFamiliy = "sans-serif",xOffset = 1,yOffset = 0,zOffset = 0,fixedSize = TRUE,radiusType = 1,radiusSize = 1.5,showBackground = FALSE)),.) %>%
      
      stageParameters(backgroundColor = "grey") %>%
      #setSize('20','20') %>%
      setQuality("high") %>%
      setFocus(0) 
    
  })
  
  
  ### Initialize everything when new tabs are clicked ############################################################################################
  # When different tabs are clicked, initialize everything in the other three panels
  initialize_everything <- function(){
    # The ab_info panel:
    ab_info$ab_info_df <- NULL
    ab_info$chain_type_message <- NULL
    ab_info$unpaired <- NULL
    # update input$file_col
    updateSelectInput(session,"file_col",label="Select the additional column(s) you want to display",
                      choices=colnames(vcab)[!colnames(vcab) %in% c('pdb','Hchain','Lchain','iden_code','Htype','Ltype','Structural.Coverage')], 
                      #multiple=TRUE,
                      selected="")
    
    # The popsComp panel:
    pops_info$mess <- NULL
    pops_info$h_df <- NULL
    pops_info$l_df <- NULL
    
    
    # The 3D structural viewer panel:
    pdb_dir_val$dir <- ""
    pdb_dir_val$mess <- NULL
  }
  observeEvent(input$tabs,{
    initialize_everything()
    
  })
  
  
  ###### Download the VCAb database ############################################################################################
  vcab_download=o_vcab
  output$download <- downloadHandler(
    filename = "VCAb.csv",
    content=function(file){
      write.csv(vcab_download,file)
    }
  )
  
  
  
  
  }


shinyApp(ui=ui,server=server)