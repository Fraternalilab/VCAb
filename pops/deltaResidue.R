#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#===============================================================================
# POPSR package
# popscomp.R: Implementation of the POPSCOMP functionality,
# i.e. processing of complex structures to compute SASA difference values.
# Returns a list of POPS output files for single-chain and pair-chain structures
#   plus a list of buried SASA values.
#
# (C) 2019-2020 Jens Kleinjung and Franca Fraternali
#===============================================================================
library("stringr")

library("bio3d");
## for atom selection mechanisms see:
## http://thegrantlab.org/bio3d/tutorials/structure-analysis

#library("parallel");

#_______________________________________________________________________________
## POPScomp function implemented in R
## The folling prefixes are used to label the sections and output files
##   for clarity (DIFF output files are called 'delta' for historic reasons):
## ID: the default '--popsr' prefix of POPS for the unmodified input PDB
##      (computed by 'input$popscomp' function in 'app.R')
## ISO: POPS on isolated chains
## PAIR: POPS on paired chains
## DIFF: difference between sum of isolated chain SASA and paired chain SASA
newpopscompR = function(inputPDB_file_name, inDir, outDir) {
  ## number of cores
  #nCore = detectCores() - 1;
  
  #________________________________________________________________________________
  ## ISO: split input PDB into chains
  inputPDB=tail(str_split(inputPDB_file_name,'/')[[1]],n=1) # get the name of the pdb file without the parent directory. e.g. 3m8o_HL_C.pdb
  inputPDB_HLs=str_split(inputPDB, "_")[[1]]
  pdb=inputPDB_HLs[1]
  hl=inputPDB_HLs[2]
  Hid=paste(pdb,substr(hl, 1, 1),sep='_')
  Lid=paste(pdb,substr(hl, 2, 2),sep='_')
  
  chain.files = pdbsplit(paste(inDir, inputPDB, sep = "/"),  ids=c(Hid,Lid),path = outDir, multi = FALSE);
  
  ## if input PDB is not a complex, return without any computations
  ##   bacause this routine "popscompR" is intended for processing protein complexes
  if (length(chain.files) < 2) {
    print("Single-chain POPScomp")
    return(0);
  }
  
  ## short names
  sink(file = "/tmp/basename1")
  print(as.character(chain.files))
  print("")
  print(basename(as.character(chain.files)))
  sink()
  
  chain.files.short = sub('\\.pdb$', '', basename(as.character(chain.files)));
  
  #________________________________________________________________________________
  ## ISO: run POPS over all single (= isolated) chains via system (= shell) call
  sapply(1:length(chain.files), function(x) {
    command = paste0("pops --outDirName ", outDir,
                     " --rout --routPrefix ", paste0(chain.files.short[x], ".iso"),
                     " --residueOut",
                     " --pdb ", chain.files[x], " 1> ", outDir, "/", chain.files.short[x], ".o",
                     " 2> ", outDir, "/", chain.files.short[x], ".e");
    system_status = system(command, wait = TRUE);
    paste("Exit code:", system_status);
  });
  
  ## Concatenate output files of single (ISO = isolated) chains.
  ## We do that here because there is only one tab on the interface for each resolution level
  ##   and the number of chains to be processed/shown will vary. Otherwise we would need
  ##   a dynamic tab structure on the interface that creates a tab for each chain.
  
  ## Residue
  command2.1 = paste0("head -n 1 ", outDir, "/*.iso.rpopsResidue > ", outDir, "/isoSASA.rpopsResidue");
  system_status2.1 = system(command2.1, wait = TRUE);
  command2.2 = paste0("tail -q -n+2 ", outDir, "/*.iso.rpopsResidue >> ", outDir, "/isoSASA.rpopsResidue");
  system_status2.2 = system(command2.2, wait = TRUE);
  
  
  #________________________________________________________________________________
  ## PAIR: create PDB files for all pairwise chain combinations
  pair.cmbn = combn(length(chain.files), 2);
  chainpair.files = vector();
  chainpair.files = sapply(1:dim(pair.cmbn)[2], function(x) {
    ## name of paired chain PDB file to create
    chainpair.files[[x]] = paste0(outDir, "/",
                                  chain.files.short[pair.cmbn[1, x]], "-",
                                  chain.files.short[pair.cmbn[2, x]], ".pdb");
    ## concatenate single chain PDB files to paired chain PDB files
    command = paste("cat", chain.files[pair.cmbn[1, x]],
                    chain.files[pair.cmbn[2, x]], ">",
                    chainpair.files[[x]]);
    system_status = system(command, wait = TRUE);
    paste("Chain pair:", x, "; Exit code:", system_status);
    return(chainpair.files[[x]]);
  });
  
  sink("/tmp/basename2")
  print(as.character(chainpair.files));
  print("")
  print(basename(as.character(chainpair.files)));
  sink()
  chainpair.files.short = sub('\\.pdb$', '', basename(as.character(chainpair.files)));
  
  #________________________________________________________________________________
  ## PAIR: run POPS over all pairwise chain combinations via system (= shell) call
  sapply(1:length(chainpair.files), function(x) {
    command = paste0("pops --outDirName ", outDir,
                     " --rout --routPrefix ", paste0(chainpair.files.short[x], ".pair"),
                     " --residueOut",
                     " --pdb ", chainpair.files[x], " 1> ", outDir, "/POPScomp_chainpair", x, ".o",
                     " 2> ", outDir, "/POPScomp_chainpair", x, ".e");
    system_status = system(command, wait = TRUE);
    paste("Exit code:", system_status);
  });
  
  #________________________________________________________________________________
  ## read SASA files
  ## the data structure will be a list (levels = 'rpopsLevel') of lists (structures)
  rpopsLevel = c("rpopsResidue");
  
  ## ISO: initialise list of lists with predefined number of output files
  iso.sasa.level.files = vector(mode = "list", length = length(rpopsLevel));
  iso.veclist = function(x) { vector(mode = "list", length = length(chain.files)) };
  iso.sasa.level.files = lapply(iso.sasa.level.files, iso.veclist);
  
  ## read ISO SASA files
  for (j in 1:length(rpopsLevel)) {
    for (i in 1:length(chain.files)) {
      ## read isolated chain output
      iso.sasa.level.files[[j]][[i]] = read.table(paste0(outDir, "/", chain.files.short[i],
                                                         ".iso.", rpopsLevel[j]),
                                                  header = TRUE, stringsAsFactors = FALSE);
    };
    names(iso.sasa.level.files[[j]]) = chain.files.short;
  };
  names(iso.sasa.level.files) = rpopsLevel;
  
  ## PAIR: initialise list of lists with predefined number of output files
  pair.sasa.level.files = vector(mode = "list", length = length(rpopsLevel));
  pair.veclist = function(x) { vector(mode = "list", length = dim(pair.cmbn)[2]) };
  pair.sasa.level.files = lapply(pair.sasa.level.files, pair.veclist);
  
  ## read PAIR SASA files
  for (j in 1:length(rpopsLevel)) {
    for (i in 1:dim(pair.cmbn)[2]) {
      ## read paired chain output
      pair.sasa.level.files[[j]][[i]] = read.table(paste0(outDir, "/", chainpair.files.short[i],
                                                          ".pair.", rpopsLevel[j]),
                                                   header = TRUE, stringsAsFactors = FALSE);
    }
    names(pair.sasa.level.files[[j]]) = chainpair.files.short;
  }
  names(pair.sasa.level.files) = rpopsLevel;
  
  #________________________________________________________________________________
  ## DIFF: compute SASA differences (POPScomp values)
  ## 'pair.cmbn' contains the order of PAIR files as column order and
  ##   the index of ISO files as column elements. That way the match between
  ##   PAIR and ISO files is reconstructed here.
  ## initialise list of lists with predefined number of SASA difference tables
  diff.sasa.level = vector(mode = "list", length = length(rpopsLevel));
  diff.veclist = function(x) { vector(mode = "list", length = dim(pair.cmbn)[2]) };
  diff.sasa.level = lapply(diff.sasa.level, diff.veclist);
  
  ## compute SASA differences
  for (j in 1:length(rpopsLevel)) {
    for (i in 1:dim(pair.cmbn)[2]) {
      ## not for level 4 (=molecule)
      if (j %in% 1:3) {
        ## rbind ISO chain SASAs
        iso.rbind.tmp = rbind(iso.sasa.level.files[[j]][[pair.cmbn[1, i]]],
                              iso.sasa.level.files[[j]][[pair.cmbn[2, i]]]);
        ## assert consistency between 'rbind' ISO files and PAIR file
        stopifnot(dim(iso.rbind.tmp) == dim(pair.sasa.level.files[[j]][[i]]));
        #print(paste(j, i, dim(iso.rbind.tmp), dim(pair.sasa.level.files[[j]][[i]])));
        ## SASA DIFF values, applies to all levels
        D_SASA.A.2 = round(iso.rbind.tmp[ , "SASA.A.2"] - pair.sasa.level.files[[j]][[i]][ , "SASA.A.2"], 2);
        ## more level-specific delta values
        if (j == 1) {
          D_Phob.A.2 = round(iso.rbind.tmp[ , "Phob.A.2"] - pair.sasa.level.files[[j]][[i]][ , "Phob.A.2"], digits = 2);
          D_Phil.A.2 = round(iso.rbind.tmp[ , "Phil.A.2"] - pair.sasa.level.files[[j]][[i]][ , "Phil.A.2"], digits = 2);
          diff.tmp.df = cbind(iso.rbind.tmp, D_Phob.A.2, D_Phil.A.2, D_SASA.A.2);
          diff.sasa.level[[j]][[i]] = diff.tmp.df[diff.tmp.df[ , "D_SASA.A.2"] > 0,
                                                  c("ResidNe", "Chain", "ResidNr", "iCode", "D_Phob.A.2", "D_Phil.A.2", "D_SASA.A.2")];
        } 
    };
    names(diff.sasa.level[[j]]) = chainpair.files.short;
  };
  names(diff.sasa.level) = rpopsLevel;
  
  #________________________________________________________________________________
  ## write DIFF SASA result files
  ## ID SASA files have been created in the App
  ## PAIR SASA files have been created here earlier
  ## that completes the set of three types of output files
  for (j in 1:length(rpopsLevel)) {
    write.table(do.call(rbind, diff.sasa.level[[j]]), paste0(outDir, "/", str_split(inputPDB, ".pdb")[[1]][1],"_deltaSASA_", rpopsLevel[j],".txt"));
  }
}
}
#===============================================================================
inputPDB_file_name=args[1]
out_dir=args[3]
inputPDB=tail(str_split(inputPDB_file_name,'/')[[1]],n=1) # get the name of the pdb file without the parent directory. e.g. 3m8o_HL_C.pdb
pops_result_file_name=paste0(out_dir, "/", str_split(inputPDB, ".pdb")[[1]][1],"_deltaSASA_rpopsResidue.txt")

if (file.exists(pops_result_file_name)==FALSE){
  newpopscompR(inputPDB_file_name,inDir=args[2],outDir = out_dir)
  system("rm ./result/!(*deltaSASA_rpopsResidue.txt)") # only keep the *deltaSASA_rpopsResidue.txt file
}

  
