# VCAb
This repository contains all the source code needed to generate VCAb database and the shiny application of it.

## VCAb database
To generate the VCAb database, go to the 'vcab_db' directory, run 'generate_db.py'.
### Note 
1. In order to run this code, BLAST must be installed in the command line. Please go to [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) for more information.
2. 'generate_db.py' would download PDB files and generate POPSComp results automatically by running bash script. Please make sure the 'pops.sh', 'deltaResidue.R' (Modified from POPSR, please go to [github](https://github.com/Fraternalilab/POPScomp/tree/master/POPSR) for more information), 'pdb_download.sh' are in the corresponding directories.
3. 'seq_db' hold all the BLAST databases. Two files are under this directory: 'ref_db' contains reference sequences for isotypes and light chain types collected from uniprot and won't be changed. 'vcab_db' includes sequences in VCAb, which would be updated when 'generate_db.py' is run.

## VCAb shiny application
To generate the VCAb shiny app, go to 'shiny' and run 'final_VCAb_shiny.R'.
### Note
In order to run this code, BLAST must be installed in the command line. Please go to [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) for more information.

