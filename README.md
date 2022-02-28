# VCAb
This repository contains all the source code needed to generate VCAb database and the shiny application of it.

## VCAb database
To generate the VCAb database, go to the `vcab_db` directory, run `generate_db.py`.
### Note 
1. `generate_db.py` should be run in python 3.
2. In order to run this code, BLAST must be installed in the command line. Please go to [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) for more information.
2. `generate_db.py` would download PDB files and generate POPSComp results automatically by running bash scripts. Please make sure the `pops.sh`, `deltaResidue.R` (Modified from POPSR, please go to [github](https://github.com/Fraternalilab/POPScomp/tree/master/POPSR) for more information), `pdb_download.sh` (Modified from the script provided by [RCSB pdb](https://www.rcsb.org/docs/programmatic-access/batch-downloads-with-shell-script)) are in the corresponding directories.
3. `seq_db` hold all the BLAST databases. Two files are under this directory: `ref_db` contains reference sequences for isotypes and light chain types collected from uniprot and won't be changed. `vcab_db` includes sequences in VCAb, which would be updated when `generate_db.py` is run.

## VCAb shiny application
To generate the VCAb shiny app, go to `shiny` and run `final_VCAb_shiny.R`.
### Note
In order to run this code, BLAST must be installed in the command line. Please go to [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) for more information.

## Packages required 
| file | Package | version |
| ---- | ------- | ------- |
|`generate_db.py`| pandas | 1.2.1|
|`generate_db.py`| numpy | 1.19.2|
|`generate_db.py`| Bio | 1.78|
|`generate_db.py`| argparse | 1.1|
|`final_VCAb_shiny.R`| shiny | 1.6.0|
|`final_VCAb_shiny.R`| DT |0.17|
|`final_VCAb_shiny.R`| rBLAST |0.99.2|
|`final_VCAb_shiny.R`| taxonomizr |0.8.0|
|`final_VCAb_shiny.R`| NGLVieweR |1.3.2|
|`final_VCAb_shiny.R`| tibble |3.1.0|
|`final_VCAb_shiny.R`| shinyhelper |0.3.2|
|`final_VCAb_shiny.R`| ggplot2 |3.3.3|
|`final_VCAb_shiny.R`| dplyr |1.0.7|
