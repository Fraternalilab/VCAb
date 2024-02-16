# VCAb
This repository contains all the source code needed to generate the web interface and assembling the available experimental antibody structure space (Antibodies containing both V and C regions).

**branch `nextflow`**: Branch for wrapping all the code to generate VCAb in a nextflow pipeline. For now the pipeline is in file `vcab.nf`. All python code it imports (i.e. aside from the ones in the `.nf` file itself) needs to be in a `bin/` folder. To compile run `nextflow run vcab.nf` on the command line after setting up nextflow.

## Access VCAb online

VCAb can be accessed [online](https://fraternalilab.cs.ucl.ac.uk/VCAb/). For documentation on how to use the website, please go to the [github wiki](https://github.com/Fraternalilab/VCAb/wiki).

## Generate local version of VCAb

### VCAb shiny application
To generate the VCAb shiny app, go to `shiny` and run `app.R`.
#### Note
In order to run this code, BLAST must be installed in the command line. Please go to [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) for more information.

### VCAb antibody structural space
To collect the available experimental structural space of antibodies with both V and C regions, go to the `vcab_db` directory, run `generate_db.py`.
#### Note 
1. `generate_db.py` should be run in python 3.
2. In order to run this code, BLAST must be installed in the command line. Please go to [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) for more information.
2. `generate_db.py` would download PDB files and generate POPSComp results automatically by running bash scripts. Please make sure the `pops.sh`, `deltaResidue.R` (Modified from POPSR, please go to [github](https://github.com/Fraternalilab/POPScomp/tree/master/POPSR) for more information), `pdb_download.sh` (Modified from the script provided by [RCSB pdb](https://www.rcsb.org/docs/programmatic-access/batch-downloads-with-shell-script)) are in the corresponding directories.
3. `seq_db` hold all the BLAST databases. Two files are under this directory: `ref_db` contains reference sequences for isotypes and light chain types collected from uniprot and won't be changed. `vcab_db` includes sequences in VCAb, which would be updated when `generate_db.py` is run.
4. anarci_vc is a package modified from anarci, used for antibody numbering for both V and C sequences. Please check [github page](https://github.com/Fraternalilab/ANARCI_vc) for more information.


### Packages required 
| file | Package | version |
| ---- | ------- | ------- |
|`generate_db.py`| pandas | 1.2.1|
|`generate_db.py`| numpy | 1.19.2|
|`generate_db.py`| Bio | 1.79|
|`generate_db.py`| argparse | 1.1|
|`get_paired_seq.py`|anarci_vc|modified from anarci, see [github page](https://github.com/Fraternalilab/ANARCI_vc)|
|`app.R`| shiny | 1.6.0|
|`app.R`| DT |0.17|
|`app.R`| rBLAST |0.99.2|
|`app.R`| taxonomizr |0.8.0|
|`app.R`| NGLVieweR |1.3.2|
|`app.R`| tibble |3.1.0|
|`app.R`| shinyhelper |0.3.2|
|`app.R`| ggplot2 |3.3.3|
|`app.R`| dplyr |1.0.7|
