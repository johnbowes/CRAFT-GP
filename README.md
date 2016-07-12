CRAFT
=====

**C**redible **R**efinement and **A**nnotation of **F**unctional **T**argets (CRAFT) is a pipeline for the calculation, annotation and visualisation of credible SNP sets. The scripts in this repoistory form the basis for a Galaxy implenentation, but can be used as a stand-alone pipeline.

Acknowledgements
----------------

CRAFT using the abf.R function written by Chris Wallace () for the calculation of credible SNPs.

Requirements
------------

__Software__
* Variants Effect Predictor (VEP) version 84  
* Perl (tested on version 5.20.2)  
* Python (tested on anaconda 2.3.0)  
    - PyVCF  
    - pandas  
    - numpy  
    - argparse  
    - re
* Tabix (tested on version 0.2.6)  
* Ruby (tested on version 2.2.3)
* R (tested on version 3.3.0)
    - dplyr  
    - tidyr  
    - stringr 
    - Gviz  

__iCSF modules__

When running on the University of Manchester computing clusters (iCSF) the software dependencies can be dealt with by loading the following modules:

module load apps/gcc/ruby/2.2.3
module load apps/binapps/anaconda/2.3.0
module load apps/gcc/R/3.3.0
module load apps/binapps/ensembl-api/84
module load apps/gcc/tabix/0.2.6

Source data
-----------

Test data
---------

Current test data is three regions from the recent PsA immunochip study (Bowes *et al*, Nature Communications 2015). A full run through of the pipeline using the test data is given in the run_all.sh script.

User input
----------

__Run variables__

NAME                # name for labelling output files
INDEX_SNP_FILE      # file name for list of index SNPs
SUMMARY_STATS_FILE  # file name for GWAS summary statistics
DISTANCE_TYPE       # how the distance should be defined (base pairs (-r) of centimorgans (-m))
DISTANCE='0.1'                                  # distance value
CASES='1962'                                    # number of cases in GWAS summary statistics
CONTROLS='8923'                                 # number of controls in GWAS summary statistics
EPIGENOMES=''                                   # list of roadmap epigenomics dataset to annotation against

__Index SNP list__

A list of index SNPs for each locus comprised of three columns and one row per index SNP.

```
rs     chr  position
rs2301888   1   17672730
rs28411352  1   38278579
rs12140275  1   38633879

```

__GWAS summary statistics__

Input file is tab, space or comma-delimited and defined as folllows:

```
SNPID      CHROM  POS       A1  A2  A1_UNAFF  PVAL
rs6823274   4     26098057  G   A   0.1484    0.4064
rs76632663  2     43568820  T   G   0.06988   0.4427
rs28719598  8     10943884  T   C   0.194     0.7702
rs78519860  6     128145388 T   C   0.07869   0.9007
```

Define genomic intervals
------------------------

__input__

__output__

__Run__

```
python scripts/filter_summary_stats.py \
    --regions $REGIONS \                        # list of genomic intervals: output from previous step
    --stats input/$SUMMARY_STATS_FILE \         # GWAS summary statistics: user supplied
    --out $GWAS_SUBSET                          # ouput path for subset of GWAS statistics
```

```
Rscript scripts/credible_snps_main.R \
    -r $REGIONS \                               # list of genomic intervals: output from previous step
    -a $CASES \                                 # number of cases samples: user supplied
    -u $CONTROLS \                              # number of controls samples: user supplied
    -s $GWAS_SUBSET                             # subset of GWAS summary statistics
```

Calculate credible SNP sets
---------------------------
