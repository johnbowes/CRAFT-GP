CRAFT
=====

Acknowledgements
----------------

Requirements
------------

__iCSF modules__

module load apps/gcc/ruby/2.2.3
module load apps/binapps/anaconda/2.1.0
module load apps/gcc/R/3.3.0
module load apps/binapps/ensembl-api/84

source("https://bioconductor.org/biocLite.R")
biocLite("Gviz")

User input
----------

NAME=''                # name for labelling output files
INDEX_SNP_FILE=''      # file name for list of index SNPs
SUMMARY_STATS_FILE=''    # file name for GWAS summary statistics
DISTANCE_TYPE='-m'   # how the distance should be defined (base pairs (-r) of centimorgans (-m))
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


__Additional variables__


Define genomic intervals
------------------------
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
