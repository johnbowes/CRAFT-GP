CRAFT
=====

**C**redible **R**efinement and **A**nnotation of **F**unctional **T**argets (CRAFT) is a pipeline for the calculation, annotation and visualisation of credible SNP sets. The scripts in this repoistory form the basis of development for a Galaxy implementation distributed using the Galaxy Tool Shed, but can be used as a stand-alone pipeline.

Acknowledgements
----------------

CRAFT uses the abf.R function written by [Chris Wallace](http://chr1swallace.github.io/) for the calculation of credible SNPs and can be found [here](http://chr1swallace.github.io/).

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

* module load apps/gcc/ruby/2.2.3  
* module load apps/binapps/anaconda/2.3.0  
* module load apps/gcc/R/3.3.0  
* module load apps/binapps/ensembl-api/84  
* module load apps/gcc/tabix/0.2.6  

Source data
-----------

__Variant Effect Predictor (VEP)__

Get the GRCh37 (release 84) core cache (not included due to size)
```
mkdir -p source_data/ensembl/{cache,plugins}
cd source_data/ensembl/cache/
wget ftp://ftp.ensembl.org/pub/release-84/variation/VEP/homo_sapiens_vep_84_GRCh37.tar.gz
tar -zxvf homo_sapiens_vep_84_GRCh37.tar.gz
cd -
```

Download CADD plugin data
```
mkdir -p source_data/ensembl/plugins/CADD
cd source_data/ensembl/plugins/CADD/
wget http://krishna.gs.washington.edu/download/CADD/v1.3/1000G.tsv.gz
wget http://krishna.gs.washington.edu/download/CADD/v1.3/1000G.tsv.gz.tbi
cd -
```

__Roadmap Epigenomics__

15 state model  
```
# get data

mkdir -p source_data/roadmap_r9/15_state_model/{raw,bed}
cd source_data/roadmap_r9/15_state_model/raw/
wget http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all.dense.browserFiles.tgz
tar -zxvf all.dense.browserFiles.tgz
cd -

# files are renamed and tabix indexed for use with VEP
python scripts/process_roadmap.py --state 15 --meta roadmap_consolidated_epigenome_ids_blood.csv --rename EDACC
```

__HapMap recombination map__

```
mkdir -p source_data/genetic_map_HapMapII_GRCh37/
cd source_data/genetic_map_HapMapII_GRCh37/
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz
tar -zxvf genetic_map_HapMapII_GRCh37.tar.gz
cd -
```

__UCSC cytogenetic bands__

```
mkdir -p source_data/ucsc/
cd source_data/ucsc/
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz
gunzip cytoBand.txt.gz
cd -
```

__CRAFT source data bundle__

Currently working on a data bundle to distribute all required source data in one package.

Test data
---------

Current test data is three regions from the recent PsA immunochip study (Bowes *et al*, Nature Communications 2015). A full run through of the pipeline using the test data is given in the run_all.sh script.

User variables
--------------

These are currently defined witihin the run_all.sh script:

*NAME*  
    Name for labelling output directory and files  

*INDEX_SNP_FILE*  
    File name for list of index SNPs, see below for format  

*SUMMARY_STATS_FILE*  
    File name for GWAS summary statistics, see below for format  

*DISTANCE_TYPE*  
    A flag to specifiy how the distance should be defined, this can either be base pairs (-r) of centimorgans (-m)  

*DISTANCE*  
    Size of region as defined above, typical values for cM are 0.1 0r 0.2  

*CASES*  
    Number of cases in GWAS summary statistics  

*CONTROLS*  
    Number of controls in GWAS summary statistics  

*EPIGENOMES*  
    List of roadmap epigenomics dataset to annotation against    
    

User input file formats
-----------------------

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

Output
------

Output is returned to the user at multiple levels:

Plots

Bed tracks/hub

SNP level annotation

Define genomic intervals
------------------------

__input__

*INDEX_SNP_FILE: user supplied file name of index SNPs  

__output__  

Details of region definitions are written to *NAME*\_results/regions/region_boundaries\_01cm.txt

```
rs2224873   chr1:197311728-197938007
rs7536201   chr1:25227739-25307539
rs9988642   chr1:67596221-67766669
```

__Run__

```
ruby scripts/define_regions_main.rb \
    -i $INDEX_SNP_FILE \                        # User supplied file name
    $DISTANCE_TYPE $DISTANCE \                  # User supplied distance parameters
    -o $OUTPUT                                  # output path for regions based on user supplied name
```

Calculate credible SNP sets
---------------------------

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

__output__

*NAME*\_results/credible_snps/credible\_snp\_list\_0.99.txt  
*NAME*\_results/credible_snps/credible\credible_snps_0.99.bed    
*NAME*\_results/credible_snps/credible\credible_snps_0.99.txt  
*NAME*\_results/credible_snps/credible\summary_table_0.99.txt  

* GWAS_SUBSET: temporary file containing subset of GWAS summary statistics for defined regions  

Annotation
----------

**Note:** the current version uses assemlby GRCh37 and VEP version 84 it therefore requires access to port 3337 on ensembldb.ensembl.org. When using iCSF this has been handled by adding ensembldb.ensembl.org:3337 to the NAT server and using the options --host=130.88.97.228 --port=3337 when calling VEP.

__Run__

```
python scripts/annotation.py \
    --input $CREDIBLE_SNPS \                    # credible snp output from previous stage
    --output $OUTPUT \                          # output path based on user suppled name
    --eid $EPIGENOMES                           # User supplied list of Roadmap Epigenomics cell types
```


__output__

*NAME*\_results/plots/rs\*

Visualisation
-------------

__Run__

```
Rscript scripts/visualisation_main.R \
    -r ${REGIONS} \                             # List of defined regions output from stage 1
    -s ${SNPS} \                                # Bed file of credible SNPs output from stage 2
    -e $EPIGENOMES \                            # User supplied list of Roadmap Epigenomics cell types
    -o ${OUTPUT}                                # output path based on user suppled name
```

__output__

Known issues
------------

* epigenome names are currently case sensitive
* gwas summary stats must be space sepatated
* No error catching implemented

TODO
----

__Before release__

*offline mode*
* test full offline mode with VCF input

*input*
* file delimiters
* strip out non-rs numbers - display warning - not necessary with offline mode?
* decide on use of chromosome label (either chr1 or 1)
*
*output*
* report non-annotated SNPs
* Are all output used/useful?
    - esp. cred SNP stage
* create data bundle
    - where to host

__Nice to have__
* wrap the whole thing in a python wrapper
* output list of credible SNPs with no annotation
* input options - PLINK
* get MAF from reference panel
* convert cache to tabix
    - https://gist.github.com/ckandoth/57d189f018b448774704d3b2191720a6