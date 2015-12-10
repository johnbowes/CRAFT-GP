# CRAFT-GP

## Define disease loci
A region around each index SNP is defined based on recombination. This process also subsets the
supplied GWAS summary statistics for SNPs within the defined region.

### Input

1. Index SNPs

A list of index SNPs for each locus comprised of three columns and one row per index SNP.

```
rs	   chr	position
rs2301888   1	17672730
rs28411352  1	38278579
rs12140275  1	38633879

```

2. Summary statistics
Input file is tab, space or comma-delimited and defined as folllows:

```
SNPID	   CHROM  POS       A1	A2  A1_UNAFF  PVAL
rs6823274   4	  26098057  G   A   0.1484    0.4064
rs76632663  2	  43568820  T   G   0.06988   0.4427
rs28719598  8	  10943884  T   C   0.194     0.7702
rs78519860  6	  128145388 T   C   0.07869   0.9007
```
 
### Run

```
$ ruby scripts/regions_main.rb -i input/test_index_snps.txt -p input/test_summary_stats.txt -m 0.1

-i input (for the index SNPs)
-p pvalues (for the pvalues of the SNPs)
-m cM (normally 0.1(default) or 0.2)
-r bp region (normally 50kb)

```

### Output

1.
2.

## Calculate credible SNP sets
