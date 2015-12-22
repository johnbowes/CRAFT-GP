# CRAFT-GP
It is widely accepted that the most strongly associated SNP at any particular risk locus, referred to as the index SNP, is unlikely to be the causal variant, but rather, a highly correlated proxy for the true causal SNP. High-density imputed genotype datasets drastically increase the chance that the true causal SNP is within this correlated set of SNPs. These dense genotype datasets can be used to fine-map genetic risk loci by statistically refining the association signal to determine the subset of most informative SNPs, this subset is referred to as the credible SNP set [1]. 

This pipeline will apply a Bayesian refinement method to all risk loci that has been previously demonstrated to be more efficient at selecting putatively functional SNPs compared to linkage disequilibrium (LD) based approaches [2,3]. For each index SNP a genomic interval will be defined where the upstream and downstream boundaries are determined by a genetic distance of 0.1 centimorgans (cM) from the index SNP using HapMap fine-scale recombination rate estimates. For each of these intervals the posterior probability that any particular SNP is the casual SNP based on the Bayes factor for the SNP as a proportion of the mean Bayes factor for all SNPs in the genomic interval will be calculated. Posterior probabilities will be aggregated to define the smallest set of SNPs with a total posterior probability of ≥ 99% [2]. The credible SNP sets can be further refined with the integration of functional genomic annotations derived from disease-relevant tissue and cell types such as those described below [3].

1 	Farh KK-H, Marson A, Zhu J, et al. Genetic and epigenetic fine mapping of causal autoimmune disease variants. Nature Published Online First: 29 October 2014.http://www.ncbi.nlm.nih.gov/pubmed/25363779 (accessed 29 Oct2014).
2 	Maller JB, McVean G, Byrnes J, et al. Bayesian refinement of association signals for 14 loci in 3 common diseases. Nat Genet 2012;44:1294301. doi:10.1038/ng.2435
3 	Onengut-Gumuscu S, Chen W-M, Burren O, et al. Fine mapping of type 1 diabetes susceptibility loci and evidence for colocalization of causal variants with lymphoid gene enhancers. Nat Genet Published Online First: 2015. doi:10.1038/ng.3245

The pipeline consists of four stages:
1. define the region around each index SNP
2. calculate creadible SNPs within each interval
3. annotate credible SNP set
4. visualise results

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
$ ruby scripts/define_regions_main.rb -i input/test_index_snps.txt -p input/test_summary_statistics.txt -m 0.1

-i input (for the index SNPs)
-p pvalues (for the pvalues of the SNPs)
-m cM (normally 0.1(default) or 0.2)
-r bp region (normally 50kb)

```

### Output

1.
2.

## Calculate credible SNP sets

### subset GWAS data
python scripts/filter_summary_stats.py --regions output/regions/region_boundaries_01cm.txt --stats input/test_summary_stats.txt --out output/credible_snps/test_summary_stats_subset.txt

### calculate credible SNP
Rscript scripts/credible_snps_main.R -r output/regions/region_boundaries_01cm.txt -a 1962 -u 8923 -s output/credible_snps/test_summary_stats_subset.txt

