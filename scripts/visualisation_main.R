library(Gviz)
library(dplyr)
library(readr)
library(biomaRt)
library(optparse)
library(stringr)

# Locate visualisation_functions.R regardless of where this script was run from
# See http://stackoverflow.com/a/6461822
source_local <- function(fname){

    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))

}
source_local('visualisation_functions.R')

# OPTIONS

# eid option has to to expect the list of epigenomes in a comma separated
# string from a file.

option_list = list(
  make_option(c("-r", "--regions_file"), type="character", default=NULL,
              help="Region summary table"),
  make_option(c("-s", "--snp_file"), type="character", default=NULL,
              help="Credible SNPs"),
  make_option(c("-e", "--eid"), type="character", default=NULL,
              help="Epigenome to plot against"),
  make_option(c("-g", "--genome_build"), type="character", default="hg19",
              help="Genome build"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="Output file path")
) 

opt <- parse_args(OptionParser(option_list=option_list))

# list of test user arguments for testing
if(interactive())
  opt <- list(genome_build = "hg19",
              regions_file = "test_results/credible_snps/summary_table_0.99.txt",
              snp_file = "test_results/credible_snps/credible_snps_0.99.txt",
              #eid = "CD4_Naive_Primary_Cells,CD8_Memory_Primary_Cells,CD8_Naive_Primary_Cells,CD4+_CD25-_IL17-_PMA-Ionomycin_stimulated_MACS_purified_Th_Primary_Cells",
              eid = "test_results/annotation/annotation.epigenomes",
              out = "test_results/plots")

# SOURCE DATA - not user defined
cytoband_file <- "source_data/ucsc/cytoBand.txt"
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
e_dir <- "source_data/roadmap_r9/15_state_model/bed/"

# LOAD DATA

# region summary table
regions <- read_delim(opt$regions_file, col_names = TRUE, delim = " ")

# cytogenetic bands
cytobands <- read_delim(cytoband_file, delim="\t",
    col_names = c("chrom", "chromStart", "chromEnd", "name", "gieStain")) %>% as.data.frame()

# credible SNPs
#credible_snps <- read_delim(opt$snp_file, delim = " ", col_names = c("chromosome","start","end","id"), skip = 1) %>%
#  as.data.frame()

credible_snps <- read_delim(opt$snp_file, delim = " ") %>%
  mutate(start = POS - 1) %>%
  mutate(chromosome = paste("chr", CHROM, sep = "")) %>%
  dplyr::select(chromosome, start, POS, SNPID, pp) %>%
  rename(end = POS) %>%
  rename(id = SNPID)

# epigenome data
epi_list <- strsplit(read_file(opt$eid), ',')[[1]]
for(e in epi_list){
   assign(e, read_eid(e, e_dir))
}

# write sessionInfo to file
writeLines(capture.output(sessionInfo()), "sessionInfo_annotation.txt")

# save image for comparing workspace objects
save.image("annotation_image.RData")

# PLOT EACH REGION - multiple epigenomes
for (row in seq_len(nrow(regions))){
  region_list <- as.list(regions[row,])
  plot_region(region_list, credible_snps, epi_list, opt$out, opt$genome_build)
}


