suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(coloc))
library(readr)
library(tidyr)
library(stringr)
library(optparse)

# Locate abf.R regardless of where this script was run from
# See http://stackoverflow.com/a/6461822
source_local <- function(fname){

    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))

}
source_local('abf.R')

option_list = list(
	make_option(c("-r", "--regions"), type="character", default=NULL,
		help="List of regions file"),
	make_option(c("-a", "--affected"), type="integer", default=NULL,
                help="Number of cases"),
	make_option(c("-u", "--unaffected"), type="integer", default=NULL,
                help="Number of controls"),
	make_option(c("-s", "--stats"), type="character", default=NULL,
		help="Summary statistics file"),
	make_option(c("-c", "--cpp"), type="double", default=0.99,
                help="Cumulative posterior probability")
) 

opt <- parse_args(OptionParser(option_list=option_list))

# read regions file
regions <- read_delim(opt$regions, delim = "\t", col_names = c('index','locus')) %>%
	separate(locus, sep = "[:-]", c('chr','start','end'))

# read summary statistics
data <- read_delim(opt$stats, delim = " ")

# calculate ABF
data$BF <- -abf(p=data$PVAL, maf=data$A1_UNAFF, n0=opt$affected, n1=opt$unaffected)

# define credible SNP sets

cred <- summ  <- vector("list",nrow(regions))

for(i in seq_along(cred)) {

        current_index <- as.character(regions[i,1])

	subset <- data %>%
			filter(index_snp==current_index) %>%
			mutate(index=SNPID==current_index) %>%
			mutate(pp=exp(BF - coloc:::logsum(BF))) %>%
			arrange(desc(pp)) %>%
			mutate(cpp=cumsum(pp))

	cred[[i]] <- subset %>%
			slice(1:which.max(cpp >= opt$cpp)) %>%
			dplyr::select(-BF, -index)

	# create summary information
	index <- current_index
	chr <- unique(subset$CHROM)
	region_n <- nrow(subset)
	region_length <- diff(range(subset$POS))
	cred_n <- nrow(cred[[i]])
	cred_length <- diff(range(cred[[i]]$POS))
	cred_start <- min(cred[[i]]$POS)
	cred_end <-max(cred[[i]]$POS) 

	summ[[i]] <- data.frame(index,chr,region_n,region_length,cred_n, cred_length,cred_start,cred_end, stringsAsFactors=FALSE)

}

# PJB need to create output dir if not present
output_dir <- "output/credible_snps/"
if (!file.exists(output_dir)) {

        print(paste("Creating output dir:",output_dir))
        dir.create(output_dir)

}

# create summary table
summary_table <- rbind_all(summ)
summary_file <- paste("output/credible_snps/summary_table_",opt$cpp,".txt",sep="")
write_delim(summary_table, delim = " ", summary_file)

# create credible snp table
credible_snps <- do.call("rbind",cred)
credible_file <- paste("output/credible_snps/credible_snps_",opt$cpp,".txt",sep="")
write_delim(credible_snps, delim = " ", credible_file)

# create data for bed file tracks
cred_region <- summary_table %>%
        dplyr::select(chr, cred_start, cred_end, index) %>%
        mutate(chr = paste('chr',chr,sep=''))

cred_snps <- credible_snps %>%
        mutate(start = POS - 1) %>%
        mutate(chr = paste('chr',CHROM,sep='')) %>%
        arrange(CHROM, start) %>%
        dplyr::select(chr, start, POS, SNPID)

# create bed file
bed_file <- paste("output/credible_snps/credible_snps_",opt$cpp,".bed",sep="")
cat("track name=\"cred\" description=\"Cred interval\" visibility=1", file = bed_file, sep = "\n")
write_delim(cred_region, bed_file, delim = " ", col_names = FALSE, append = TRUE)
cat("track name=\"credSNPs\" description=\"Cred SNPs\" visibility=1", file = bed_file, sep = "\n", append = TRUE)
write_delim(cred_snps, bed_file, delim = " ", col_names = FALSE, append = TRUE)
