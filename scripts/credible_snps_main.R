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
                help="Cumulative posterior probability"),
	make_option(c("-o", "--out"), type="character", default=NULL,
                help="Directory path for saving results"),
	make_option(c("-b", "--bed"), type="character", default=NULL,
	            help="Directory path for saving bed files")
) 

opt <- parse_args(OptionParser(option_list=option_list))

if(interactive())
  opt <- list(regions = "test_results/regions/region_boundaries_01cm.txt",
  			      affected = 1962,
              unaffected = 8923,
              stats = "test_results/credible_snps/summary_stats_subset.txt",
              cpp = 0.99,
              out = "test_results/credible_snps/",
              bed = "test_results/bed/")

# read regions file
regions <- read_delim(opt$regions, delim = "\t", col_names = c('index','locus'), col_types=cols(index=col_character())) %>%
	separate(locus, sep = "[:-]", c('chr','start','end'))

# read summary statistics
data <- read_delim(opt$stats, delim = " ", col_types=cols(SNPID=col_character()))

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
	chr <- paste("chr", unique(subset$CHROM), sep = "")
	region_n <- nrow(subset)
	region_length <- diff(range(subset$POS))
  region_start <- regions[i,]$start
  region_end <- regions[i,]$end
	cred_n <- nrow(cred[[i]])
	cred_length <- diff(range(cred[[i]]$POS))
	cred_start <- min(cred[[i]]$POS)
	cred_end <-max(cred[[i]]$POS) 

	summ[[i]] <- data.frame(index,chr,region_n,region_length,region_start,region_end,cred_n, cred_length,cred_start,cred_end, stringsAsFactors=FALSE)

}

# PJB need to create output dir if not present
#output_dir <- "output/credible_snps/"
#if (!file.exists(output_dir)) {
#
#        print(paste("Creating output dir:",output_dir))
#        dir.create(output_dir)
#
#}

# create summary table
summary_table <- rbind_all(summ)
summary_file <- file.path(opt$out, paste("summary_table_", opt$cpp,".txt",sep=""))
write_delim(summary_table, delim = " ", summary_file)

# create credible snp table
credible_snps <- do.call("rbind",cred)
credible_file <- file.path(opt$out, paste("credible_snps_", opt$cpp,".txt",sep=""))
write_delim(credible_snps, delim = " ", credible_file)

# create list of credible SNP rs numbers for VEP input
<<<<<<< HEAD
#credible_snp_list <- credible_snps %>%
#	dplyr::select(SNPID)
#
#credible_snp_list_file = paste(opt$out, "credible_snp_list_", opt$cpp,".txt",sep="")
#write_delim(credible_snp_list, delim = " ", credible_snp_list_file, col_names = FALSE)
=======
credible_snp_list <- credible_snps %>%
	dplyr::select(SNPID)

credible_snp_list_file = file.path(opt$out, paste("credible_snp_list_", opt$cpp,".txt",sep=""))
write_delim(credible_snp_list, delim = " ", credible_snp_list_file, col_names = FALSE)
>>>>>>> c27396323a0486fc42d0e749bbb88dded54282e1

# create data for bed file tracks
cred_region <- summary_table %>%
        dplyr::select(chr, cred_start, cred_end, index)

cred_snps <- credible_snps %>%
        mutate(start = as.integer(POS - 1)) %>%
        mutate(chr = paste('chr',CHROM,sep='')) %>%
        arrange(CHROM, start) %>%
        dplyr::select(chr, start, POS, SNPID)

# write sessionInfo to file
writeLines(capture.output(sessionInfo()), "sessionInfo_cred_snps.txt")

# save image for comparing workspace objects
save.image("credible_snps_image.RData")

# create bed file
<<<<<<< HEAD
bed_file <- paste(opt$bed, "credible_snps_", opt$cpp,".bed",sep="")
cat("track name=\"cred\" description=\"Cred interval\" visibility=1", file = bed_file, sep = "\n")
write_delim(cred_region, bed_file, delim = " ", col_names = FALSE, append = TRUE)
cat("track name=\"credSNPs\" description=\"Cred SNPs\" visibility=1", file = bed_file, sep = "\n", append = TRUE)
=======
bed_file <- file.path(opt$out, paste("credible_snps_", opt$cpp,".bed",sep=""))
#cat("track name=\"cred\" description=\"Cred interval\" visibility=1", file = bed_file, sep = "\n")
#write_delim(cred_region, bed_file, delim = " ", col_names = FALSE, append = TRUE)
cat("track name=\"credSNPs\" description=\"Cred SNPs\" visibility=1", file = bed_file, sep = "\n")
>>>>>>> c27396323a0486fc42d0e749bbb88dded54282e1
write_delim(cred_snps, bed_file, delim = " ", col_names = FALSE, append = TRUE)

