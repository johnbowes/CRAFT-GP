suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(coloc))
library(readr)
library(tidyr)
library(stringr)
library(optparse)

source('scripts/credible_snps_functions.R')
source('scripts/abf_corrected.R')

option_list = list(
	make_option(c("-r", "--regions"), type="character", default=NULL,
		help="List of regions file"),
	make_option(c("-a", "--affected"), type="integer", default=NULL,
                help="Number of cases"),
	make_option(c("-u", "--unaffected"), type="integer", default=NULL,
                help="Number of controls"),
	make_option(c("-s", "--stats"), type="character", default=NULL,
		help="Summary statistics file")
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

cred99 <- cred95 <- summ  <- vector("list",nrow(regions))

for(i in seq_along(cred99)) {

        current_index <- as.character(regions[i,1])
        psub <- pproc(subset(data,index_snp==current_index), current_index)

        sm <- dosumm(psub)
	cred95[[i]] <- psub[1:sm$n95, c("SNPID","CHROM","POS","A1_UNAFF","PVAL","index_snp","pp","cpp")]
        cred99[[i]] <- psub[1:sm$n99, c("SNPID","CHROM","POS","A1_UNAFF","PVAL","index_snp","pp","cpp")]
        summ[[i]] <- sm

}

# create summary table
summary <- rbind_all(summ) %>%
	mutate(cred95.start = sapply(cred95,function(x) min(x$POS))) %>%
	mutate(cred95.end = sapply(cred95,function(x) max(x$POS))) %>%
	mutate(cred99.start = sapply(cred99,function(x) min(x$POS))) %>%
	mutate(cred99.end = sapply(cred99,function(x) max(x$POS))) %>%
	mutate(cred95.interval = cred95.end - cred95.start) %>%
	mutate(cred99.interval = cred99.end - cred99.start)

write_delim(summary, delim = " ", "output/credible_snps/summary_table.txt")

# create credible snp table
cred99.long <- do.call("rbind",cred99)
cred95.long <- do.call("rbind",cred95)
write_delim(cred99.long, delim = " ", "output/credible_snps/credible_snps.txt")

# create data for bed file tracks
cred95_region <- summary %>% dplyr::select(chr, cred95.start, cred95.end, index)
cred99_region <- summary %>% dplyr::select(chr, cred99.start, cred99.end, index)

cred95_snps <- cred95.long %>%
	mutate(start = POS - 1) %>%
	mutate(chr = as.integer(str_replace(CHROM, 'chr', ''))) %>%
	arrange(chr, start) %>%
	dplyr::select(CHROM, start, POS, SNPID)

cred99_snps <- cred99.long %>%
        mutate(start = POS - 1) %>%
        mutate(chr = as.integer(str_replace(CHROM, 'chr', ''))) %>%
        arrange(chr, start) %>%
        dplyr::select(CHROM, start, POS, SNPID)

# create bed file
bed_file <- "output/credible_snps/credible_snps.bed"
cat("track name=\"cred99\" description=\"Cred99 interval\" visibility=1", file = bed_file, sep = "\n")
write_delim(cred99_region, bed_file, delim = " ", col_names = FALSE, append = TRUE)
cat("track name=\"cred95\" description=\"Cred95 interval\" visibility=1", file = bed_file, sep = "\n", append = TRUE)
write_delim(cred95_region, bed_file, delim = " ", append = TRUE)
cat("track name=\"cred99SNPs\" description=\"Cred99 SNPs\" visibility=1", file = bed_file, sep = "\n", append = TRUE)
write_delim(cred99_snps, bed_file, delim = " ", col_names = FALSE, append = TRUE)
cat("track name=\"cred95SNPs\" description=\"Cred99 SNPs\"  visibility=1", file = bed_file, sep = "\n", append = TRUE)
write_delim(cred99_snps, bed_file, delim = " ", col_names = FALSE, append = TRUE)
