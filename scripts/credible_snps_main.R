suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(coloc))
library(readr)
library(tidyr)
library(stringr)
library(optparse)

source('credible_snps_functions.R')

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
data <- read_delim(opt$stats, delim = "\t")

# calculate ABF
data$BF <- -abf(p=data$PVAL, maf=data$MAF, n0=opt$affected, n1=opt$unaffected)

# define credible SNP sets

cred99 <- cred95 <- summ  <- vector("list",nrow(regions))

for(i in seq_along(cred99)) {

	current_index <- as.character(regions[i,1])
	psub <- pproc(subset(data,indexSNP==current_index), current_index)

	sm <- dosumm(psub)
	cred95[[i]] <- psub[1:sm["n95"], c("SNPID","CHROM","POS","MAF","PVAL","indexSNP","pp","cpp")]
	cred99[[i]] <- psub[1:sm["n99"], c("SNPID","CHROM","POS","MAF","PVAL","indexSNP","pp","cpp")]
	summ[[i]] <- c(narrow=sm)

}

# create summary table
result <- data.frame(indexSNP = un_inSNP$indexSNP, do.call("rbind",summ), stringsAsFactors =FALSE)
result$cred95.start <- sapply(cred95,function(x) min(x$POS))
result$cred95.end <- sapply(cred95,function(x) max(x$POS))
result$cred99.start <- sapply(cred99,function(x) min(x$POS))
result$cred99.end <- sapply(cred99,function(x) max(x$POS))
result$cred95.interval <- result$cred95.end - result$cred95.start
result$cred99.interval <- result$cred99.end - result$cred99.start

write_delim(result, "output/credible_snps/summary_table.txt")

# create credible snp table
cred99.long <- do.call("rbind",cred99)

write_delim(cred99.long, "output/credible_snps/credible_snps.txt")

# create bed files

#write_delim(cred99.long, "output/credible_snps/credible_snps.bed")
