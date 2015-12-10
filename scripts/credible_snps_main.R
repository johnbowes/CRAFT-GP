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
	separate(locus, sep = "[:-]", c('chr','start','end')

# read summary statistics
data <- read_delim(opt$stats, delim = "\t")

	
# calculate ABF
data$BF <- -abf(p=data$PVAL, maf=data$MAF, n0=opt$affected, n1=opt$unaffected)

# define credible SNP sets
cred99 <- cred95 <- summ  <- vector("list",1)

psub <- pproc(data, opt$snp) %>%
	mutate(index_snp=opt$snp)

sm <- dosumm(psub)

cred95 <- psub[1:sm["n95"], c("SNPID","index_snp","CHR","POS","MAF","PVAL","index","pp","cpp")]
cred99 <- psub[1:sm["n99"], c("SNPID","index_snp","CHR","POS","MAF","PVAL","index","pp","cpp")]
summ <- c(narrow=sm)

# write results to file
cred95_file <- paste(opt$snp, "credible_snps_95.txt", sep='_')
cred99_file <- paste(opt$snp, "credible_snps_99.txt", sep='_')

write_tsv(cred95, cred95_file)
write_tsv(cred99, cred99_file)
