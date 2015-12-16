dosumm <- function(psub) {
	index <- unique(psub$index_snp)
	chr <- unique(psub$CHROM)
        n95 <- which(psub$cpp>=0.95)[1]
        n99 <- which(psub$cpp>=0.99)[1]
	region.snps=nrow(psub)
	region.length=diff(range(psub$POS))

	df <- data.frame(index,chr,n95,n99,region.snps,region.length, stringsAsFactors=FALSE)
        return(df)
}

pproc <- function(psub, index) {
	psub$index <- psub$SNPID==index
	# if an indexSNP can not be found print an error
	if(sum(psub$index)!=1){
	return(FALSE)
	# stop()
	} else {
	# calculate posterior probability of the SNPs for each indexSNP
	psub$pp <- exp(psub$BF - coloc:::logsum(psub$BF))
	# sort the SNPs based on the posterior probability
	psub <- psub[order(psub$pp,decreasing=TRUE),]
	# take all SNPs which are in 99% cumulative sum 
	psub$cpp <- cumsum(psub$pp)
	}

	return(psub)
}
