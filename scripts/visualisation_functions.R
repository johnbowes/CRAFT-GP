plot_region <- function(region,snps,list, out_dir, gen){
  
  # buffer base pairs to add to plot
  buffer <- 1000
  t_width <- 3
  
  # create output file names
  #pdf_file <- paste(out_dir, "/", region$index, ".pdf", sep="")
  png_file <- paste(out_dir, "/", region$index, ".png", sep="")
  
  # base tracks
  ideogram_track <- IdeogramTrack(genome = gen, chromosome = region$chr, bands = cytobands)
  genome_track <- GenomeAxisTrack(range=IRanges(start=region$cred_start, end=region$cred_end, names="credible interval"), showId=T)
  
  # data track
  pp <- dplyr::select(snps, chromosome, start, end, pp)
  pp_track <- DataTrack(pp, genome = gen, start = region$cred_start, end = region$cred_end, name = "posterior probability")
  displayPars(pp_track) = list(background.title="transparent", col.title="black", col.axis="black", cex.axis=0.8, cex.title=0.8)
                                 
  # credible region tracks - incorporated into genome axis track
  #interval_track <- AnnotationTrack(start = region$cred_start, end = region$cred_end,
  #                                  chromosome = region$chr, name = "Interval", ucscChromosomeNames = TRUE)
  
  
  snps <- dplyr::select(snps, chromosome, start, end, id)
  snp_track <- AnnotationTrack(snps, stacking="dense", name = "SNPs")
  displayPars(snp_track) = list(rotation.title=360, background.title="transparent", col.title="black")
  
  # gene information
  gene_track <- BiomartGeneRegionTrack(genome = gen, chromosome = region$chr,
                                       start = region$region_start, end = region$region_end, name = "ENSEMBL", biomart = ensembl,
                                       transcriptAnnotation = "symbol", collapseTranscripts = "longest", ucscChromosomeNames = TRUE)
  
  displayPars(gene_track) = list(rotation.title=360, showId=TRUE,stackHeight=0.5, background.title="transparent", col.title="black", rotation.tile=360)
  
  # epigenome data
  track_list <- list(ideogram_track, genome_track, pp_track, snp_track)
  
  
  for(x in list){
    epi <- filter(get(x), chromosome == region$chr &
                    (start >= region$cred_start | end >= region$cred_start) &
                    (start <= region$cred_end | end <= region$cred_end))
    
    
    # crude fix to deal with long track names - split names with 6 or more words
    track_name <- x
    num_words <- lengths(str_split(track_name, "_"))
    
    if (num_words > 6){
      words <- str_split(track_name, "_")[[1]]
      split <- ceiling(num_words/2)
      string1 <- str_c(words[1:split], collapse="_")
      string2 <- str_c(words[(split+1):num_words], collapse="_")
      track_name <- paste(string1, string2, sep=" ")
      t_width <- 4
    }
      
    #y <- AnnotationTrack(epi, stacking = "dense", name = track_name, legend=TRUE, col=epi$feature, groups=epi$id)
    y <- AnnotationTrack(epi, name = track_name)
    displayPars(y) = list(rotation.title = 360, cex.title=0.5, stacking="dense", shape="box",stackHeight=1, background.title="transparent", col.title="black", col=NULL)

    # colour features
    feat <- unique(feature(y))
    featCol <- setNames(as.list(rgb(t(sapply(strsplit(feat, ","), as.numeric)), maxColorValue=255)), feat)
    displayPars(y) <- featCol
    
    # assign epigenome name to track
    track_name <- paste(x,"_track", sep="")
    assign(track_name, y)
    track_list <- c(track_list, get(track_name))
    
  }
  
  # add gene track after epigenome tracks
  track_list <- c(track_list, gene_track)
  
  # plot
  #pdf(pdf_file)
  png(png_file)
  plotTracks(track_list,from = (region$cred_start - buffer), to = (region$cred_end + buffer), title.width=t_width)
  dev.off()
  
}

read_eid <- function(eid,dir){
  
  file <- paste(dir, "/",eid,".bed.gz", sep="")
  
  df <- read_delim(file, delim="\t", col_names=FALSE, col_types = cols(X9 = "c")) %>%
    dplyr::select(-X5, -X6, -X7, -X8) %>%
    as.data.frame()
  
  colnames(df) <- c("chromosome","start","end","id","feature")

  return(df)

}
