plot_region <- function(region,snps,list, out_dir, gen){
  
  # create output file names
  pdf_file <- paste(out_dir, "/", region$index, ".pdf", sep="")
  
  # base tracks
  ideogram_track <- IdeogramTrack(genome = gen, chromosome = region$chr, bands = cytobands)
  genome_track <- GenomeAxisTrack()
  
  # credible region tracks
  interval_track <- AnnotationTrack(start = region$cred_start, end = region$cred_end,
                                    chromosome = region$chr, name = "Interval", ucscChromosomeNames = TRUE)
  
  snp_track <- AnnotationTrack(snps, stacking="dense", name = "Credible SNPs")
  
  # gene information
  gene_track <- BiomartGeneRegionTrack(genome = gen, chromosome = region$chr,
                                       start = region$region_start, end = region$region_end, name = "ENSEMBL", biomart = ensembl,
                                       transcriptAnnotation = "symbol", collapseTranscripts = "longest", ucscChromosomeNames = TRUE)
  
  displayPars(gene_track) = list(showId=TRUE,stackHeight=0.5)
  
  # epigenome data
  track_list <- list(ideogram_track, genome_track, interval_track, snp_track)
  
  
  for(x in list){
    epi <- filter(get(x), chromosome == region$chr &
                    (start >= region$cred_start | end >= region$cred_start) &
                    (start <= region$cred_end | end <= region$cred_end))
    
    track_name <- str_replace_all(x, "_", " ")
    y <- AnnotationTrack(epi, stacking = "dense", name = track_name)
    displayPars(y) = list(rotation.title = 360, cex.title = 0.3)

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
  pdf(pdf_file)
  plotTracks(track_list,from = region$region_start, to = region$region_end)
  dev.off()
  
}

read_eid <- function(eid,dir){
  
  file <- paste(dir, "/",eid,".bed.gz", sep="")
  
  df <- read_delim(file, delim="\t", col_names=FALSE) %>%
    dplyr::select(-X5, -X6, -X7, -X8) %>%
    as.data.frame()
  
  colnames(df) <- c("chromosome","start","end","id","feature")

  return(df)

}
