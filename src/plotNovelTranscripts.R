# Find novel ZFX transcripts
library(tidyverse)
library(ggtranscript) # devtools::install_github("dzhang32/ggtranscript")
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(GenomicFeatures)
library(BSgenome)

cat("Reading genome FASTA\n")
chicken.genome <- Biostrings::readDNAStringSet("genomes/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz")
names(chicken.genome) <- c(1:214) # replace FA heaeder with chr names. We only need chr1 anyway

# Read gtf file and convert to tibble. Export FASTA of each transcript
read.novel.transcripts <- function(gtf.file){
  chicken.novel.gtf <- rtracklayer::import(gtf.file)
  # for(transcript.id in unique(chicken.novel.gtf$transcript_id)){
  #   exon.gtf <- GenomicRanges::GRangesList(chicken.novel.gtf[chicken.novel.gtf$transcript_id==transcript.id & chicken.novel.gtf$type=="exon",])
  #   if(any(strand(exon.gtf)=="*")){
  #     cat("Cannot get transcript sequence with ambiguous strand\n")
  #     next
  #   }
  #   chicken.seqs <- GenomicFeatures::extractTranscriptSeqs(chicken.genome, exon.gtf)
  #   names(chicken.seqs) <- paste0(gtf.file, transcript.id)
  #   writeXStringSet(chicken.seqs, filepath = paste0(gtf.file, ".", transcript.id, ".fa"))
  # }
  as_tibble(chicken.novel.gtf) %>%
    dplyr::mutate(file = gtf.file)
}

# Read gtf file andplot. Export FASTA of each transcript
plot.novel.transcripts <- function(gtf.file){
  cat("Plotting", gtf.file, "\n")
  chicken.novel.gtf <- rtracklayer::import(gtf.file)
  chicken.novel.gtf.tbl <- as_tibble(chicken.novel.gtf)
  chicken.novel.gtf.exon <- chicken.novel.gtf.tbl %>% dplyr::filter(type=="exon") 
  chicken.novel.gtf.intron <- to_intron(chicken.novel.gtf.exon, "transcript_id")
  
  transcript.plot <- chicken.novel.gtf.tbl %>% dplyr::filter(type=="exon") %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = transcript_id
    )) +
    geom_intron(
      data = chicken.novel.gtf.intron,
      aes(strand = strand, col=strand)
    )+
    geom_range(aes(fill = as.numeric(TPM))
    ) +
    scale_color_manual(values = c("-"="black", "+"="red"))+
    coord_cartesian(xlim = c(118319000, 118296000))+
    labs(x = "Position", y="Assembled transcript", title=gtf.file, fill="TPM", col="Strand")+
    geom_text(
      data = add_exon_number(chicken.novel.gtf.exon, "transcript_id"),
      aes(
        x = (start + end) / 2, # plot label at midpoint of exon
        label = exon_number
      ),
      size = 3.5,
      nudge_y = 0.4
    )+
    theme(title = element_text(size=3))
  
  ggsave(paste0(gtf.file, ".transcripts.png"), last_plot())
  
  transcript.plot
  
}

files <- list.files(path = "data/stringtie",pattern = "*.gtf$", full.names = TRUE)

gtf.data <- do.call(rbind, lapply(files, read.novel.transcripts))
plots <- lapply(files,plot.novel.transcripts )
combined.plot <- patchwork::wrap_plots(plots, ncol = 3)+ patchwork::plot_layout(guides = "collect", axes = "collect", axis_titles = "collect")
ggsave("data/stringtie/chicken.transcripts.all.png", combined.plot, units = "mm", height = 230, width = 170, dpi = 300)
# BLAST of the sequence hits https://www.ncbi.nlm.nih.gov/gene/121109097
# It's not in the Ensembl genebuild because it's a predicted transcript only, so the reads were aggregated into Zfx during feature mapping
