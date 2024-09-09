# Find novel ZFX transcripts
library(tidyverse)
library(ggtranscript) # devtools::install_github("dzhang32/ggtranscript")
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(GenomicFeatures)
library(BSgenome)
library(patchwork)
source("src/functions.R")

#### Define functions ####
# Given a gzipped fasta file, read and ensure seqnames are just the chr name
read.reference.genome <- function(fa.gz.file){
  genome <-  Biostrings::readDNAStringSet(fa.gz.file)
  # replace FA header with only the chr names
  # Anything up to the first space
  names(genome) <- str_extract(names(genome), "^[\\w\\d]+")
  genome
}

# Read gtf file and convert to tibble. 
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

# Read gtf file and export FASTA of each transcript. Requires a reference genome in Biostrings DNAStringSet format
export.novel.transcript.sequences <- function(novel.gtf.file, reference.fasta){
  novel.gtf.data <- rtracklayer::import(novel.gtf.file)
  for(transcript.id in unique(novel.gtf.data$transcript_id)){
    exon.gtf <- GenomicRanges::GRangesList(novel.gtf.data[novel.gtf.data$transcript_id==transcript.id & novel.gtf.data$type=="exon",])
    if(any(strand(exon.gtf)=="*")){
      cat("Cannot get transcript sequence with ambiguous strand\n")
      next
    }
    trascript.seqs <- GenomicFeatures::extractTranscriptSeqs(reference.fasta, exon.gtf)
    names(trascript.seqs) <- paste0(novel.gtf.file, transcript.id)
    writeXStringSet(trascript.seqs, filepath = paste0(novel.gtf.file, ".", transcript.id, ".fa"))
  }
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

#### Process data ####

novel.gtf.files <- list.files(path = "data/stringtie",pattern = "*.gtf$")
novel.gtf.genomes <- sapply(novel.gtf.files, \(x) {
  species <- str_extract(x, "^\\w+")
  
  genome.data <- get.genome.data() %>%
    as.data.frame %>%
    dplyr::filter(CommonName==species) %>%
    dplyr::select(GTF)
  genome.file <- gsub("\\.112.gtf", ".dna.toplevel.fa.gz", genome.data$GTF)
} )

# TODO: read, plot and export fasta for all transcript files

#### Manual analyses #####

cat("Reading genome FASTA\n")
chicken.genome <- read.reference.genome("genomes/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz")

files <- list.files(path = "data/stringtie",pattern = "*.gtf$", full.names = TRUE)

gtf.data <- do.call(rbind, lapply(files, read.novel.transcripts))
plots <- lapply(files,plot.novel.transcripts )
combined.plot <- patchwork::wrap_plots(plots, ncol = 3)+ patchwork::plot_layout(guides = "collect", axes = "collect", axis_titles = "collect")
ggsave("data/stringtie/chicken.transcripts.all.png", combined.plot, units = "mm", height = 230, width = 170, dpi = 300)

# BLAST of the sequence hits https://www.ncbi.nlm.nih.gov/gene/121109097
# It's not in the Ensembl genebuild because it's a predicted transcript only, so the reads were aggregated into Zfx during feature mapping


# Check stranded junctions
chicken.plus <- read.rds.file("data/merged/chicken.testis.adult.ENSGALG00010003052.sense.Rds_+")
chicken.minus <- read.rds.file("data/merged/chicken.testis.adult.ENSGALG00010003052.sense.Rds_-")

patchwork::wrap_plots(list(make.sashimi.panel(chicken.plus, label="Sense")$plot, 
                make.gene.track.sashimi.panel(chicken.plus, label="Sense",show.x.axis=FALSE, is.collapse.introns=TRUE)$plot, 
                make.sashimi.panel(chicken.minus, label="Antisense")$plot, 
                make.gene.track.sashimi.panel(chicken.minus, label="Antisense",show.x.axis=FALSE, is.collapse.introns=TRUE)$plot), ncol = 2)+
  plot_annotation(title = "Chicken adult testis")+
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect")

ggsave("data/stringtie/chicken.transcripts.stranded.png", last_plot(), units = "mm", height = 170, width = 170, dpi = 300)


platy.plus <- read.rds.file("data/merged/platypus.testis.adult.ENSOANG00000046710.sense.Rds_+")
make.gene.track.sashimi.panel(platy.plus, label="Sense",show.x.axis=FALSE, is.collapse.introns=TRUE)$plot
