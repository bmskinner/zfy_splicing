#!/bin/Rscript
# Merge bam files for splice junction counting
# Invoke from the base directory
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(parallel)
library(fs)
source("src/functions.R")
# Read the filtered samples, match folder names
filtered.samples <- read.filtered.samples()

# Create command to merge bams in groups
groups <- filtered.samples %>% 
  dplyr::group_by(Organism, Organism_part, Timepoint, CommonName) %>% # not by sex - no difference seen in first pass
  dplyr::mutate(bam.file = paste0("data/", CommonName, "/", Run, ".bam"),
                lock.file = paste0("data/", CommonName, "/", Run, ".lck")) %>%
  dplyr::summarise(bams = paste(bam.file, collapse = " "),
                   all.bams.present = all(file.exists(bam.file)),
                   lock.files.exist  = any(file.exists(lock.file)), # bam may be in process of being written
                   count = n()) %>%
  dplyr::mutate(merged.bam = paste0("data/merged/", CommonName, ".", Organism_part, ".", Timepoint, ".bam"),
                samtools.merge.arguments = paste("merge -@ 7 -r -o", merged.bam,  bams))  %>%
  dplyr::mutate(merged.bam.exists = file.exists(merged.bam),
                index.exists = file.exists(paste0(merged.bam, ".csi"))) 

# Merge the bams
cat("Merging bams\n")

fs::dir_create("data/merged")
to.merge <- groups %>% # don't repeat merging
  dplyr::filter(all.bams.present & !merged.bam.exists & !lock.files.exist) # ensure we only try to merge when all bams of a group are available and complete
if(nrow(to.merge)>0){
  mapply(system2, command="samtools", args=to.merge$samtools.merge.arguments)
}

# Index the bams
cat("Indexing bams\n")
to.index <-  groups %>% # only index if the bam is present and there is no index
  dplyr::mutate(merged.bam.exists = file.exists(merged.bam)) %>% # update if merged bam exists
  dplyr::filter(merged.bam.exists & !index.exists)

# Index with CSI since opossum chromosomes are longer than the max for bai
if(nrow(to.index)>0){
  mapply(system2, command="samtools", args=paste("index -@ 7 -c ", to.index$merged.bam))
}

# Read GTF files from ./genomes
gtf.data <- read.gtf.data()

# Identify regions to select
select.range <- function(common.name, gene.id, size=1000){
  gtf <- gtf.data[[common.name]]

  data <- gtf[gtf$gene_id==gene.id]
  # Expand by size on each flank
  if(size>0) data <- GenomicRanges::resize(data, size*2, fix = "center")
  paste0(unique(GenomicRanges::seqnames(data)), ":", min(GenomicRanges::start(data)), "-", max(GenomicRanges::end(data)))
}

# Select gene locations from the GTF files
GENE.LOCATIONS$FlankedLocations <- mapply(select.range, common.name=GENE.LOCATIONS$CommonName, gene.id=GENE.LOCATIONS$EnsemblId, size=0, SIMPLIFY = TRUE)

# Combine the gene locations with bams
groups <- merge(groups, GENE.LOCATIONS, by="CommonName")

create.xlsx(groups, "report/grouped.bams.xlsx")

# Extract splice junctions for each combination
to.sashimi <-  groups %>% # don't repeat merging
  dplyr::mutate(merged.bam.exists = file.exists(merged.bam)) %>%
  dplyr::filter(merged.bam.exists) %>% # ensure we only try to sashimi when all bams of a group are available
  dplyr::mutate(junctions.file = paste0("data/merged/", CommonName, ".", Organism_part, ".", Timepoint, ".", EnsemblId, ".sense.Rds"))

if(nrow(to.sashimi)==0){
  stop("No valid samples to extract")
}

# ggsashimi requires a conda environment with pysam installed
cat("Extracting splice sites from", nrow(to.sashimi), "samples\n")
for(i in 1:nrow(to.sashimi)){
  data <- to.sashimi[i,]
  if(!file.exists(data$junctions.file)){
    cmd <-  paste0("activate ggsashimi && python src/ggsashimi.py --bam ", data$merged.bam, 
                   " --coordinates ", data$FlankedLocations, 
                   " --gtf genomes/", data$GTF, 
                   " --out-prefix ", data$junctions.file, 
                   " --strand SENSE  --out-format png")
    cat("source", cmd, "\n")
    system2("source",cmd)
  }
}

# Zip the results
fs::file_delete("values.tar.gz")
system2("tar", "-czf values.tar.gz data/merged/*.Rds*")
cat("Done!\n")

