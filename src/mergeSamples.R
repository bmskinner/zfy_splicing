#!/bin/Rscript
# Merge bam files for splice junction counting. This aggregates the data plotted
# from each sample individually.
cat("Merge samples: Beginning\n")
source("src/functions.R")

cat("Merge samples: selecting bams for merging\n")

fs::dir_create("data/merged")

#### Create merge command for each gene, tissue, species, timepoint, sex ####

# If a merged bam does not exist, create it. If a merged bam exists, and is more
# recent than any of the individual bams, don't overwrite it.

# Ditto for StringTie - run if we need to

groups <- SELECTED.SAMPLES |>
  merge(GENE.LOCATIONS, by=c("CommonName", "GTF_FILE", "GeneId")) |>
  dplyr::group_by(Organism, Tissue, Timepoint, CommonName, Sex, GTF_FILE, GeneId, Gene) |>
  dplyr::mutate(
    bam.file = paste0("data/", CommonName, "/", Run, ".", GeneId, ".bam"),
    lock.file = paste0("data/", CommonName, "/", Run, ".lck"),
    bam.mtime = file.mtime(bam.file)
  ) |>
  dplyr::summarise(
    bams = paste(bam.file, collapse = " "),
    bams.exist = all(file.exists(bam.file)),
    lock.files.exist = any(file.exists(lock.file)), # bam may be in process of being written
    count = n(),
    max.bam.mtime = max(bam.mtime), # what is the most recent bam creation time?
    .groups = "drop_last"
  ) |>
  dplyr::mutate(
    merged.bam = paste0("data/merged/", CommonName, ".", Tissue, ".", Timepoint, ".", Sex,".", GeneId,".bam"),
    merged.bam.exists = file.exists(merged.bam),
    is.replace.merged.bam = merged.bam.exists & file.mtime(merged.bam) < max.bam.mtime,
    is.merge  = !lock.files.exist & bams.exist & (!merged.bam.exists | is.replace.merged.bam), # Do not overwrite if no bams have changed
    samtools.merge.arguments = paste("merge -@ 7 -r -o", merged.bam, bams),
    merged.gtf = stringr::str_replace(merged.bam, ".bam", ".gtf"),
    merged.gtf.exists = file.exists(merged.gtf),
    is.replace.merged.gtf = merged.gtf.exists & file.mtime(merged.gtf) < max.bam.mtime,
    is.stringtie  = !lock.files.exist & bams.exist & (!merged.gtf.exists | is.replace.merged.gtf),
    GFF.FILE = stringr::str_replace(GTF_FILE, ".gtf", ".gff"),
    annotation.file = ifelse(file.exists(GFF.FILE), GFF.FILE, GTF_FILE), # only some replaced 
    stringtie.arguments = paste0("-o ", merged.gtf, " -p 1 -l ", CommonName," -G ",annotation.file, " -f 0.01 ", merged.bam),
  )

#### Merge the bams ####

# Any existing merged files will be overwritten
to.merge <- groups |>
  dplyr::filter(is.merge) # ensure we only try to merge when all bams of a group are available and complete

if (nrow(to.merge) > 0) {
  cat("Merge samples: Merging bams\n")
  mcmapply(system2, command = "samtools", 
           args = to.merge$samtools.merge.arguments, 
           mc.cores = DEFAULT.MC.CORES)
} 

#### Index the bams ####
cat("Merge samples: Indexing bams\n")

# Index with CSI since opossum chromosomes are longer than the max for bai
invisible(mcmapply(system2, command = "samtools", 
         args = paste("index -@ 7 -c ", to.merge$merged.bam),
         mc.cores = DEFAULT.MC.CORES))

#### Run StringTie on the merged bam if needed ####
cat("Merge samples: Running stringtie on merged bams\n")
to.stringtie <- groups |>
  dplyr::filter(is.stringtie) 

invisible(mcmapply(system2, command = "stringtie", 
         args = to.stringtie$stringtie.arguments,
         mc.cores = DEFAULT.MC.CORES))

#### Make a summary table of which bams were merged ####
create.xlsx(groups, "report/merged.bams.xlsx")
cat("Merge samples: Exported merged bam table to report/merged.bams.xlsx\n")

#### Zip the results ####
cat("Merge samples: Done!\n")
