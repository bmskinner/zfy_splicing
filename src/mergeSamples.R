#!/bin/Rscript
# Merge bam files for splice junction counting. This aggregates the data plotted
# from each sample individually.
cat("Merge samples: Beginning\n")
source("src/functions.R")

cat("Merge samples: selecting bams for merging\n")

fs::dir_delete("data/merged")
fs::dir_create("data/merged")

#### Create merge command for each gene, tissue, species, timepoint, sex ####
groups <- SELECTED.SAMPLES |>
  merge(GENE.LOCATIONS, by=c("CommonName", "GTF_FILE", "GeneId")) |>
  dplyr::group_by(Organism, Tissue, Timepoint, CommonName, Sex, GeneId, Gene) |>
  dplyr::mutate(
    bam.file = paste0("data/", CommonName, "/", Run, ".", Gene, ".bam"),
    lock.file = paste0("data/", CommonName, "/", Run, ".lck")
  ) |>
  dplyr::summarise(
    bams = paste(bam.file, collapse = " "),
    lock.files.exist = any(file.exists(lock.file)), # bam may be in process of being written
    count = n(),
    .groups = "drop_last"
  ) |>
  dplyr::mutate(
    merged.bam = paste0("data/merged/", CommonName, ".", Tissue, ".", Timepoint, ".", Sex,".", GeneId,".bam"),
    samtools.merge.arguments = paste("merge -@ 7 -r -o", merged.bam, bams)
  )

#### Merge the bams ####

# Any existing merged files will be overwritten
to.merge <- groups |>
  dplyr::filter(!lock.files.exist) # ensure we only try to merge when all bams of a group are available and complete

if (nrow(to.merge) > 0) {
  cat("Merge samples: Merging bams\n")
  mapply(system2, command = "samtools", args = to.merge$samtools.merge.arguments)
} 

#### Index the bams ####
cat("Merge samples: Indexing bams\n")

# Index with CSI since opossum chromosomes are longer than the max for bai
mapply(system2, command = "samtools", args = paste("index -@ 7 -c ", to.merge$merged.bam))


# Make a summary table of which bams were merged
create.xlsx(groups, "report/merged.bams.xlsx")
cat("Merge samples: Exported merged bam table to report/grouped.bams.xlsx\n")

#### Zip the results ####
cat("Merge samples: Done!\n")
