#!/bin/Rscript
# Merge bam files for splice junction counting. This aggregates the data plotted
# from each sample individually.
cat("Merge samples: Beginning\n")
source("src/functions.R")

cat("Merge samples: selecting bams for merging\n")

fs::dir_create("data/merged")

#### Create merge command for each gene, tissue, species, timepoint, sex ####
groups <- SELECTED.SAMPLES |>
  merge(GENE.LOCATIONS, by=c("CommonName", "GTF_FILE")) |>
  dplyr::group_by(Organism, Tissue, Timepoint, CommonName, Sex, GeneId, Gene) |>
  dplyr::mutate(
    bam.file = paste0("data/", CommonName, "/", Run, ".", Gene, ".bam"),
    lock.file = paste0("data/", CommonName, "/", Run, ".lck")
  ) |>
  dplyr::summarise(
    bams = paste(bam.file, collapse = " "),
    all.bams.present = all(file.exists(bam.file)),
    lock.files.exist = any(file.exists(lock.file)), # bam may be in process of being written
    count = n(),
    .groups = "drop_last"
  ) |>
  dplyr::mutate(
    merged.bam = paste0("data/merged/", CommonName, ".", Tissue, ".", Timepoint, ".", Sex,".", Gene,".bam"),
    samtools.merge.arguments = paste("merge -@ 7 -r -o", merged.bam, bams)
  ) |>
  dplyr::mutate(
    merged.bam.exists = file.exists(merged.bam),
    index.exists = file.exists(paste0(merged.bam, ".csi"))
  )

#### Merge the bams ####

to.merge <- groups |>
  dplyr::filter(all.bams.present & !merged.bam.exists & !lock.files.exist) # ensure we only try to merge when all bams of a group are available and complete

if (nrow(to.merge) > 0) {
  cat("Merge samples: Merging bams\n")
  mapply(system2, command = "samtools", args = to.merge$samtools.merge.arguments)
} else {
  cat("Merge samples: No unmerged bams\n")
}

#### Index the bams ####
cat("Merge samples: Indexing bams\n")
to.index <- groups |> # only index if the bam is present and there is no index
  dplyr::mutate(merged.bam.exists = file.exists(merged.bam)) |> # update if merged bam exists
  dplyr::filter(merged.bam.exists & !index.exists)

# Index with CSI since opossum chromosomes are longer than the max for bai
if (nrow(to.index) > 0) {
  mapply(system2, command = "samtools", args = paste("index -@ 7 -c ", to.index$merged.bam))
} else {
  cat("Merge samples: No unindexed bams\n")
}

# Make a summary table of which bams were merged
create.xlsx(groups, "report/merged.bams.xlsx")
cat("Merge samples: Exported merged bam table to report/grouped.bams.xlsx\n")

#### Zip the results ####
cat("Merge samples: Done!\n")
