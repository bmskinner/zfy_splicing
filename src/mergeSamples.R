#!/bin/Rscript
# Merge bam files for splice junction counting
# Invoke from the base directory
cat("Merge samples: Beginning\n")
source("src/functions.R")

cat("Merge samples: selecting bams for merging\n")

#### Create command to merge bams in groups ####
groups <- SELECTED.SAMPLES %>%
  dplyr::group_by(Organism, Organism_part, Timepoint, CommonName) %>% # not by sex - no difference seen in first pass
  dplyr::mutate(
    bam.file = paste0("data/", CommonName, "/", Run, ".bam"),
    lock.file = paste0("data/", CommonName, "/", Run, ".lck")
  ) %>%
  dplyr::summarise(
    bams = paste(bam.file, collapse = " "),
    all.bams.present = all(file.exists(bam.file)),
    lock.files.exist = any(file.exists(lock.file)), # bam may be in process of being written
    count = n(),
    .groups = "drop_last"
  ) %>%
  dplyr::mutate(
    merged.bam = paste0("data/merged/", CommonName, ".", Organism_part, ".", Timepoint, ".bam"),
    samtools.merge.arguments = paste("merge -@ 7 -r -o", merged.bam, bams)
  ) %>%
  dplyr::mutate(
    merged.bam.exists = file.exists(merged.bam),
    index.exists = file.exists(paste0(merged.bam, ".csi"))
  )

#### Merge the bams ####

fs::dir_create("data/merged")
to.merge <- groups %>% # don't repeat merging
  dplyr::filter(all.bams.present & !merged.bam.exists & !lock.files.exist) # ensure we only try to merge when all bams of a group are available and complete
if (nrow(to.merge) > 0) {
  cat("Merge samples: Merging bams\n")
  mapply(system2, command = "samtools", args = to.merge$samtools.merge.arguments)
} else {
  cat("Merge samples: No unmerged bams\n")
}

#### Index the bams ####
cat("Merge samples: Indexing bams\n")
to.index <- groups %>% # only index if the bam is present and there is no index
  dplyr::mutate(merged.bam.exists = file.exists(merged.bam)) %>% # update if merged bam exists
  dplyr::filter(merged.bam.exists & !index.exists)

# Index with CSI since opossum chromosomes are longer than the max for bai
if (nrow(to.index) > 0) {
  mapply(system2, command = "samtools", args = paste("index -@ 7 -c ", to.index$merged.bam))
} else {
  cat("Merge samples: No unindexed bams\n")
}

# Combine the gene locations with bams
groups <- merge(groups, GENE.LOCATIONS, by = "CommonName")

create.xlsx(groups, "report/grouped.bams.xlsx")
cat("Merge samples: Exported grouped bam table to report/grouped.bams.xlsx\n")

#### Zip the results ####
cat("Merge samples: Done!\n")
