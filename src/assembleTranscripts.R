# Assemble transcripts usign StringTie
# This expects StringTie binary in ./bin

library(fs)
library(ggplot2)
source("src/functions.R")


fs::dir_create("data/stringtie")
if (file.exists("data/stringtie/ratios.txt")) file.remove("data/stringtie/ratios.txt")

# Get the distinct groups
cat("Reading sample groups\n")
sample.groups <- merge(make.sample.groups(), GENE.LOCATIONS, by = "CommonName")

#### Assemble transcripts ####

assemble.transcript <- function(common.name, tissue, timepoint, gene.id, gene, coordinates, full.bam.file, full.gtf.file) {
  cat("Assembling transcripts from", common.name, tissue, timepoint, gene.id, "\n")

  gene.bam.file <- paste0("data/stringtie/", common.name, ".", tissue, ".", timepoint, ".", gene.id, ".", gene, ".bam")
  gtf.out.file <- paste0("data/stringtie/", common.name, ".", tissue, ".", timepoint, ".", gene.id, ".", gene, ".gtf")

  # Write the reads covering the gene
  system2("samtools", paste0("view -o ", gene.bam.file, " ", full.bam.file, " '", coordinates, "'"))

  # Run stringtie usin the reference genome to guide assembly
  system2("bin/stringtie", paste(
    "-o ", gtf.out.file, # output file name
    "-p 1 -l", common.name, # label for novel transcripts
    "-G", paste0("genomes/", full.gtf.file), # genome annotation
    "-f 0.01", # min fraction of reads supporting splices
    gene.bam.file
  )) # input file to analyse

  # Count the reads on sense and antisense strands
  system2("bash", paste("src/countStrandRatio.sh", gene.bam.file))
}
cat("Assembling transcripts from", nrow(sample.groups), "sample groups\n")
mapply(
  assemble.transcript,
  sample.groups$CommonName,
  sample.groups$Organism_part,
  sample.groups$Timepoint,
  sample.groups$GeneId,
  sample.groups$Gene,
  sample.groups$FlankedLocations,
  sample.groups$merged.bam,
  sample.groups$GTF
)

#### Assess strandedness of reads ####
cat("Reading strand ratios\n")

read.ratios <- read.delim("data/stringtie/ratios.txt", sep = " ", header = FALSE)
colnames(read.ratios) <- c("sample", "forward", "reverse", "total")
read.ratios <- tidyr::separate_wider_delim(read.ratios, sample,
  delim = ".",
  names = c("species", "tissue", "timepoint", "gene.id", "gene"),
  too_many = "debug"
) %>%
  dplyr::mutate(
    total.stranded = forward + reverse,
    f.stranded = total.stranded / total,
    f.forward = forward / total.stranded
  ) %>%
  dplyr::filter(timepoint %in% c("adult", "mid-meiosis", "birth"))


ggplot(read.ratios) +
  annotate("rect", xmax = Inf, xmin = -Inf, ymax = 0.1, ymin = -Inf, fill = "lightgreen") +
  annotate("rect", xmax = Inf, xmin = -Inf, ymax = Inf, ymin = 0.9, fill = "lightgreen") +
  geom_point(aes(x = species, y = f.forward)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "Strand ratio") +
  facet_grid(tissue ~ timepoint) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

save.double.width("report/strandedness.png", last_plot())

ggplot(read.ratios) +
  geom_point(aes(x = species, y = f.stranded)) +
  labs(y = "Stranded read fraction of total reads") +
  facet_grid(tissue ~ timepoint) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

save.double.width("report/stranded_fraction.png", last_plot())
