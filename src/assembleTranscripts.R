# "Assemble transcripts usign StringTie
# This expects StringTie binary in ./bin
cat("Assemble transcripts: Beginning\n")
source("src/functions.R")

invisible(fs::dir_create("data/stringtie"))
if (file.exists("data/stringtie/ratios.txt")) file.remove("data/stringtie/ratios.txt")

# Get the distinct groups
cat("Assemble transcripts: Reading sample groups\n")
sample.groups <- merge(make.sample.groups(), GENE.LOCATIONS, by = "CommonName")

#### "Assemble transcripts ####


assemble.transcript <- function(common.name, tissue, timepoint, sex, gene.id, gene, coordinates, full.bam.file, full.gtf.file) {
  cat("Assemble transcripts: Running Stringtie on", common.name, tissue, timepoint, sex, gene.id, "\n")

  gene.bam.file <- paste0("data/stringtie/", common.name, ".", tissue, ".", timepoint, ".", sex, ".", gene.id, ".", gene, ".bam")

  if (!file.exists(gene.bam.file)) {
    # Write the reads covering the gene
    system2("samtools", paste0("view -o ", gene.bam.file, " ", full.bam.file, " '", coordinates, "'")) # keep paste0 for location quoting
    system2("samtools", paste("index -c ", gene.bam.file))
  } else {
    cat("Assemble transcripts: Gene bam already exists for", common.name, tissue, timepoint, gene.id, "\n")
  }

  # Count the reads on sense and antisense strands
  system2("bash", paste("src/countStrandRatio.sh", gene.bam.file))

  gtf.out.file <- paste0("data/stringtie/", common.name, ".", tissue, ".", timepoint, ".", sex, ".", gene.id, ".", gene, ".gtf")

  if (!file.exists(gtf.out.file)) {
    # Run stringtie using the reference genome to guide assembly
    system2("bin/stringtie", paste(
      "-o ", gtf.out.file, # output file name
      "-p 1 -l", common.name, # label for novel transcripts
      "-G", full.gtf.file, # genome annotation
      "-f 0.01", # min fraction of reads supporting splices
      gene.bam.file
    )) # input file to analyse
  } else {
    cat("Assemble transcripts: Stringtie gtf already exists for", common.name, tissue, timepoint, sex, gene.id, "\n")
  }
}

cat("Assemble transcripts: Assembling transcripts with Stringtie from", nrow(sample.groups), "sample groups\n")
invisible(mapply(
  assemble.transcript,
  sample.groups$CommonName,
  sample.groups$Organism_part,
  sample.groups$Timepoint,
  sample.groups$sex,
  sample.groups$GeneId,
  sample.groups$Gene,
  sample.groups$Location,
  sample.groups$merged.bam,
  sample.groups$GTF_FILE
))

#### Assess strandedness of reads ####
cat("Assemble transcripts: Reading strand ratios\n")

if (!file.exists("data/stringtie/ratios.txt")) {
  stop("Assemble transcripts: Strand ratio file from Stringtie does not exist: data/stringtie/ratios.txt")
}

read.ratios <- read.delim("data/stringtie/ratios.txt", sep = " ", header = FALSE)
colnames(read.ratios) <- c("sample", "forward", "reverse", "total")

cat("Assemble transcripts: Parsing strand ratios\n")
read.ratios <- tidyr::separate_wider_delim(read.ratios, sample,
  delim = ".",
  names = c("species", "tissue", "timepoint", "sex", "gene.id", "gene"),
  too_many = "debug"
) %>%
  dplyr::mutate(
    total.stranded = forward + reverse,
    f.stranded = total.stranded / total,
    f.forward = forward / total.stranded
  ) %>%
  dplyr::filter(timepoint %in% c("adult", "mid-meiosis", "birth"))

create.xlsx(read.ratios, "report/strand_ratios.xlsx")


ggplot(read.ratios) +
  ggplot2::annotate("rect", xmax = Inf, xmin = -Inf, ymax = 0.1, ymin = -Inf, fill = "lightgreen") +
  ggplot2::annotate("rect", xmax = Inf, xmin = -Inf, ymax = Inf, ymin = 0.9, fill = "lightgreen") +
  ggplot2::geom_point(aes(x = species, y = f.forward, col = f.forward <= 0.1 | f.forward >= 0.9)) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "lightgrey")) +
  ggplot2::labs(y = "Strand ratio") +
  ggplot2::facet_grid(tissue ~ timepoint) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

save.double.width("report/strandedness.png", last_plot(), height = 230)

ggplot(read.ratios) +
  geom_point(aes(x = species, y = f.stranded)) +
  labs(y = "Stranded read fraction of total reads") +
  facet_grid(tissue ~ timepoint) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

save.double.width("report/stranded_fraction.png", last_plot(), height = 230)
cat("Assemble transcripts: Done!\n")
