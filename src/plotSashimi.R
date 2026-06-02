# Make sashimi plot based on ggsashimi.py modified to write data objects to Rds
# when run This custom sashimi plot ensures the transcripts are always left to
# right irrespective of strand.
source("src/functions.R")
source("src/ggsashimi.R")

cat("Plot sashimi: running shashimi plotting\n")

# Ensure output dirs exist
fs::dir_create(c("report/species", "report/timepoints", "report/tissues", "report/raw_sashimi"))

# Read all GTF files once, since we have multiple genes/tissues per species
GTF.DATA <- read_gtf_data(GENOME.DATA$GTF_FILE, GENOME.DATA$CommonName)

cat("Plot sashimi: Making figures\n")

bam.files <- data.frame(path = list.files(path = "data/stringtie", pattern = ".*.bam", full.names = TRUE)) |>
  dplyr::mutate(file = basename(path)) |>
  tidyr::separate_wider_delim(file,
    delim = ".", names = c("species", "tissue", "timepoint", "gene_id", "gene_name", "ext"),
    too_few = "debug", too_many = "debug"
  )

for (i in 1:nrow(bam.files)) {
  bam.row <- bam.files[i, ]
  species <- bam.row$species
  tissue <- bam.row$tissue
  timepoint <- bam.row$timepoint
  gene_id <- bam.row$gene_id
  gene_name <- bam.row$gene_name

  cat("Detecting splice junctions for", species, tissue, timepoint, gene_id, "\n")

  gene.data <- GENE.LOCATIONS[GENE.LOCATIONS$GeneId == gene_id, ]
  coords <- parse_coordinates(gene.data$Location)

  tryCatch(
    {
      sashimi.data <- read_sashimi_data(
        bam.file = bam.row$path,
        gtf.data = GTF.DATA[[species]],
        chr = coords$coord.chr, start = coords$coord.start, end = coords$coord.end,
        reference.gene.id = gene.data$GeneId,
        reference.transcript.id = gene.data$CanonicalTranscriptId
      )

      # Create plot with collapsed introns
      sashimi.plot.collapsed <- make_sashimi_plot(sashimi.data,
        is.collapse.introns = TRUE, show.x.axis = FALSE,
        min.spanning.reads = 5, label = paste0(species, "\n", tissue, "\n", timepoint, "\n", gene_name)
      )

      save.double.width(paste0("report/raw_sashimi/", species, ".", tissue, ".", timepoint, ".", gene_id, ".", gene_name, ".condensed.png"),
        sashimi.plot.collapsed$plot,
        height = 80
      )
      # Create plot with expanded introns
      sashimi.plot.expanded <- make_sashimi_plot(sashimi.data,
        is.collapse.introns = FALSE,
        min.spanning.reads = 5, label = paste0(species, "\n", tissue, "\n", timepoint, "\n", gene_name)
      )

      save.double.width(paste0("report/raw_sashimi/", species, ".", tissue, ".", timepoint, ".", gene_id, ".", gene_name, ".expanded.png"),
        sashimi.plot.expanded$plot,
        height = 80
      )

      # Create a plot with coverage info. Only expanded available
      coverage.plot <- make_coverage_plot(sashimi.data, label = paste0(species, "\n", tissue, "\n", timepoint, "\n", gene_name))
      save.double.width(paste0("report/raw_sashimi/", species, ".", tissue, ".", timepoint, ".", gene_id, ".", gene_name, ".coverage.png"),
        coverage.plot,
        height = 80
      )
    },
    error = \(e) {
      cat("Error plotting sashimi data from", bam.row$path, "\n", paste(e))
      e
    }
  )
}

cat("Plot sashimi: Done!\n")
