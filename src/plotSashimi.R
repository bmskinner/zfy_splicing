# Make sashimi plot based on ggsashimi.py modified to write data objects to Rds
# when run This custom sashimi plot ensures the transcripts are always left to
# right irrespective of strand.
source("src/functions.R")
source("src/ggsashimi.R")

cat("Plot sashimi: running shashimi plotting\n")

# Ensure output dirs exist
fs::dir_create(c("report/species", "report/timepoints", "report/tissues", "report/raw_sashimi"))

# Read all GTF files once, since we have multiple genes/tissues per species
GTF.DATA <- read_gtf_data(GENOME.DATA$GTF_FILE[c(1, 10)], GENOME.DATA$CommonName[c(1, 10)])

cat("Plot sashimi: Making figures\n")

bam.files <- data.frame(path = list.files(path = "data/stringtie", pattern = ".*.bam", full.names = TRUE)) |>
  dplyr::mutate(file = basename(path)) |>
  tidyr::separate_wider_delim(file,
    delim = ".", names = c("species", "tissue", "timepoint", "gene_id", "gene_name", "ext"),
    too_few = "debug", too_many = "debug"
  )

species.aggregate <- list()
tissue.aggregate <- list()
timepoint.aggregate <- list()

for (i in 1:nrow(bam.files)) {
  bam.row <- bam.files[i, ]
  species <- bam.row$species
  tissue <- bam.row$tissue
  timepoint <- bam.row$timepoint
  gene_id <- bam.row$gene_id

  # Create aggregates of all plots generated for later combination
  if (is.null(species.aggregate[[species]])) {
    species.aggregate[[species]] <- list()
  }

  if (is.null(tissue.aggregate[[tissue]])) {
    tissue.aggregate[[tissue]] <- list()
  }

  if (is.null(timepoint.aggregate[[timepoint]])) {
    timepoint.aggregate[[timepoint]] <- list()
  }

  cat("Detecting splice junctions for", species, tissue, timepoint, gene_id, "\n")

  gene.data <- GENE.LOCATIONS[GENE.LOCATIONS$GeneId == gene_id, ]
  coords <- parse_coordinates(gene.data$Location)

  sashimi.data <- read_sashimi_data(
    bam.file = bam.row$path,
    gtf.data = GTF.DATA[[species]],
    chr = coords$coord.chr, start = coords$coord.start, end = coords$coord.end,
    reference.gene.id = gene.data$GeneId,
    reference.transcript.id = gene.data$CanonicalTranscriptId
  )

  # Create plot with collapsed introns
  sashimi.plot.collapsed <- make_sashimi_plot(sashimi.data,
    is.collapse_introns = TRUE, show.x.axis = FALSE,
    min.spanning.reads = 5, label = paste0(species, "\n", tissue, "\n", timepoint)
  )

  save.double.width(paste0("report/raw_sashimi/", species, ".", tissue, ".", timepoint, ".", gene_id, ".condensed.png"),
    sashimi.plot.collapsed$plot,
    height = 80
  )
  # Create plot with expanded introns
  sashimi.plot.expanded <- make_sashimi_plot(sashimi.data,
    is.collapse_introns = FALSE,
    min.spanning.reads = 5, label = paste0(species, "\n", tissue, "\n", timepoint)
  )

  save.double.width(paste0("report/raw_sashimi/", species, ".", tissue, ".", timepoint, ".", gene_id, ".expanded.png"),
    sashimi.plot.expanded$plot,
    height = 80
  )

  species.aggregate[[species]] <- append(species.aggregate[[species]], sashimi.plot.collapsed$plot)
  tissue.aggregate[[tissue]] <- append(tissue.aggregate[[tissue]], sashimi.plot.collapsed$plot)
  timepoint.aggregate[[timepoint]] <- append(timepoint.aggregate[[timepoint]], sashimi.plot.collapsed$plot)
}

for (species in names(species.aggregate)) {
  png.file <- paste0("report/species/", species, ".png")
  wrapped.plots <- patchwork::wrap_plots(species.aggregate[[species]], ncol = 1)
  save.double.width(filename = png.file, plot = wrapped.plots, height = 50 * length(wrapped.plots))
}

for (tissue in names(tissue.aggregate)) {
  png.file <- paste0("report/tissues/", tissue, ".png")
  wrapped.plots <- patchwork::wrap_plots(tissue.aggregate[[tissue]], ncol = 1)
  save.double.width(filename = png.file, plot = wrapped.plots, height = 50 * length(wrapped.plots))
}

for (timepoint in names(timepoint.aggregate)) {
  png.file <- paste0("report/timepoints/", timepoint, ".png")
  wrapped.plots <- patchwork::wrap_plots(timepoint.aggregate[[timepoint]], ncol = 1)
  save.double.width(filename = png.file, plot = wrapped.plots, height = 50 * length(wrapped.plots))
}


cat("Plot sashimi: Done!\n")
