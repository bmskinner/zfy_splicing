# Make sashimi plot based on ggsashimi.py modified to write data objects to Rds
# when run This custom sashimi plot ensures the transcripts are always left to
# right irrespective of strand.
source("src/functions.R")
source("src/ggsashimi.R")

cat("Plot sashimi: running shashimi plotting\n")

#### Create individual plots  ####

# Ensure output dirs exist
fs::dir_create(c("report/species", "report/timepoints", "report/tissues", "report/raw_sashimi"))

# Read all GTF files once, since we have multiple genes/tissues per species
GTF.DATA <- read_gtf_data(GENOME.DATA$GTF_FILE, GENOME.DATA$CommonName)

cat("Plot sashimi: Making figures\n")

bam.files <- data.frame(path = list.files(path = "data/stringtie", pattern = ".*.bam$", full.names = TRUE)) |>
  dplyr::mutate(file = basename(path)) |>
  tidyr::separate_wider_delim(file,
    delim = ".", names = c("species", "tissue", "timepoint", "sex", "gene_id", "gene_name", "ext"),
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
  sex <- bam.row$sex
  gene_id <- bam.row$gene_id
  gene_name <- bam.row$gene_name

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

  # Skip completed files for testing
  final.out.file <- paste0("report/raw_sashimi/", paste(c(species, sex, tissue, timepoint, gene_id, gene_name), collapse = "."), ".expanded.png")
  if (file.exists(final.out.file)) next

  cat("Detecting splice junctions for", i, ": ", species, tissue, timepoint, gene_id, "\n")

  gene.data <- GENE.LOCATIONS[GENE.LOCATIONS$GeneId == gene_id & GENE.LOCATIONS$CommonName == species, ] # filter on species too - some genomes do not have an accession for geneid
  coords <- parse_coordinates(gene.data$Location)

  sashimi.data <- read_sashimi_data(
    bam.file = bam.row$path,
    gtf.data = GTF.DATA[[species]],
    chr = coords$coord.chr, start = coords$coord.start, end = coords$coord.end,
    reference.gene.id = gene_id,
    reference.transcript.id = gene.data$CanonicalTranscriptId
  )

  if (sashimi.data$total.reads == 0) next

  # Create plot with collapsed introns
  sashimi.plot.collapsed <- make_sashimi_coverage_plot(sashimi.data,
    is.collapse.introns = TRUE, show.x.axis = FALSE,
    min.spanning.reads = 2, label = paste0(species, "\n", sex, "\n", tissue, "\n", timepoint, "\n", gene_name)
  )

  species.aggregate[[species]] <- append(species.aggregate[[species]], sashimi.plot.collapsed$plot)
  tissue.aggregate[[tissue]] <- append(tissue.aggregate[[tissue]], sashimi.plot.collapsed$plot)
  timepoint.aggregate[[timepoint]] <- append(timepoint.aggregate[[timepoint]], sashimi.plot.collapsed$plot)

  save.double.width(
    paste0(
      "report/raw_sashimi/",
      paste(c(species, sex, tissue, timepoint, gene_id, gene_name), collapse = "."),
      ".condensed.png"
    ),
    sashimi.plot.collapsed$plot,
    height = 50
  )
  # Create plot with expanded introns
  sashimi.plot.expanded <- make_sashimi_coverage_plot(sashimi.data,
    is.collapse.introns = FALSE,
    min.spanning.reads = 2, paste0(species, "\n", sex, "\n", tissue, "\n", timepoint, "\n", gene_name)
  )

  save.double.width(
    paste0(
      "report/raw_sashimi/",
      paste(c(species, sex, tissue, timepoint, gene_id, gene_name), collapse = "."),
      ".expanded.png"
    ),
    sashimi.plot.expanded$plot,
    height = 50
  )
}

#### Create combined plots ####

# Combine plots for each species, developmental stage and timepoint
mapply(\(species, plot.list){
  cat("Plots for", species, "\n")
  # cat("There are ", length(plot.list), "plots\n")

  plots <- patchwork::wrap_plots(plot.list, ncol = 1)
  save.double.width(
    paste0("report/species/", species, ".png"),
    plots,
    height = 50 * length(plots)
  )
}, names(species.aggregate), species.aggregate)

cat("Plot sashimi: Done!\n")
