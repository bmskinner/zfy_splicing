.packages <- c(
  "parallel", "installr", "openxlsx2", "tidyverse", "GenomicRanges",
  "fs", "data.table", "patchwork", "grid", "scales", "ggbeeswarm",
  "rtracklayer", "Rsamtools", "bitops", "rlang", "R.utils", "gridExtra", "png",
  "Biostrings", "GenomicFeatures", "BSgenome", "ggtranscript", # devtools::install_github("dzhang32/ggtranscript")
  "futile.logger"
)

suppressPackageStartupMessages({
  is.installed <- sapply(.packages, require, character.only = TRUE)
  if (!all(is.installed)) stop("The following packages are required:", paste(.packages[!is.installed], collapse = ", "))
  rm(is.installed)
})

invisible(flog.layout(layout.format('~t ~m')))

flog.info("Setup: Defining global functions and variables\n")

# The file that will contain samples to be processed. Created de novo on each run
# from selectSamples.R
MAPPING.FILE <- "data/mapping.samples.csv"
CDS.MAPPING.FILE <- "cds/cds.samples.csv"

# Use multiple cores on Unix via the parallel package
DEFAULT.MC.CORES <- ifelse(installr::is.windows(), 1, 6)

# Chart display order for timepoints
TIME.ORDER <- factor(c("birth", "adolescence", "adult", "Day_00-06", "Day_07-13", "Day_14-20", "Day_21-27"),
  levels = c("birth", "adolescence", "adult", "Day_00-06", "Day_07-13", "Day_14-20", "Day_21-27")
)

#### Common functions ####

#' Save a ggplot with print resolution
#'
#' @param filename 
#' @param plot 
#' @param ... 
#'
#' @returns
#' @export
#'
#' @examples
save.plot <- function(filename, plot, ...) {
  ggsave(filename, plot, units = "mm", dpi = 300, , create.dir = TRUE, ...)
}

#' Save a ggplot as a double column width image at 300 dpi.
#'
#' @param filename where to save. Missing directories will be created.
#' @param plot the ggplot to save
#' @param width width in mm
#' @param height height in mm
#'
#' @returns
#' @export
#'
#' @examples
save.double.width <- function(filename, plot, width = 170, height = 170) {
  ggsave(filename, plot, units = "mm", height = height, width = width, dpi = 300, create.dir = TRUE)
}


# add the given data frame to an Excel workbook
add.and.freeze <- function(workbook, data, sheet.name) {
  oldOpt <- options()
  options("openxlsx2.dateFormat" = "yyyy-mm-dd")
  options("openxlsx2.datetimeFormat" = "yyyy-mm-dd hh:mm:ss")
  options("openxlsx2.na" = "_openxlsx_NULL")
  options("openxlsx2.maxWidth" = 40)

  wb <- workbook |>
    openxlsx2::wb_add_worksheet(sheet.name) |>
    openxlsx2::wb_add_data(sheet.name, data, row_names = FALSE) |>
    openxlsx2::wb_freeze_pane(sheet = sheet.name, first_row = TRUE, first_col = TRUE) |>
    openxlsx2::wb_add_filter(sheet = sheet.name, rows = 1, cols = 1:ncol(data)) |>
    openxlsx2::wb_set_col_widths(sheet = sheet.name, cols = 1:ncol(data), widths = "auto")

  options(oldOpt)
  wb
}

# Create an Excel file with column filtering
# data - the date frame to export
# file.name - the name of the file to export to
create.xlsx <- function(data, file.name) {
  data <- as.data.frame(data) # ensure not a tibble

  oldOpt <- options()
  options("openxlsx2.dateFormat" = "yyyy-mm-dd")
  options("openxlsx2.datetimeFormat" = "yyyy-mm-dd hh:mm:ss")
  options("openxlsx2.na" = "_openxlsx_NULL")

  wb <- openxlsx2::wb_workbook() |>
    add.and.freeze(data, "Sheet 1")
  # openxlsx2::wb_add_worksheet("Sheet 1") |>
  # openxlsx2::wb_add_data("Sheet 1", data, row_names = FALSE) |>
  # openxlsx2::wb_freeze_pane(first_row = TRUE, first_col = TRUE) |>
  # openxlsx2::wb_add_filter(rows = 1, cols = 1:ncol(data)) |>
  # openxlsx2::wb_set_col_widths(cols = 1:ncol(data), widths = "auto")

  openxlsx2::wb_save(wb, file = file.name)
  options(oldOpt)
}

#### Reading metadata files ####


# Read the file of gene and transcript ids we want to study.
# Find their location in the appropriate genome GTF.
# Allows us to update genome assembly version without manually updating coordinates.
get.gene.locations <- function(genome.data) {
  # We don't want to constantly reload the GTFs if they have not changed.
  # But, if a GTF is updated, redo everything from scratch and write a new data
  # file.
  gene.ids <- readr::read_csv("metadata/gene_ids.csv", show_col_types = FALSE)

  if (file.exists("./data/gene_coordinates.csv")) {
    flog.info("Setup: An existing gene coordinate file was found\n")
    existing.coordinates <- readr::read_csv("./data/gene_coordinates.csv", show_col_types = FALSE)

    has.raw.gtfs <- all(genome.data$GTF_FILE %in% existing.coordinates$GTF_FILE)
    has.saved.gtfs <- all(existing.coordinates$GTF_FILE %in% genome.data$GTF_FILE)
    has.raw.genes <- all(gene.ids$GeneId %in% existing.coordinates$GeneId)
    has.saved.genes <- all(existing.coordinates$GeneId %in% gene.ids$GeneId)

    if (!has.raw.gtfs) {
      flog.info(
        "Setup: Missing coordinates from a GTF file in ./genomes : ",
        paste(genome.data$GTF_FILE[!genome.data$GTF_FILE %in% existing.coordinates$GTF_FILE], collapse = ", "),
        "\n"
      )
    }
    if (!has.saved.gtfs) {
      flog.info(
        "Setup: Saved coordinates do not have a GTF file in ./genomes : ",
        paste(existing.coordinates$GTF_FILE[!existing.coordinates$GTF_FILE %in% genome.data$GTF_FILE], collapse = ", "),
        "\n"
      )
    }
    if (!has.raw.genes) {
      flog.info(
        "Setup: Missing coordinates from a gene in metadata/gene_locations.csv : ",
        paste(gene.ids$GeneId[!gene.ids$GeneId %in% existing.coordinates$GeneId], collapse = ", "),
        "\n"
      )
    }
    if (!has.saved.genes) {
      flog.info("Setup: Saved coordinates from a gene are not found in in metadata/gene_locations.csv\n")
      flog.info(
        "Setup: Saved coordinates from a gene are not found in in metadata/gene_locations.csv : ",
        paste(existing.coordinates$GeneId[!existing.coordinates$GeneId %in% gene.ids$GeneId], collapse = ", "),
        "\n"
      )
    }
    if (has.raw.gtfs & has.saved.gtfs & has.raw.genes & has.saved.genes) {
      return(existing.coordinates)
    }
  }

  # Find these ids in the relevant GTF file and extract location
  flog.info("Setup: Finding gene coordinates in GTF\n")
  # Identify the coordinates of a given gene id from GTF. Expand by size on each flank if desired
  get.gene.coordinates <- function(common.name, gtf.file, size = 1000) {
    flog.info("Setup: Reading GTF file for", common.name, "\n")
    if (!file.exists(gtf.file)) stop("Missing GTF file", gtf.file)
    gtf.data <- rtracklayer::import(gtf.file)

    # Get all genes for this gtf file
    species.ids <- gene.ids[gene.ids$CommonName == common.name, ]
    species.ids$GTF_FILE <- gtf.file


    species.ids$Location <- sapply(species.ids$GeneId, \(gene.id){
      # The transcript ids in the GTF
      filt.data <- gtf.data[gtf.data$gene_id == gene.id]
      filt.transcripts <- unique(filt.data[filt.data$type == "transcript", ])

      # The reference transcript id for the gene we are checking
      transcript.id <- species.ids |>
        dplyr::filter(GeneId == gene.id) |>
        dplyr::select(CanonicalTranscriptId) |>
        dplyr::pull()


      # Sanity check that the transcript id is in the GTF
      if (!any((transcript.id %in% filt.transcripts$transcript_id))) {
        stop(
          "Unable to find a transcript with id ", transcript.id, " for gene ",
          gene.id, " in ", common.name, ". The found transcripts were: ",
          paste(filt.transcripts$transcript_id, collapse = ", ")
        )
      }

      # Expand on each flank
      # Return a location string
      location.string <- paste0(
        unique(GenomicRanges::seqnames(filt.data)), ":",
        min(GenomicRanges::start(filt.data) - size), "-",
        max(GenomicRanges::end(filt.data) + size)
      )
      flog.info("Setup: Found", gene.id, "in", common.name, "at", location.string, "\n")
      location.string
    })

    species.ids
  }

  # Bind each gene into a data frame
  locations <- do.call(rbind, parallel::mcmapply(get.gene.coordinates,
    common.name = genome.data$CommonName,
    gtf.file = genome.data$GTF_FILE,
    MoreArgs = list(
      size = 1000
    ), # ensure flanking lncRNAs will be detected
    SIMPLIFY = FALSE,
    mc.cores = DEFAULT.MC.CORES
  ))
  readr::write_csv(locations, file = "./data/gene_coordinates.csv")

  locations
}


# Annotate which exons contain interesting features for labelling plots
get.annotated.exons <- function() {
  readr::read_csv("metadata/gene_features.csv", show_col_types = FALSE)
}

# Read the metadata to find samples. Aggregate to groups based on tissue type
# and note which samples still need processing
make.sample.groups <- function() {
  # Create command to merge bams in groups
  SELECTED.SAMPLES %>%
    dplyr::group_by(Organism, Tissue, Timepoint, Sex, CommonName) %>%
    dplyr::mutate(bam.file = paste0("data/", CommonName, "/", Run, ".bam")) %>%
    dplyr::summarise(
      bams = paste(bam.file, collapse = " "),
      all.bams.present = all(file.exists(bam.file)),
      count = n(),
      .groups = "drop_last"
    ) %>%
    dplyr::mutate(
      merged.bam = paste0("data/merged/", CommonName, ".", Tissue, ".", Timepoint, ".", Sex, ".bam"),
      samtools.merge.arguments = paste("merge -@ 7 -r -o", merged.bam, bams)
    ) %>%
    dplyr::mutate(
      merged.bam.exists = file.exists(merged.bam),
      index.exists = file.exists(paste0(merged.bam, ".csi"))
    )
}

#### Global variables ####

GENOME.DATA <- readr::read_csv("metadata/genomes.csv", show_col_types = FALSE) |>
  dplyr::mutate(Clade = as.factor(Clade))

# Download GTF files if missing so we can find gene coordinates
download.gtf <- function(file, url) {
  if (!file.exists(file)) {
    gz.file <- paste0(file, ".gz")
    cat("Downloading ", url, "to ", gz.file, "\n")
    download.file(url, destfile = gz.file)
    dl.file <- R.utils::gunzip(gz.file)
    fs::file_move(dl.file, file)
  }
}
invisible(mapply(download.gtf, GENOME.DATA$GTF_FILE, GENOME.DATA$GTF_URL))

# Match the gene ids to coordinates in the genome version downloaded
GENE.LOCATIONS <- get.gene.locations(GENOME.DATA)

# Filter SRR samples and merge the genome metadata
source("src/selectSamples.R")

# Global object with samples being analysed
SELECTED.SAMPLES <- read.csv(MAPPING.FILE) |>
  dplyr::mutate(Timepoint = factor(Timepoint, levels = TIME.ORDER)) |>
  dplyr::arrange(CommonName, Run)

# Locations of exon junctions for exon 2 splice detection
JUNCTION.COORDINATES <- read.csv("metadata/exon_junctions.csv")

flog.info("Setup: Common functions and global variables loaded\n")
