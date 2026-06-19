cat("Setup: Loading packages\n")
packages <- c(
  "parallel", "installr", "xlsx", "tidyverse", "GenomicRanges",
  "fs", "data.table", "patchwork", "grid", "scales", "ggbeeswarm",
  "rtracklayer", "Rsamtools", "bitops", "rlang", "R.utils"
)

suppressPackageStartupMessages({
  is.installed <- sapply(packages, require, character.only = TRUE)
  if (!all(is.installed)) stop("The following packages are required:", paste(packages[!is.installed], collapse = ", "))
})

cat("Setup: Defining global functions\n")

TIME.ORDER <- factor(c("birth", "mid-meiosis", "adult", "Day_00-06", "Day_07-13", "Day_14-20", "Day_21-27"),
  levels = c("birth", "mid-meiosis", "adult", "Day_00-06", "Day_07-13", "Day_14-20", "Day_21-27")
)

#### Common functions ####

save.double.width <- function(filename, plot, width = 170, height = 170) {
  ggsave(filename, plot, units = "mm", height = height, width = width, dpi = 300)
}

# Write the given data frame to an Excel file
create.xlsx <- function(data, file.name) {
  data <- as.data.frame(data) # ensure not a tibble

  oldOpt <- options()
  options(xlsx.date.format = "yyyy-mm-dd") # change date format
  wb <- xlsx::createWorkbook(type = "xlsx")
  sh <- xlsx::createSheet(wb)
  xlsx::addDataFrame(data, sh, row.names = F)

  # How many columns in the data frame? Convert to letters in base 26
  endColPart1 <- LETTERS[ncol(data) / 26]
  endColPart2 <- LETTERS[ncol(data) %% 26]
  endCol <- paste0(endColPart1, endColPart2)

  xlsx::addAutoFilter(sh, paste0("A1:", endCol, "1"))
  xlsx::createFreezePane(sh, 2, 2, 2, 2) # freeze top row and first column
  xlsx::autoSizeColumn(sh, 1:ncol(data))
  xlsx::saveWorkbook(wb, file = file.name)
  options(oldOpt)
}

#### Reading metadata files ####


# Get the identifiers for genes of interest
get.gene.locations <- function(genome.data) {
  # We don't want to constantly reload the GTFs if they have not changed.
  # But, if a GTF is updated, redo everything from scratch and write a new data
  # file
  gene.locations <- readr::read_csv("metadata/gene_locations.csv", show_col_types = FALSE)

  if (file.exists("./data/gene_coordinates.csv")) {
    cat("Setup: An existing gene coordinate file was found\n")
    existing.coordinates <- readr::read_csv("./data/gene_coordinates.csv", show_col_types = FALSE)

    has.raw.gtfs <- all(genome.data$GTF_FILE %in% existing.coordinates$GTF_FILE)
    has.saved.gtfs <- all(existing.coordinates$GTF_FILE %in% genome.data$GTF_FILE)
    has.raw.genes <- all(gene.locations$GeneId %in% existing.coordinates$GeneId)
    has.saved.genes <- all(existing.coordinates$GeneId %in% gene.locations$GeneId)

    if (!has.raw.gtfs) {
      cat(
        "Setup: Missing coordinates from a GTF file in ./genomes : ",
        paste(genome.data$GTF_FILE[!genome.data$GTF_FILE %in% existing.coordinates$GTF_FILE], collapse = ", "),
        "\n"
      )
    }
    if (!has.saved.gtfs) {
      cat(
        "Setup: Saved coordinates do not have a GTF file in ./genomes : ",
        paste(existing.coordinates$GTF_FILE[!existing.coordinates$GTF_FILE %in% genome.data$GTF_FILE], collapse = ", "),
        "\n"
      )
    }
    if (!has.raw.genes) {
      cat(
        "Setup: Missing coordinates from a gene in metadata/gene_locations.csv : ",
        paste(gene.locations$GeneId[!gene.locations$GeneId %in% existing.coordinates$GeneId], collapse = ", "),
        "\n"
      )
    }
    if (!has.saved.genes) {
      cat("Setup: Saved coordinates from a gene are not found in in metadata/gene_locations.csv\n")
      cat(
        "Setup: Saved coordinates from a gene are not found in in metadata/gene_locations.csv : ",
        paste(existing.coordinates$GeneId[!existing.coordinates$GeneId %in% gene.locations$GeneId], collapse = ", "),
        "\n"
      )
    }
    if (has.raw.gtfs & has.saved.gtfs & has.raw.genes & has.saved.genes) {
      return(existing.coordinates)
    }
  }

  # Find these ids in the relevant GTF file and extract location
  cat("Setup: Finding gene coordinates in GTF\n")
  # Identify the coordinates of a given gene id from GTF. Expand by size on each flank if desired
  get.gene.coordinates <- function(common.name, gtf.file, size = 1000) {
    cat("Setup: Reading GTF file for", common.name, "\n")
    if (!file.exists(gtf.file)) stop("Missing GTF file", gtf.file)
    gtf.data <- rtracklayer::import(gtf.file)

    # Get all genes for this gtf file
    species.zfxy.location.data <- gene.locations[gene.locations$CommonName == common.name, ]
    species.zfxy.location.data$GTF_FILE <- gtf.file
    species.zfxy.location.data$Location <- sapply(species.zfxy.location.data$GeneId, \(gene.id){
      filt.data <- gtf.data[gtf.data$gene_id == gene.id]
      # Expand on each flank
      # Return a location string
      location.string <- paste0(
        unique(GenomicRanges::seqnames(filt.data)), ":",
        min(GenomicRanges::start(filt.data) - size), "-",
        max(GenomicRanges::end(filt.data) + size)
      )
      cat("Setup: Found", gene.id, "in", common.name, "at", location.string, "\n")
      location.string
    })

    species.zfxy.location.data
  }

  # Bind each gene into a data frame
  locations <- do.call(rbind, mapply(get.gene.coordinates,
    common.name = genome.data$CommonName,
    gtf.file = genome.data$GTF_FILE,
    MoreArgs = list(
      size = 1000
    ), # ensure flanking lncRNAs will be detected
    SIMPLIFY = FALSE
  ))
  readr::write_csv(locations, file = "./data/gene_coordinates.csv")

  locations
}


# Annotate which exons contain interesting features for labelling plots
get.annotated.exons <- function() {
  readr::read_csv("metadata/gene_features.csv", show_col_types = FALSE)
}

# Read all selected samples from ./metadata
# i.e. all files with .filt. in the name
read.selected.samples <- function() {
  cat("Setup: Reading selected samples\n")
  # Make a factor of times to allow ordering of plots
  do.call(rbind, lapply(list.files(path = "metadata", pattern = "*.filt.csv", full.names = TRUE), read.csv)) %>%
    dplyr::mutate(Timepoint = factor(Timepoint, levels = TIME.ORDER)) %>%
    dplyr::arrange(CommonName, Run)
}


# Read the metadata to find samples. Aggregate to groups based on tissue type
# and note which samples still need processing
make.sample.groups <- function() {
  # Create command to merge bams in groups
  SELECTED.SAMPLES %>%
    dplyr::group_by(Organism, Organism_part, Timepoint, CommonName) %>% # not by sex - no difference seen in first pass
    dplyr::mutate(bam.file = paste0("data/", CommonName, "/", Run, ".bam")) %>%
    dplyr::summarise(
      bams = paste(bam.file, collapse = " "),
      all.bams.present = all(file.exists(bam.file)),
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
}

#### Global variables ####

cat("Setup: Defining global variables\n")

GENOME.DATA <- readr::read_csv("metadata/genomes.csv", show_col_types = FALSE)
# |>
# dplyr::mutate(
#   # Keep .gz extension for FASTA file - only used once uncompressed
#   FASTA_FILE = case_when(basename(FASTA_URL) == "unmasked.fa.gz" ~ paste0("./genomes/", Genome, ".fa.gz"),
#     .default = paste0("./genomes/", basename(FASTA_URL))
#   ),
#   # Keep GTF uncompressed, used several times in pipeline
#   GTF_FILE = case_when(basename(GTF_URL) == "genes.gtf.gz" ~ paste0("./genomes/", Genome, ".gtf"),
#     .default = paste0("./genomes/", stringr::str_remove(basename(GTF_URL), ".gz"))
#   )
# )

# Download GTF files if missing
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

# Read the filtered SRR samples and merge the genome and gene metadata
SELECTED.SAMPLES <- read.selected.samples()

cat("Setup: Common functions and global variables loaded\n")
