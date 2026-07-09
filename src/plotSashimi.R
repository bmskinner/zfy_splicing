# Make sashimi plot based on ggsashimi.py modified to write data objects to Rds
# when run This custom sashimi plot ensures the transcripts are always left to
# right irrespective of strand.
source("src/functions.R")
source("src/ggsashimi.R")
library(gridExtra)
library(png)
library(grid)
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

for (i in 1:nrow(bam.files)) {
  bam.row <- bam.files[i, ]
  species <- bam.row$species
  tissue <- bam.row$tissue
  timepoint <- bam.row$timepoint
  sex <- bam.row$sex
  gene_id <- bam.row$gene_id
  gene_name <- bam.row$gene_name

  gene.data <- GENE.LOCATIONS[GENE.LOCATIONS$GeneId == gene_id & GENE.LOCATIONS$CommonName == species, ] # filter on species too - some genomes do not have an accession for geneid
  coords <- parse_coordinates(gene.data$Location)
  group <- gene.data$Group

  # Skip missing data or genes we don't need splice data from
  if (length(group) == 0) next
  if (group %in% c("RBMX", "RBMY")) next

  # Skip completed files for testing
  final.out.file <- paste0("report/raw_sashimi/", paste(c(species, sex, tissue, timepoint, gene_id, gene_name, group), collapse = "."), ".condensed.png")
  final.junction.file <- paste0(
    "report/junctions/",
    paste(c(species, sex, tissue, timepoint, gene_id, gene_name, group), collapse = "."),
    ".junctions.csv"
  )

  if (file.exists(final.out.file) & file.exists(final.junction.file)) next

  cat("Detecting splice junctions for", i, ": ", species, tissue, timepoint, gene_id, "in group", group, "\n")

  sashimi.data <- read_sashimi_data(
    bam.file = bam.row$path,
    gtf.data = GTF.DATA[[species]],
    chr = coords$coord.chr, start = coords$coord.start, end = coords$coord.end,
    reference.gene.id = gene_id,
    reference.transcript.id = gene.data$CanonicalTranscriptId
  )

  if (sashimi.data$total.reads == 0) next

  junction.data <- sashimi.data$junctions |>
    dplyr::mutate(
      CommonName = species, Sex = sex, Tissue = tissue, Timepoint = timepoint, GeneId = gene_id,
      Gene = gene_name, Group = group
    )

  readr::write_csv(junction.data,
    file = final.junction.file,
    quote = "needed"
  )

  if (file.exists(final.out.file)) next

  # Create plot with collapsed introns
  sashimi.plot.collapsed <- make_sashimi_coverage_plot(sashimi.data,
    is.collapse.introns = TRUE, show.x.axis = FALSE,
    min.spanning.reads = 2, label = paste0(species, "\n", sex, "\n", tissue, "\n", timepoint, "\n", gene_name)
  )

  save.double.width(
    paste0(
      "report/raw_sashimi/",
      paste(c(species, sex, tissue, timepoint, gene_id, gene_name, group), collapse = "."),
      ".condensed.png"
    ),
    sashimi.plot.collapsed$plot,
    height = 50
  )
  # Create plot with expanded introns
  # sashimi.plot.expanded <- make_sashimi_coverage_plot(sashimi.data,
  #   is.collapse.introns = FALSE,
  #   min.spanning.reads = 2, paste0(species, "\n", sex, "\n", tissue, "\n", timepoint, "\n", gene_name)
  # )
  #
  # save.double.width(
  #   paste0(
  #     "report/raw_sashimi/",
  #     paste(c(species, sex, tissue, timepoint, gene_id, gene_name), collapse = "."),
  #     ".expanded.png"
  #   ),
  #   sashimi.plot.expanded$plot,
  #   height = 50
  # )
}

#### Match the junction coordinates found with the coding exon 2 splice sites ####

# Coding exon 2 splice location manually determined from the sashimi plots and
# junction location file. Look for any matching junctions in each gene. This
# must be updated if genome assembly version changes.
junction.coordinates <- tribble(
  ~CommonName, ~GeneId, ~start, ~end,
  "anole", "ENSACAG00000007227", 125817818, 125821463,
  "cattle", "ENSBTAG00000007730", 119118397, 119133682,
  "cattle", "ENSBTAG00000078450", 7993697, 8010836,
  "chicken", "ENSGALG00010003052", 118304535, 118313129,
  "duck", "ENSAPLG00000011089", 116945669, 116953843,
  "echidna", "ZFX", 2887543, 2888139,
  "human", "ENSG00000005889", 24172800, 24207326,
  "human", "ENSG00000067646", 2953997, 2975095,
  "koala", "ENSPCIG00000024018", 12169119, 12204163,
  "macaque", "ENSMMUG00000009801", 24000889, 24035078,
  "macaque", "ENSMMUG00000046378", 258566, 281083,
  "mouse", "ENSMUSG00000079509", 93125905, 93145897,
  "mouse", "ENSMUSG00000053211", 735168, 759850,
  "mouse", "ENSMUSG00000000103", 2117212, 2133385,
  "opossum", "ZFX", 46374628, 46418348,
  "pig", "ENSSSCG00000021322", 20243490, 20266919,
  "pig", "ENSSSCG00000012179", 9911038, 9928932,
  "platypus", "ENSOANG00000046710", 15377677, 15398814,
  "rat", "ENSRNOG00000005624", 62806207, 62824310,
  "rat", "ENSRNOG00000053042", 156128, 165382,
  "tasmaniandevil", "ENSSHAG00000004983", 269241175, 269277566,
  "turkey", "ZFX", 116183643, 116192451,
  "wallaby", "ZFX", 452798276, 452836092,
  "xenopus", "ENSXETG00000007785", 42980939, 42986352,
  "zebrafinch", "ENSTGUG00000007219", 100739793, 100749515,
  "zebrafish", "ENSDARG00000074453", 17007505, 17011774
)

# Only coding exon 2 splice junctions are retained here
junction.data <- do.call(rbind, lapply(list.files(path = "report/junctions", pattern = ".*.csv", full.names = TRUE), read.csv)) |>
  merge(junction.coordinates, by = c("CommonName", "GeneId", "start", "end"))

# Save for combination with gene expression levels in featureCounts analysis
readr::write_tsv(junction.data, "report/coding_exon_2_splice_junctions.tsv")

#### Create combined plots ####

tissues <- unique(bam.files$tissue)
timepoints <- unique(bam.files$timepoint)
groups <- unique(GENE.LOCATIONS$Group)


for (tissue in tissues) {
  for (timepoint in timepoints) {
    for (group in groups) {
      plots <- lapply(
        list.files(path = "report/raw_sashimi", pattern = paste0(".*", tissue, ".*", timepoint, ".*", group, ".*png"), full.names = TRUE),
        \(x){
          img <- as.raster(readPNG(x))
          rasterGrob(img, interpolate = FALSE, default.units = "mm", width = 170, height = 50)
        }
      )
      ggsave(
        paste0("report/tissues/", tissue, ".", timepoint, ".pdf"),
        arrangeGrob(grobs = plots, nrow = length(plots), ncol = 1),
        dpi = 300, units = "mm", width = 170, height = 50 * length(plots)
      )
    }
  }
}



cat("Plot sashimi: Done!\n")
