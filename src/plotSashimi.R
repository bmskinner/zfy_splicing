# Make sashimi plot based on ggsashimi.py modified to write data objects to Rds
# when run This custom sashimi plot ensures the transcripts are always left to
# right irrespective of strand.
source("src/functions.R")
source("src/ggsashimi.R")
library(gridExtra)
library(png)
library(grid)
cat("Plot sashimi: running shashimi plotting\n")

#### Read the bed files of pair spanning reads ####

bed.files <- list.files(path = "data", pattern = "[SDE]RR.*.bed$",
                        full.names = TRUE, recursive=TRUE)

# Read the bed file. Keep only reads with pairs on opposite strands of the
# same chromosome.
read.bed <- function(file){
  if(file.size(file)==0) return()
  readr::read_tsv(file, col_names = c("Seqname", "Mate1Start", "Mate1End", "Seqname2",
                               "Mate2Start", "Mate2End", "ReadName", "Score",
                               "Mate1Strand", "Mate2Strand"),
           col_types = c("ciiciicicc")) |>
    dplyr::filter(Mate1Strand!=Mate2Strand,
                  Seqname==Seqname2,
                  Mate1Start!=Mate2Start,
                  Mate1End != Mate2End) |>
    dplyr::mutate(File = file,
                  Run = str_extract(File, "([SDE]RR\\d+)", group = 1),
                  GeneId = str_extract(File, "[SDE]RR\\d+\\.([\\w\\d]+)\\.bam.bed$", group = 1)
    ) |>
    dplyr::select(-Seqname2, -File, -Score)
}

bed.data <- do.call(dplyr::bind_rows, lapply(bed.files, read.bed) ) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    Mate1Size = Mate1End - Mate1Start + 1,
    Mate2Size = Mate2End - Mate2Start + 1,
    InsertStart = min(Mate1Start, Mate2Start),
    InsertEnd = max(Mate1End, Mate2End),
    InsertSize = InsertEnd - InsertStart + 1,
    JunctionStart = min(Mate1End, Mate2End),
    JunctionEnd = max(Mate1Start, Mate2Start)
  ) |>
  merge(SELECTED.SAMPLES, by = c("Run", "GeneId")) |>
  merge(GENE.LOCATIONS, by = c("CommonName", "GeneId", "GTF_FILE", "Location"))


pair.spanning.data <- bed.data |>
  dplyr::filter(Group %in% c("ZFX", "ZFY")) |>
  dplyr::group_by(Run)|>
  dplyr::mutate(RunMedianInsertSize = median(InsertSize)) |>
  merge(JUNCTION.COORDINATES, by = c("CommonName", "GeneId")) |>
  dplyr::mutate(SpansExon2 = JunctionStart <= start & JunctionEnd >= end) |>
  dplyr::filter(SpansExon2)
 


#### Create individual plots  ####

# Ensure output dirs exist
fs::dir_create(c(
  "report/species", "report/timepoints", "report/tissues",
  "report/raw_sashimi", "report/junctions", "report/merged_sashimi"
))

# Read all GTF files once, since we have multiple genes/tissues per species
GTF.DATA <- read_gtf_data(GENOME.DATA$GTF_FILE, GENOME.DATA$CommonName)

cat("Plot sashimi: Making figures\n")

# Read all single sample bam files and bind in the complete metadata
bam.files <- data.frame(path = list.files(path = "data", pattern = "[SDE]RR.*.bam$", 
                                          full.names = TRUE, recursive=TRUE)) |>
  dplyr::mutate(file = basename(path)) |>
  tidyr::separate_wider_delim(file,
                              delim = ".", names = c("Run", "GeneId", "ext")
  ) |>
  merge(GENE.LOCATIONS, by = "GeneId") |>
  merge(GENOME.DATA, by=c("CommonName", "GTF_FILE")) |>
  merge(SELECTED.SAMPLES, by = c("CommonName", "Run",  "Genome", "GTF_FILE", "GeneId", "Location")) |>
  dplyr::filter( !(Group %in% c("RBMX", "RBMY"))) # skip genes we don't need splice data from

for (i in 1:nrow(bam.files)) {
  bam.row <- bam.files[i, ]
  run <- bam.row$Run
  species <- bam.row$CommonName
  tissue <- bam.row$Tissue
  timepoint <- bam.row$Timepoint
  sex <- bam.row$Sex
  gene_id <- bam.row$GeneId
  gene_name <- bam.row$Gene

  gene.data <- GENE.LOCATIONS[GENE.LOCATIONS$GeneId == gene_id & GENE.LOCATIONS$CommonName == species, ] # filter on species too - some genomes do not have an accession for geneid
  coords <- parse_coordinates(gene.data$Location)
  group <- gene.data$Group
  
  # Skip missing data or genes we don't need splice data from
  if (length(group) == 0) next
  
  # Skip completed files for testing
  final.out.file <- paste0("report/raw_sashimi/", paste(c(run, gene_id), collapse = "."), ".condensed.png")
  final.junction.file <- paste0(
    "report/junctions/",
    paste(c(run, gene_id), collapse = "."),
    ".junctions.csv"
  )
  
  if (file.exists(final.out.file) & file.exists(final.junction.file)) next
  
  # Only read the bam file if needed

  cat("Detecting splice junctions for", i, ": ", run, gene_id, "in group", group, "\n")

  sashimi.data <- read_sashimi_data(
    bam.file = bam.row$path,
    gtf.data = GTF.DATA[[species]],
    chr = coords$coord.chr, start = coords$coord.start, end = coords$coord.end,
    reference.gene.id = gene_id,
    reference.transcript.id = gene.data$CanonicalTranscriptId
  )

  if (sashimi.data$total.reads == 0){
    fs::file_touch(final.pair.spanning.file)
    fs::file_touch(final.junction.file)
    fs::file_touch(final.out.file)
    next
  } 

  junction.data <- sashimi.data$junctions |>
    dplyr::mutate(Run = run, GeneId = gene_id)

  readr::write_csv(junction.data,
    file = final.junction.file,
    quote = "needed"
  )
  
  #### Create plot with collapsed introns ####
  
  if(!file.exists(final.out.file) ){

    sashimi.plot.collapsed <- make_sashimi_coverage_plot(sashimi.data,
                                                         is.collapse.introns = TRUE, show.x.axis = FALSE,
                                                         min.spanning.reads = 1, label = paste0(run,"\n", species, "\n", sex, "\n", tissue, "\n", timepoint, "\n", gene_name)
    )
    
    save.double.width(
      paste0(
        "report/raw_sashimi/",
        paste(c(run, gene_id), collapse = "."),
        ".condensed.png"
      ),
      sashimi.plot.collapsed$plot,
      height = 50
    )
  }
}

#### Create aggregate plots for merged samples ####

merged.bam.files <- data.frame(path = list.files(path = "data/merged", pattern = ".*.bam$", 
                                          full.names = TRUE, recursive=TRUE)) |>
  dplyr::mutate(file = basename(path)) |>
  tidyr::separate_wider_delim(file,
                              delim = ".", names = c("CommonName", "Tissue", "Timepoint", "Sex",  "Gene", "ext")
  ) |>
  merge(GENE.LOCATIONS, by =c("CommonName",  "Gene")) |>
  merge(GENOME.DATA, by=c("CommonName", "GTF_FILE")) |>
  dplyr::filter(Group %in% c("ZFX", "ZFY"))

for (i in 1:nrow(merged.bam.files)) {
  bam.row <- merged.bam.files[i, ]
  species <- bam.row$CommonName
  tissue <- bam.row$Tissue
  timepoint <- bam.row$Timepoint
  sex <- bam.row$Sex
  gene_id <- bam.row$GeneId
  gene_name <- bam.row$Gene
  
  gene.data <- GENE.LOCATIONS[GENE.LOCATIONS$GeneId == gene_id & GENE.LOCATIONS$CommonName == species, ] # filter on species too - some genomes do not have an accession for geneid
  coords <- parse_coordinates(gene.data$Location)
  group <- gene.data$Group
  if (length(group) == 0) next
  cat("Detecting splice junctions for", i, ": ", species, gene_id, "in group", group, "\n")
  
  final.out.file <- paste0("report/merged_sashimi/", paste(c(species,tissue, 
                                                             timepoint, sex, gene_id, gene_name), 
                                                           collapse = "."), ".condensed.png")
  
  sashimi.data <- read_sashimi_data(
    bam.file = bam.row$path,
    gtf.data = GTF.DATA[[species]],
    chr = coords$coord.chr, start = coords$coord.start, end = coords$coord.end,
    reference.gene.id = gene_id,
    reference.transcript.id = gene.data$CanonicalTranscriptId
  )
  
  if (sashimi.data$total.reads == 0){
    fs::file_touch(final.out.file)
    next
  } 
  
  if(!file.exists(final.out.file) ){
    
    sashimi.plot.collapsed <- make_sashimi_coverage_plot(sashimi.data,
                                                         is.collapse.introns = TRUE, show.x.axis = FALSE,
                                                         min.spanning.reads = 2, label = paste0(species, "\n", sex, "\n", tissue, "\n", timepoint, "\n", gene_name)
    )
    
    save.double.width(final.out.file,
      sashimi.plot.collapsed$plot,
      height = 50
    )
  }
  
}

#### Match the junction coordinates found with the coding exon 2 splice sites ####

# Some junction file may be zero size, exclude these
junction.files <- list.files(path = "report/junctions", 
                             pattern = "[EDS]RR.*junctions.csv", full.names = TRUE)
junction.files <- junction.files[file.size(junction.files)>50] # 50 bytes is the col headers only

# Filters all splice junctions to those in e1, e2, e3
junction.data <- do.call(bind_rows, lapply(junction.files, read.csv)) |>
  merge(JUNCTION.COORDINATES, by = c("GeneId", "start", "end")) |>
  dplyr::select(-type) |>
  tidyr::pivot_wider(id_cols = c(CommonName, Run, GeneId),  
                     names_from = Junction, values_from = count, values_fill = 0) |>
  merge(SELECTED.SAMPLES, by=c("CommonName", "Run"))|>
  dplyr::group_by(CommonName, Sex, Tissue, Timepoint, GeneId)|>
  dplyr::summarise(E1E2 = sum(E1E2),
                   E2E3 = sum(E2E3),
                   E1E3 = sum(E1E3),
                   .groups = "drop_last")|>
  dplyr::mutate(MeanCanonical = (E1E2 + E2E3)/2,
                pctE2Spliced = E1E3 / MeanCanonical * 100)

# Save for combination with gene expression levels in featureCounts analysis
readr::write_tsv(junction.data, "report/coding_exon_splice_junctions.tsv")

cat("Plot sashimi: Done!\n")
