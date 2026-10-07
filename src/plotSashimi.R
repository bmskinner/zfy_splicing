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
fs::dir_create(c(
  "report/species", "report/timepoints", "report/tissues",
  "report/raw_sashimi", "report/junctions", "report/merged_sashimi"
))

# Read all GTF files once, since we have multiple genes/tissues per species
GTF.DATA <- read_gtf_data(GENOME.DATA$GTF_FILE, GENOME.DATA$CommonName)

cat("Plot sashimi: Making figures\n")

# Read all single sample bam files and bind in the complete metadata
bam.files <- data.frame(path = list.files(path = "data", pattern = "[SDE]RR\\d+\\..*\\.bam$", 
                                          full.names = TRUE, recursive=TRUE)) |>
  dplyr::mutate(file = basename(path),
                mtime = file.mtime(path)) |>
  tidyr::separate_wider_delim(file,
                              delim = ".", names = c("Run", "GeneId", "ext")
  ) |>
  merge(GENE.LOCATIONS, by = "GeneId") |>
  merge(GENOME.DATA, by=c("CommonName", "GTF_FILE")) |>
  merge(SELECTED.SAMPLES, by = c("CommonName", "Run",  "Genome", "GTF_FILE", "GeneId", "Location")) |>
  dplyr::filter( !(Group %in% c("RBMX", "RBMY"))) # skip genes we don't need splice data from

create.sashimi.plot <- function(i){

  # Error handling
  skip.file <- FALSE
  
  tryCatch({
    
    bam.row <- bam.files[i, ]
    run <- bam.row$Run
    species <- bam.row$CommonName
    tissue <- bam.row$Tissue
    timepoint <- bam.row$Timepoint
    sex <- bam.row$Sex
    gene_id <- bam.row$GeneId
    gene_name <- bam.row$Gene
    bam.mtime <- bam.row$mtime
    
    # Skip completed files for testing
    final.out.file <- paste0("report/raw_sashimi/", paste(c(run, gene_id), collapse = "."), ".condensed.png")
    final.junction.file <- paste0(
      "report/junctions/",
      paste(c(run, gene_id), collapse = "."),
      ".junctions.csv"
    )
    
    # If the output image file exists and was created after the bam file (in
    # case a sample was remapped), we can skip
    if (file.exists(final.out.file) & file.exists(final.junction.file) & file.mtime(final.out.file)>bam.mtime) return(TRUE)
    
    gene.data <- GENE.LOCATIONS[GENE.LOCATIONS$GeneId == gene_id & GENE.LOCATIONS$CommonName == species, ] # filter on species too - some genomes do not have an accession for geneid
    coords <- parse_coordinates(gene.data$Location)
    group <- gene.data$Group
    
    # Skip missing data or genes we don't need splice data from
    if (length(group) == 0) return(TRUE)

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
      fs::file_touch(final.junction.file)
      fs::file_touch(final.out.file)
      return(TRUE)
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
    
  }, error=function(e){
    cat("Error making shashimi plot\n")
    print(e)
    # tryCatch will return an error object which will quit the loop by default
    # Set a skip variable instead
    skip.file <<- TRUE
  })
  
  return(skip.file)
  # Skip to next loop iteration if an error was caught
  # if(skip.file) { next }     
}

mclapply(1:nrow(bam.files), create.sashimi.plot, mc.cores = DEFAULT.MC.CORES)

#### Create aggregate plots for merged samples ####

merged.bam.files <- data.frame(path = list.files(path = "data/merged", pattern = ".*.bam$", 
                                          full.names = TRUE, recursive=TRUE)) |>
  dplyr::mutate(file = basename(path),
                mtime = file.mtime(path)) |>
  tidyr::separate_wider_delim(file,
                              delim = ".", names = c("CommonName", "Tissue", "Timepoint", "Sex",  "GeneId", "ext")
  ) |>
  merge(GENE.LOCATIONS, by =c("CommonName",  "GeneId")) |>
  merge(GENOME.DATA, by=c("CommonName", "GTF_FILE")) |>
  dplyr::filter(Group %in% c("ZFX", "ZFY"))

for (i in 1:nrow(merged.bam.files)) {
  
  # Error handling
  skip.file <- FALSE
  
  tryCatch({
    bam.row <- merged.bam.files[i, ]
    species <- bam.row$CommonName
    tissue <- bam.row$Tissue
    timepoint <- bam.row$Timepoint
    sex <- bam.row$Sex
    gene_id <- bam.row$GeneId
    gene_name <- bam.row$Gene
    bam.mtime <- bam.row$mtime
    
    gene.data <- GENE.LOCATIONS[GENE.LOCATIONS$GeneId == gene_id & GENE.LOCATIONS$CommonName == species, ] # filter on species too - some genomes do not have an accession for geneid
    coords <- parse_coordinates(gene.data$Location)
    group <- gene.data$Group
    if (length(group) == 0) next
    cat("Detecting splice junctions for", i, ": ", species, gene_id, "in group", group, "\n")
    
    final.out.file <- paste0("report/merged_sashimi/", paste(c(species,tissue, 
                                                               timepoint, sex, gene_id, gene_name), 
                                                             collapse = "."), ".condensed.png")
    
    if(file.exists(final.out.file) & file.mtime(final.out.file)>bam.mtime){ next }
    
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
    
    sashimi.plot.collapsed <- make_sashimi_coverage_plot(sashimi.data,
                                                         is.collapse.introns = TRUE, show.x.axis = FALSE,
                                                         min.spanning.reads = 2, label = paste0(species, "\n", sex, "\n", tissue, "\n", timepoint, "\n", gene_name)
    )
    
    save.double.width(final.out.file,
                      sashimi.plot.collapsed$plot,
                      height = 50
    )
    
    
  }, error=function(e){
    cat("Error making shashimi plot\n")
    print(e)
    # tryCatch will return an error object which will quit the loop by default
    # Set a skip variable instead
    skip.file <<- TRUE
  })
  
  
  # Skip to next loop iteration if an error was caught
  if(skip.file) { next }     
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
