# Visualise transcripts from StringTie
# This expects StringTie binary in ./bin
source("src/functions.R")
flog.info("Visualise transcripts: Beginning\n")
#### Define functions ####
# Given a gzipped fasta file, read and ensure seqnames are just the chr name
read.reference.genome <- function(fa.gz.file) {
  genome <- Biostrings::readDNAStringSet(fa.gz.file)
  # replace FA header with only the chr names
  # Anything up to the first space
  names(genome) <- str_extract(names(genome), "^[\\w\\d]+")
  genome
}

# Read gtf file and convert to tibble.
read.novel.transcripts <- function(novel.gtf.file) {
  novel.gtf <- rtracklayer::import(novel.gtf.file)
  as_tibble(novel.gtf) %>%
    dplyr::mutate(file = novel.gtf.file)
}

# Read gtf file and export FASTA of each transcript. Requires a reference genome in Biostrings DNAStringSet format
export.novel.transcript.sequences <- function(novel.gtf.file, reference.fasta) {
  novel.gtf.data <- rtracklayer::import(novel.gtf.file)
  for (transcript.id in unique(novel.gtf.data$transcript_id)) {
    exon.data <- novel.gtf.data[novel.gtf.data$transcript_id == transcript.id & novel.gtf.data$type == "exon", ]
    exon.gtf <- GenomicRanges::GRangesList(exon.data)

    # note that GenomicFeatures::extractTranscriptSeqs does not consider orientation
    # it uses the rank order of the exons
    # Hence we need to reverse the exons if they are on the - strand
    if (any(strand(exon.gtf) == "-")) {
      exon.gtf <- sort(exon.gtf, decreasing = TRUE)
    }

    if (any(strand(exon.gtf) == "*")) {
      cat("Cannot get transcript sequence with ambiguous strand for ", novel.gtf.file, transcript.id, "\n")
      next
    }

    transcript.seqs <- GenomicFeatures::extractTranscriptSeqs(reference.fasta, exon.gtf)
    names(transcript.seqs) <- paste0(novel.gtf.file, transcript.id)
    writeXStringSet(transcript.seqs, filepath = paste0(novel.gtf.file, ".", transcript.id, ".fa"))
  }
}

# Read novel transcript gtf file and plot.
plot.novel.transcripts <- function(gtf.file, stringtie.plot.file) {
  flog.info("Plotting", gtf.file, "\n")
  
  novel.gtf <- read.novel.transcripts(gtf.file) %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::arrange(strand, transcript_id, start)
  
  n.transcripts <- length(unique(novel.gtf$transcript_id))

  novel.gtf.exon <- novel.gtf %>% dplyr::filter(type == "exon")

  novel.gtf.intron <- novel.gtf %>%
    dplyr::group_by(transcript_id, strand) %>%
    dplyr::summarise(
      start = min(start),
      end = max(end),
      .groups = "drop_last"
    ) |>
    dplyr::mutate(strand = ifelse(as.character(strand) == "*", "-", as.character(strand)))


  if (!is.null(novel.gtf.exon$ref_gene_id)) {
    novel.gene.ids <- novel.gtf.exon %>%
      dplyr::group_by(transcript_id, reference_id) %>%
      dplyr::summarise(
        start = min(start),
        end = max(end),
        mid = (start + end) / 2,
        .groups = "drop_last"
      ) %>%
      na.omit()
  } else {
    novel.gene.ids <- data.frame()
  }

  transcript.plot <- novel.gtf %>%
    dplyr::filter(type == "exon") %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = transcript_id
    ))

  # Only plot if we have introns - ignore for single exon genes
  if (nrow(novel.gtf.intron) > 0) {
    transcript.plot <- transcript.plot +
      geom_intron(
        data = novel.gtf.intron,
        aes(strand = strand), col = "black"
      )
  }

  transcript.plot <- transcript.plot +
    geom_range(aes(fill = as.numeric(cov))) +
    scale_fill_viridis_c() +
    labs(
      x = "Position", y = "Assembled transcript",
      title = stringtie.gtf, fill = "TPM"
    ) +
    geom_text(
      data = add_exon_number(novel.gtf.exon, "transcript_id"),
      aes(
        x = (start + end) / 2, # plot label at midpoint of exon
        label = exon_number
      ),
      nudge_y = 0.4
    )

  if (nrow(novel.gene.ids) > 0) {
    # Add transcript id if present
    transcript.plot <- transcript.plot +
      geom_text(
        data = novel.gene.ids,
        aes(
          x = mid,
          y = transcript_id,
          label = reference_id
        ),
        nudge_y = -0.4
      )
  }

  transcript.plot <- transcript.plot +
    theme_bw() +
    theme(
      title = element_text(),
      legend.position = "top"
    )

  save.plot(stringtie.plot.file, transcript.plot,
    width = 170, height = min(1000, 20 + (n.transcripts * 50))
  )

  TRUE
}

#### Find and read merged StringTie GTF files ####

# Find the transcript GTFs

stringtie.gtf.files <- SELECTED.SAMPLES |>
  merge(GENE.LOCATIONS, by=c("CommonName", "GTF_FILE", "GeneId", "Location")) |>
  dplyr::filter(!(Group %in% c("RBMX", "RBMY"))) |>
  dplyr::select(Organism, Tissue, Timepoint, CommonName, Sex, GTF_FILE, GeneId, Gene) |>
  dplyr::distinct() |>
  dplyr::mutate(
    stringtie.gtf = paste0("data/merged/", CommonName, ".", Tissue, ".", Timepoint, ".", Sex,".", GeneId,".gtf"),
    stringtie.gtf.exists = file.exists(stringtie.gtf),
    stringtie.plot.file = paste0("report/stringtie/", CommonName, ".", Tissue, ".", Timepoint, ".", Sex,".", Gene,".png")
  )

stringtie.to.plot <- stringtie.gtf.files |>
  dplyr::filter(stringtie.gtf.exists)

invisible(mcmapply(\(x, ...) tryCatch(plot.novel.transcripts(x,...), 
                                      error = \(e) {
                                        print(e)
                                        return(FALSE)
                                      }), 
                   stringtie.to.plot$stringtie.gtf, 
                   stringtie.to.plot$stringtie.plot.file,
                   mc.cores = DEFAULT.MC.CORES))

#### Export FASTA sequence of novel transcripts ####

# Read the genome FASTA files
# cat("Reading genome FASTA files\n")
# genome.fastas <- lapply(GENOME.DATA$FASTA_FILE, read.reference.genome)
# names(genome.fastas) <- GENOME.DATA$CommonName

# # Extract the transcript sequence for each
# cat("Extracting transcript FASTA sequences\n")
# sapply(novel.gtf.files, \(x){
#   common.name <- stringr::str_split(basename(x), "\\.")[[1]][1]
#   cat("Extracting from genome", common.name, "\n")
#   fasta <- genome.fastas[[common.name]]
#   tryCatch(export.novel.transcript.sequences(x, fasta), error = \(e) print(e))
# })

#### Manual analyses #####
# read.novel.transcripts("data/stringtie/rat.testis.adult.ENSRNOG00000053042.Zfy2.gtf")


# cat("Reading genome FASTA\n")
# chicken.genome <- read.reference.genome("genomes/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz")
# export.novel.transcript.sequences("data/stringtie/chicken.testis.adult.ENSGALG00010003052.gtf", chicken.genome)
flog.info("Visualise transcripts: Done\n")
