# Make sashimi plot based on ggsashimi output.
# ggsashimi.py modified to write data objects to Rds when run
# This custom sashimi plot ensures the transcripts are always left to right
# irrespective of strand
cat("Plot sashimi: Beginning\n")
source("src/functions.R")
source("src/ggsashimi.R")

cat("Plot sashimi: running shashimi plotting\n")

# ANNOTATED.EXONS <- get.annotated.exons()

# Ensure output dirs exist
fs::dir_create(c("report/species", "report/timepoints", "report/tissues", "report/raw_sashimi"))

# Read all the GTF files to a global variable
read.gtf.data <- function() {
  cat("Plot sashimi: Reading full genome GTF files\n")
  gtf.data <- mclapply(GENOME.DATA$GTF_FILE, \(x) as.data.frame(rtracklayer::import(x)),
    mc.cores = ifelse(installr::is.windows(), 1, 6)
  )
  names(gtf.data) <- dplyr::pull(GENOME.DATA[, "CommonName"])
  cat("Plot sashimi: Read full genome GTF files\n")
  gtf.data
}

GTF.DATA <- read.gtf.data()

#### Main functions ####


#' Test if the given transcript is on the forward or reverse strand
#'
#' @param transcript.id the transcript to test
#' @param gtf.data the GTF data for the genome as GenomicRanges e.g. as read by rtracklayer
#'
#' @returns true if any exons of the given transcript are on the reverse strand, false otherwise
#' @export
#'
#' @examples
transcript.is.reverse.strand <- function(transcript.id, gtf.data) {
  exons <- gtf.data[gtf.data$type == "exon" & gtf.data$transcript_id == transcript.id, ]
  any(exons$strand == "-")
}

#
#' Given a canonical transcript id, find the splice junctions
#'
#' @param transcript.id the transcript to test
#' @param gtf.data the GTF data for the genome as GenomicRanges e.g. as read by rtracklayer
#'
#' @returns
#' @export
#'
#' @examples
get.transcript.junctions <- function(transcript.id, gtf.data) {
  # Find the exons for this transcript
  exons <- gtf.data[gtf.data$type == "exon" & gtf.data$transcript_id == transcript.id, ]

  if (transcript.is.reverse.strand(transcript.id, gtf.data)) {
    return(na.omit(data.frame(
      j1 = dplyr::lead(exons$end) + 1, # end of intron is start of next exon (end coordinate rev strand)
      j2 = exons$start # first base of intron is start of exon coordinate (rev strand)
    )))
  }

  data.frame(
    j1 = exons$end + 1, # start of the intron is the end of the exon
    j2 = dplyr::lead(exons$start) # end of the intron is the start of the next exon
  )
}

#' Is the given junction found in a genome GTF file?
#'
#' Given a splice junction, determine if the junction is already known within
#' the given transcript.
#'
#' @param transcript.id the transcript id to test
#' @param gtf.data the GTF data for the genome as GenomicRanges e.g. as read by rtracklayer
#' @param start the start coordinate of the junction
#' @param end the end coordinate of the junction
#'
#' @returns true if the junctions overlap any intron/exon boundaries in the transcript, false otherwise
#' @export
#'
#' @examples
#' junction.is.in.GTF("ENSGALT00010007119", chicken.gtf, 118299992, 118301364)
junction.is.in.GTF <- function(transcript.id, gtf.data, junction.start, junction.end) {
  junctions <- get.transcript.junctions(transcript.id, gtf.data)
  any(junctions$j1 == junction.start & junctions$j2 == junction.end)
}

#' Read ggshasimi output
#'
#' Read junction and coverage data from ggsashimi R implementation. Build exon
#' and intron annotations.
#'
#' @param coverage a ggsashimi coverage output
#' @param junctions a ggsashimi junctions output
#' @param gene.id the gene of interest
#' @param gtf.data the full genome annotation
#' @param location.string the location to restrict the search
#'
#' @returns
#' @export
#'
#' @examples
read.sashimi.data <- function(bam.file, gtf.data, chr, start, end,
                              reference.gene.id, reference.transcript.id = NA) {
  tryCatch(
    {
      cat("Reading sashimi data in region '", paste0(chr, ":", start, "-", end), "'\n")
      sashimi.data <- list()
      sashimi.data$input.file <- bam.file
      bam.data <- read_bam(bam.file, paste0(chr, ":", start, "-", end), "SENSE")

      sashimi.data$reference.gene.id <- reference.gene.id
      # Get the longest transcript in the gene if none specified
      if (is.na(reference.transcript.id)) {
        cat("Plot sashimi: No reference transcript given, selecting longest for gene id", reference.gene.id, "\n")
        reference.transcript.id <- gtf.data |>
          dplyr::filter(gene_id == reference.gene.id, type == "transcript") |>
          dplyr::mutate(length = end - start + 1) |>
          dplyr::arrange(length) |>
          dplyr::slice_tail(n = 1) |>
          dplyr::select(transcript_id) |>
          dplyr::pull()
      }
      sashimi.data$reference.transcript.id <- reference.transcript.id
      cat("Plot sashimi: Reference transcript is", reference.transcript.id, "\n")

      gene.gtf <- gtf.data[gtf.data$gene_id == reference.gene.id, ]
      if (nrow(gene.gtf) == 0) {
        stop("Plot sashimi: Unable to detect reference gene id in genome GTF\n")
      }

      # Is the gene on the forward or reverse strand? Note - this is the gene, not the junctions or reads
      sashimi.data$reference.transcript.strand <- gtf.data |>
        dplyr::filter(
          transcript_id == reference.transcript.id,
          type == "exon"
        ) |>
        dplyr::select(strand) |>
        dplyr::distinct() |>
        dplyr::pull(strand)

      cat("Plot sashimi: Reference transcript is on strand '", paste(sashimi.data$reference.transcript.strand, collapse = ","), "'\n")

      sashimi.data$reference.transcript.boundaries <- get_exon_boundaries(gtf.data, chr, start, end)

      cat("Plot sashimi: Reference exon/intron bounds detected:\n")

      sashimi.data$reference.gene.name <- gtf.data |>
        dplyr::filter(gene_id == reference.gene.id) |>
        dplyr::select(gene_name) |>
        dplyr::distinct() |>
        dplyr::pull()

      sashimi.data$coverage <- bam.data$coverage
      sashimi.data$junctions <- bam.data$junctions
      sashimi.data$junctions.reference.strand <- bam.data$reference.strand

      sashimi.data$junctions$matches.known.exon <- mapply(junction.is.in.GTF,
        junction.start = bam.data$junctions$start,
        junction.end = bam.data$junctions$end,
        MoreArgs = list(
          transcript.id = reference.transcript.id,
          gtf.data = gtf.data
        )
      )


      return(sashimi.data)
    },
    error = \(e) {
      cat("Plot sashimi: ", "Error reading sashimi data from", bam.file, "\n", paste(e))
      e
    }
  )
}

#' Read ggshasimi output
#'
#' Read a junction file and note if each junction matches transcript GTF coordinates
#'
#' @param rds.file a ggsashimi Rds output file
#'
#' @returns
#' @export
#'
#' @examples
# read.ggsashimi.py.rds <- function(rds.file) {
#' #  tryCatch(
#' #    {
#' #      rds.data <- readRDS(rds.file)
#' #      file.name.parts <- str_split_1(basename(rds.file), "\.")
#' #      rds.data$junction.strand <- ifelse(str_detect(rds.file, "_\+$"), "+",
#' #        ifelse(str_detect(rds.file, "_-$"), "-", "*")
#' #      )
#' #      cat("Plot sashimi: ", rds.file, "is for junctions on strand", rds.data$junction.strand, "\n")
#' #      rds.data$species <- file.name.parts[1]
#' #      rds.data$tissue <- file.name.parts[2]
#' #      rds.data$timepoint <- file.name.parts[3]
#' #      rds.data$gene.id <- file.name.parts[4]
#' #      rds.data$gene.name <- GENE.LOCATIONS |>
#' #        dplyr::filter(GeneId == rds.data$gene.id) |>
#' #        dplyr::select(Gene) |>
#' #        dplyr::pull()
#' #      rds.data$filename <- basename(rds.file)
#' #      rds.data$density_list <- rds.data$density_list[[1]]
#' #      rds.data$junction_list <- rds.data$junction_list[[1]]
#' #      rds.data$canonical.transcript.id <- GENE.LOCATIONS |>
#' #        dplyr::filter(GeneId == rds.data$gene.id) |>
#' #        dplyr::select(CanonicalTranscriptId) |>
#' #        dplyr::pull()
#' #      rds.data$gtf.data <- GTF.DATA[[rds.data$species]]
#' #      # Is the gene on the forward or reverse strand? Note - this is the gene, not the junctions or reads
#' #      rds.data$transcript.is.reverse.strand <- transcript.is.reverse.strand(rds.data$canonical.transcript.id, rds.data$gtf.data)
#' #      rds.data$transcript.strand <- ifelse(rds.data$transcript.is.reverse.strand, "-", "+")
#' #      rds.data$junction_list$matches.known.exons <- mapply(junction.is.in.GTF,
#' #        junction.start = rds.data$junction_list$x,
#' #        junction.end = rds.data$junction_list$xend,
#' #        MoreArgs = list(
#' #          transcript.id = rds.data$canonical.transcript.id,
#' #          gtf.data = rds.data$gtf.data
#' #        )
#' #      )
#' #      return(rds.data)
#' #    },
#' #    error = \(e) {
#' #      cat("Plot sashimi: ", "Error reading Rds data from", rds.file, "\n", paste(e))
#' #      e
#' #    }
#' #  )
# }

#' Make a conversion table that can rescale coordinates to collapse introns.
#'
#' @param exon.data data frame of exon coordinates with columns start and end
#' @param intron.data data frame of intron coordinates with columns start and end
#' @param strand the strand with the transcript to convert against
#' @param max.intron.length the maximum length for an intron
#'
#' @returns a function that will convert coordinates collapsing introns
#' @export
#'
#' @examples
create.intron.collapser <- function(exon.data, intron.data, strand, max.intron.length = 500) {
  if (is.null(exon.data)) stop("No exon data provided")

  cat("Plot sashimi: Creating intron collapser\n")
  # Reduce any overlapping exons if we have multiple transcripts
  exon.data <- exon.data[exon.data$strand == strand, ]
  intron.data <- intron.data[intron.data$strand == strand, ]

  cat("Plot sashimi: Detected", nrow(exon.data), "exons and", nrow(intron.data), "introns\n")

  exon.ranges <- GenomicRanges::reduce(GenomicRanges::GRanges(
    seqnames = rep("test", nrow(exon.data)),
    ranges = IRanges::IRanges(
      start = exon.data$start,
      end = exon.data$end
    ),
    strand = strand
  ))

  # Break introns apart, since they can be part of an exon for a different transcript
  intron.ranges <- GenomicRanges::disjoin(GenomicRanges::GRanges(
    seqnames = rep("test", nrow(intron.data)),
    ranges = IRanges::IRanges(
      start = intron.data$start,
      end = intron.data$end
    ),
    strand = strand
  ))

  # Remove introns that overlap an exon due to multiple transcripts
  intron.ranges$overlappingExons <- GenomicRanges::countOverlaps(intron.ranges, exon.ranges, minoverlap = 1)
  intron.ranges <- intron.ranges[intron.ranges$overlappingExons == 0, ]

  # Some introns may be missing. Fill in gaps from min start to max end that are
  # not covered by intron or exons
  missing.introns <- GenomicRanges::gaps(GenomicRanges::reduce(c(intron.ranges, exon.ranges)),
    start = min(exon.data$start)
  )

  # Convert back to data frames
  intron.ranges <- c(intron.ranges, missing.introns) |>
    as.data.frame() |>
    dplyr::select(start, end, strand) |>
    dplyr::mutate(Type = "intron")


  exon.ranges <- exon.ranges |>
    as.data.frame() |>
    dplyr::select(start, end, strand) |>
    dplyr::mutate(Type = "exon")

  full.ranges <- rbind(intron.ranges, exon.ranges)

  full.ranges <- full.ranges |>
    dplyr::arrange(start, end) |>
    dplyr::distinct() |>
    dplyr::mutate(
      original.length = end - start + 1,
      new.length = ifelse(original.length > 500 & Type == "intron", 500, original.length),
      new.end = min(start) + cumsum(new.length) - 1,
      new.start = new.end - new.length + 1, # how much offset to apply
      is.ordered = new.start < new.end,
      is.contiguous = new.start == dplyr::lag(new.end) + 1,
      is.full.coverage = start == dplyr::lag(end) + 1
    ) |>
    dplyr::select(everything(), new.start, new.end)

  # Create a function that uses the above tables to convert a coordinate to the new ranges
  calculate <- function(coordinate) {
    # If a coordinate is out of min bounds, don't adjust it
    if (coordinate < min(full.ranges$start)) {
      return(coordinate)
    }
    # If a coordinate is out of max bounds, reduce as needed
    if (coordinate > max(full.ranges$end)) {
      difference.from.end <- coordinate - max(full.ranges$end)
      return(max(full.ranges$new.end) + difference.from.end)
    }

    # Find the range within which to fit
    feature <- full.ranges %>%
      dplyr::filter(start <= coordinate & end >= coordinate) %>%
      dplyr::slice_head(n = 1)

    # If we find nothing, do not adjust
    if (nrow(feature) == 0) {
      warning("Plot sashimi: No range table entry covering", coordinate)
      return(coordinate)
    }

    # How far along the feature are we?
    fractional.distance <- (coordinate - feature$start) / feature$original.length

    # Nearest integer to the same fraction of the new coordinate space
    result <- (fractional.distance * feature$new.length) + feature$new.start
    return(result)
  }
  cat("Plot sashimi: Created intron collapser\n")
  list(
    calculate = calculate,
    full.ranges = full.ranges
  )
}


collapse.introns <- function(sashimi.data, exon.data, intron.data) {
  # Calculate offsets to make all introns at most 500bp
  intron.collapser <- create.intron.collapser(exon.data, intron.data, sashimi.data$reference.transcript.strand)

  cat("Plot sashimi: Collapsing introns\n")
  # Apply offsets to coordinates
  sashimi.data$reference.transcript.boundaries$exons <- sashimi.data$reference.transcript.boundaries$exons |>
    dplyr::rowwise() |>
    dplyr::mutate(
      old.start = start, old.end = end,
      start = intron.collapser$calculate(old.start),
      end = intron.collapser$calculate(old.end)
    )
  sashimi.data$reference.transcript.boundaries$introns <- sashimi.data$reference.transcript.boundaries$introns |>
    dplyr::rowwise() |>
    dplyr::mutate(
      old.start = start, old.end = end,
      start = intron.collapser$calculate(old.start),
      end = intron.collapser$calculate(old.end)
    )

  sashimi.data$junctions <- sashimi.data$junctions |>
    dplyr::rowwise() |>
    dplyr::mutate(
      old.start = start, old.end = end,
      start = intron.collapser$calculate(old.start),
      end = intron.collapser$calculate(old.end)
    )

  sashimi.data
}


#### Functions to build chart components ####

#' Plot splice junctions on a gene exon track
#'
#' @param sashimi.data the sashimi data from ggsashimi
#' @param min.spanning.reads  the minimum number of spanning reads to include
#' @param label the label for the gene track
#' @param show.x.axis if true, display the x axis
#' @param is.collapse.introns if true, make each intron at most 500bp wide
#'
#' @returns
#' @export
#'
#' @examples
make.gene.track.sashimi.panel <- function(sashimi.data, min.spanning.reads = 5, label = "tissue",
                                          show.x.axis = TRUE, is.collapse.introns = FALSE) {
  if (is.collapse.introns) {
    sashimi.data <- collapse.introns(
      sashimi.data,
      sashimi.data$reference.transcript.boundaries$exons,
      sashimi.data$reference.transcript.boundaries$introns
    )
  }

  junctions <- sashimi.data$junctions
  anns <- sashimi.data$reference.transcript.boundaries
  is.x.reverse <- sashimi.data$reference.transcript.strand == "-"

  if (any(!is.numeric(anns$exons$end)) | any(!is.numeric(anns$exons$start))) {
    cat("Plot sashimi: Error in annotations: at least one start or end is NA\n")
    print(anns$exons)
    str(anns$exons)
    # stop("Error in annotations: at least one start or end is NA")
  }

  # Only plot exons from the reference transcript
  reference.exons <- anns$exons |>
    dplyr::filter(
      transcript_id == sashimi.data$reference.transcript.id
    )

  non.reference.exons <- anns$exons |>
    dplyr::filter(
      gene_id != sashimi.data$reference.gene.id,
      strand != sashimi.data$reference.transcript.strand
    )

  reference.introns <- anns$introns |>
    dplyr::filter(
      start > min(reference.exons$start),
      end < max(reference.exons$start),
      strand == sashimi.data$reference.transcript.strand
    )

  non.reference.introns <- anns$introns |>
    dplyr::filter(
      start > min(reference.exons$start),
      end < max(reference.exons$start),
      strand != sashimi.data$reference.transcript.strand
    )

  # Set coordinate range for the x axis
  xmin <- min(c(reference.exons$start, reference.exons$end), na.rm = T) - 500
  xmax <- max(c(reference.exons$start, reference.exons$end), na.rm = T) + 500

  # Make the gene track
  splot <- ggplot() +
    # Introns
    geom_segment(data = reference.introns, aes(x = start, xend = end, y = 0.5, yend = 0.5), linewidth = 0.3) +
    geom_segment(data = non.reference.introns, aes(x = start, xend = end, y = -0.5, yend = -0.5), linewidth = 0.3) +
    geom_hline(yintercept = 0) +
    # annotate("text",
    #          x = ifelse(is.x.reverse, xmax, xmin),
    #          y = 0,
    #          label = label,
    #          size = 5, col = "black"
    # ) +
    annotate("text",
      x = ifelse(is.x.reverse, xmax - 250, xmin + 250),
      y = 0.5,
      label = sashimi.data$reference.transcript.strand,
      size = 2, col = "black"
    ) +
    annotate("text",
      x = ifelse(is.x.reverse, xmax - 250, xmin + 250),
      y = -0.5,
      label = ifelse(sashimi.data$reference.transcript.strand == "+", "-", "+"),
      size = 2, col = "black"
    ) +

    # Reference transcript exons
    geom_rect(data = reference.exons, aes(xmin = start, xmax = end, ymin = 0, ymax = 1), fill = "blue", alpha = 1) +
    geom_rect(data = non.reference.exons, aes(xmin = start, xmax = end, ymin = 0, ymax = -1), fill = "grey", alpha = 1) +
    scale_y_discrete(expand = c(-2, 2)) +
    coord_cartesian(xlim = c(xmin, xmax), ylim = c(-2, 2)) +
    scale_x_continuous(expand = c(0, 0.25)) +
    theme_minimal() +
    labs(y = label) +
    theme(
      axis.line.y = element_blank(),
      axis.line.x = element_line(),
      axis.ticks.x = element_line(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )

  if (is.x.reverse) {
    splot <- splot +
      scale_x_reverse(expand = c(0, 0.25))
  }

  # Hide x-axis if needed
  if (!show.x.axis) {
    splot <- splot +
      theme(
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
      )
  }

  # Add the sashimi splines. Each sashimi arc is made of two splines that can be
  # above or below the transcript
  add.junction <- function(splot, xmin, xmax, is.above, ymin, ymax, count, matches.known.exons) {
    # Define the spline line styles that make the junction lines
    spline.color <- ifelse(matches.known.exons, "blue", "black")
    spline.size <- ifelse(matches.known.exons, 1, 1)
    spline.alpha <- ifelse(matches.known.exons, 1, 1)

    # Define splice shapes.
    l.spline.btm <- grid::xsplineGrob(x = c(0, 0, 1, 1), y = c(1, 0, 0, 0), shape = 1, gp = gpar(lwd = spline.size, col = spline.color, alpha = spline.alpha))
    l.spline.top <- grid::xsplineGrob(x = c(0, 0, 1, 1), y = c(0, 1, 1, 1), shape = 1, gp = gpar(lwd = spline.size, col = spline.color, alpha = spline.alpha))
    r.spline.btm <- grid::xsplineGrob(x = c(1, 1, 0, 0), y = c(1, 0, 0, 0), shape = 1, gp = gpar(lwd = spline.size, col = spline.color, alpha = spline.alpha))
    r.spline.top <- grid::xsplineGrob(x = c(1, 1, 0, 0), y = c(0, 1, 1, 1), shape = 1, gp = gpar(lwd = spline.size, col = spline.color, alpha = spline.alpha))

    # What is the separation from the zero axis to the spline? This should match
    # the top of the exon box y position
    spline.y.offset <- 0
    label.y.offset <- 1.1 # multiplier for separation between spline and label

    # Determine which splines to use for minus strand versus plus strand transcripts
    l.grob.btm <- l.spline.btm
    if (is.x.reverse) l.grob.btm <- r.spline.btm

    l.grob.top <- l.spline.top
    if (is.x.reverse) l.grob.top <- r.spline.top

    r.grob.btm <- r.spline.btm
    if (is.x.reverse) r.grob.btm <- l.spline.btm

    r.grob.top <- r.spline.top
    if (is.x.reverse) r.grob.top <- l.spline.top

    xmid <- (xmin + xmax) / 2

    # Left arc
    if (is.above) { # Junctions below zero
      splot <- splot + annotation_custom(
        grob = l.grob.btm,
        xmin = xmin,
        xmax = xmid,
        ymin = -ymin,
        ymax = -ymax - spline.y.offset
      )
    } else {
      splot <- splot + annotation_custom(
        grob = l.grob.top,
        xmin = xmin,
        xmax = xmid,
        ymax = ymin,
        ymin = ymax + spline.y.offset
      )
    }

    # Right arc
    if (is.above) { # Junctions below zero
      splot <- splot + annotation_custom(
        grob = r.grob.btm,
        xmin = xmid,
        xmax = xmax,
        ymin = -ymin,
        ymax = -ymax - spline.y.offset
      )
    } else {
      splot <- splot + annotation_custom(
        grob = r.grob.top,
        xmin = xmid,
        xmax = xmax,
        ymax = ymin,
        ymin = ymax + spline.y.offset
      )
    }

    splot <- splot + annotate("label",
      x = xmid,
      y = ifelse(is.above, -ymax * label.y.offset, ymax * label.y.offset),
      label = as.character(count),
      size = 1, col = spline.color, fill = NA, label.size = NA
    )

    return(splot)
  }

  # Check junctions meet plot criteria
  if (nrow(junctions) > 0) {
    # junctions$x <- as.vector(as.numeric(junctions$x))

    junctions <- junctions %>%
      dplyr::filter(count >= min.spanning.reads) %>%
      na.omit() %>%
      dplyr::arrange(start, end)


    if (any(!is.numeric(junctions$end)) | any(!is.numeric(junctions$start))) {
      print(junctions)
      str(junctions)
      stop("Plot sashimi: Error in junctions: at least one x or xend is NA")
    }

    junctions <- junctions %>%
      dplyr::mutate(length = end - start + 1)
  }

  # Add the junctions, adjusting for plus vs minus strand
  if (nrow(junctions) > 0) {
    # Calculate charting coordinates
    # junctions$isEven <- sapply(1:nrow(junctions), \(x) x %% 2 == 0)
    junctions$y.offset <- rep(seq(0, 2, 0.5), length.out = nrow(junctions)) # give each junction a separate y offset


    for (i in 1:nrow(junctions)) {
      jrow <- junctions[i, ]

      splot <- add.junction(
        splot = splot,
        xmin = jrow$start, xmax = jrow$end,
        is.above = jrow$strand == sashimi.data$reference.transcript.strand,
        ymin = 1, ymax = 1.5 + jrow$y.offset,
        count = jrow$count,
        matches.known.exons = jrow$matches.known.exon
      )
    }
  }

  list(plot = splot, junctions = junctions)
}


#### Functions to create multi-panel charts ####
#
# # Create a sashimi panel plot for all tissues of the given species and timepoint
# make.species.panel <- function(species, timepoint, gene.id) {
#   data.files <- list.files(path = "data/merged", pattern = paste0(species, ".*\\.", timepoint, "\\.", gene.id, "\\..*Rds_*"), full.names = TRUE)
#   if (length(data.files) == 0) {
#     return()
#   }
#   data <- lapply(data.files, read.ggsashimi.py.rds)
#
#   gene.name <- GENE.LOCATIONS[GENE.LOCATIONS$GeneId == gene.id, ]$Gene
#
#   out.png.file <- paste0("report/species/", species, ".", timepoint, ".", gene.id, ".", gene.name, ".png")
#
#   plots <- lapply(data, \(x)  make.sashimi.panel(x, label = paste0(x$tissue, " ", x$junction.strand))$plot)
#   track <- make.gene.track(data[[1]]) # only one gene, only need one track
#   plots[[length(plots) + 1]] <- track
#
#   patchwork::wrap_plots(plots, nrow = length(plots))
#   save.double.width(filename = out.png.file, plot = last_plot(), height = 300)
# }
#
# # Create a sashimi panel plot for all species of the given tissue and timepoint
# make.tissue.panel <- function(tissue, timepoint) {
#   cat("Plot sashimi: Making", tissue, "at", timepoint, "\n")
#   data.files <- list.files(path = "data/merged", pattern = paste0(".*\\.", tissue, "\\.", timepoint, "\\..*Rds_*"), full.names = TRUE)
#   if (length(data.files) == 0) {
#     return()
#   }
#   data <- lapply(data.files, read.ggsashimi.py.rds)
#
#   out.png.file <- paste0("report/tissues/", tissue, ".", timepoint, ".png")
#
#   plots <- lapply(data, \(x) make.sashimi.panel(x, label = paste0(x$species, "\n", x$gene.name, " ", x$junction.strand))$plot)
#   tracks <- lapply(data, make.gene.track)
#   plots <- c(rbind(plots, tracks))
#
#   patchwork::wrap_plots(plots, nrow = length(plots))
#   save.double.width(filename = out.png.file, plot = last_plot(), height = 300)
# }
#
# # Create a sashimi panel plot for all timepoint of the given tissue and species
# make.timepoint.panel <- function(species, tissue, gene.id) {
#   data.files <- list.files(path = "data/merged", pattern = paste0(species, ".*\\.", tissue, "\\..*", gene.id, "\\..*Rds_*"), full.names = TRUE)
#   if (length(data.files) == 0) {
#     return()
#   }
#
#   # Ensure files are plotted in time order
#   ordered.files <- list()
#   for (t in TIME.ORDER) {
#     for (f in data.files) {
#       if (str_detect(f, t)) {
#         ordered.files <- c(ordered.files, f)
#       }
#     }
#   }
#
#   data <- lapply(ordered.files, read.ggsashimi.py.rds)
#
#   gene.name <- GENE.LOCATIONS[GENE.LOCATIONS$GeneId == gene.id, ]$Gene
#
#   out.png.file <- paste0("report/timepoints/", species, ".", tissue, ".", gene.id, ".", gene.name, ".png")
#
#   plots <- lapply(data, \(x) make.sashimi.panel(x, label = paste0(x$species, " ", gene.name, "\n", x$timepoint, " ", x$junction.strand))$plot)
#   track <- make.gene.track(data[[1]])
#   plots[[length(plots) + 1]] <- track
#
#   patchwork::wrap_plots(plots, nrow = length(plots))
#   save.double.width(filename = out.png.file, plot = last_plot(), height = 300)
# }

#### Select groups for plotting ####

# find all species combinations for plotting

# all.samples <- merge(make.sample.groups(), GENE.LOCATIONS, by = "CommonName")
# cat("Plot sashimi: Selected samples for plotting\n")
#
# species.groups <- all.samples %>%
#   dplyr::group_by(CommonName, Timepoint, GeneId) %>%
#   dplyr::summarise(SampleCount = n(), .groups = "drop_last")
#
# tissue.groups <- all.samples %>%
#   dplyr::group_by(Organism_part, Timepoint) %>%
#   dplyr::summarise(SampleCount = n(), .groups = "drop_last")
#
# timepoint.groups <- all.samples %>%
#   dplyr::group_by(CommonName, Organism_part, GeneId) %>%
#   dplyr::summarise(SampleCount = n(), .groups = "drop_last")

#### Functions to condense introns for neater plotting ####

# Make species plot using combined panels
# make.condensed.species.panels <- function(species, timepoint, gene.id) {
#   data.files <- list.files(path = "data/merged", pattern = paste0(species, ".*\\.", timepoint, "\\.", gene.id, "\\..*Rds_*"), full.names = TRUE)
#   if (length(data.files) == 0) {
#     return()
#   }
#   data <- lapply(data.files, read.ggsashimi.py.rds)
#
#   gene.name <- GENE.LOCATIONS[GENE.LOCATIONS$GeneId == gene.id, ]$Gene
#
#   out.png.file <- paste0("report/species/", species, ".", timepoint, ".", gene.id, ".", gene.name, ".condensed.png")
#   cat("Plot sashimi: Making", out.png.file, "\n")
#
#   plots <- lapply(data, \(x)  make.gene.track.sashimi.panel(x, label = paste0(species, " ", gene.name, "\n", x$tissue, " ", x$junction.strand), show.x.axis = FALSE, is.collapse.introns = TRUE)$plot)
#
#   patchwork::wrap_plots(plots, nrow = length(plots))
#
#   save.double.width(filename = out.png.file, plot = last_plot(), height = 300)
# }
#
# # Make tissue plot using condensed panels
# make.condensed.tissue.panels <- function(tissue, timepoint) {
#   cat("Plot sashimi: Making", tissue, "at", timepoint, "\n")
#   data.files <- list.files(path = "data/merged", pattern = paste0(".*\\.", tissue, "\\.", timepoint, "\\..*Rds_*"), full.names = TRUE)
#   if (length(data.files) == 0) {
#     return()
#   }
#   data <- lapply(data.files, read.ggsashimi.py.rds)
#
#   out.png.file <- paste0("report/tissues/", tissue, ".", timepoint, ".condensed.png")
#
#   plots <- lapply(data, \(x) make.gene.track.sashimi.panel(x,
#     label = paste0(x$species, "\n", x$gene.name, " ", x$junction.strand),
#     show.x.axis = FALSE, is.collapse.introns = TRUE
#   )$plot)
#
#   patchwork::wrap_plots(plots, nrow = length(plots))
#   save.double.width(filename = out.png.file, plot = last_plot(), height = 300)
# }
#
# # Make species plot using combined panels
# make.condensed.timepoint.panels <- function(species, tissue, gene.id) {
#   cat("Plot sashimi: Making", species, tissue, "for", gene.id, "\n")
#   data.files <- list.files(path = "data/merged", pattern = paste0(species, ".*\\.", tissue, "\\..*", gene.id, "\\..*Rds_*"), full.names = TRUE)
#   if (length(data.files) == 0) {
#     return()
#   }
#
#   # Ensure files are plotted in time order
#   ordered.files <- list()
#   for (t in TIME.ORDER) {
#     for (f in data.files) {
#       if (str_detect(f, t)) {
#         ordered.files <- c(ordered.files, f)
#       }
#     }
#   }
#
#   data <- lapply(ordered.files, read.ggsashimi.py.rds)
#
#   gene.name <- GENE.LOCATIONS[GENE.LOCATIONS$GeneId == gene.id, ]$Gene
#
#   out.png.file <- paste0("report/timepoints/", species, ".", tissue, ".", gene.id, ".", gene.name, ".condensed.png")
#
#   plots <- lapply(data, \(x) make.gene.track.sashimi.panel(x,
#     label = paste0(x$species, " ", gene.name, "\n", x$timepoint, " ", x$junction.strand),
#     show.x.axis = FALSE, is.collapse.introns = TRUE
#   )$plot)
#
#   patchwork::wrap_plots(plots, nrow = length(plots))
#   save.double.width(filename = out.png.file, plot = last_plot(), height = 300)
# }

#### Make condensed figures ####
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
  coords <- parse.coordinates(gene.data$Location)

  cat("Reference location:", coords$coord.chr, coords$coord.start, coords$coord.end, "\n")

  sashimi.data <- read.sashimi.data(
    bam.file = bam.row$path,
    gtf.data = GTF.DATA[[species]],
    chr = coords$coord.chr, start = coords$coord.start, end = coords$coord.end,
    reference.gene.id = gene.data$GeneId,
    reference.transcript.id = gene.data$CanonicalTranscriptId
  )

  # Create plot with collapsed introns
  sashimi.plot.collapsed <- make.gene.track.sashimi.panel(sashimi.data,
    is.collapse.introns = TRUE, show.x.axis = FALSE,
    min.spanning.reads = 5, label = paste0(species, "\n", tissue, "\n", timepoint)
  )

  save.double.width(paste0("report/raw_sashimi/", species, ".", tissue, ".", timepoint, ".", gene_id, ".condensed.png"),
    sashimi.plot.collapsed$plot,
    height = 80
  )
  # Create plot with expanded introns
  sashimi.plot.expanded <- make.gene.track.sashimi.panel(sashimi.data,
    is.collapse.introns = FALSE,
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
