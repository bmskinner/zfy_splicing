# Make sashimi plot based on ggsashimi output.
# ggsashimi.py modified to write data objects to Rds when run
# This custom sashimi plot ensures the transcripts are always left to right
# irrespective of strand
cat("Plot sashimi: Beginning\n")
source("src/functions.R")

cat("Plot sashimi: running shashimi plotting\n")

ANNOTATED.EXONS <- get.annotated.exons()

# Ensure output dirs exist
fs::dir_create(c("report/species", "report/timepoints", "report/tissues"))

# Read all the GTF files to a global variable
read.gtf.data <- function() {
  cat("Plot sashimi: Reading full genome GTF files\n")
  gtf.data <- mclapply(GENOME.DATA$GTF_FILE, rtracklayer::import,
    mc.cores = ifelse(installr::is.windows(), 1, 6)
  )
  names(gtf.data) <- dplyr::pull(GENOME.DATA[, "CommonName"])
  cat("Plot sashimi: Read full genome GTF files\n")
  gtf.data
}
GTF.DATA <- read.gtf.data()
# GTF.DATA <- list()
# GTF.DATA[["chicken"]]<-chicken.gtf
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
#' Read a junction file and note if each junction matches transcript GTF coordinates
#'
#' @param rds.file a ggsashimi Rds output file
#'
#' @returns
#' @export
#'
#' @examples
read.rds.file <- function(rds.file) {
  tryCatch(
    {
      rds.data <- readRDS(rds.file)

      file.name.parts <- str_split_1(basename(rds.file), "\\.")

      rds.data$junction.strand <- ifelse(str_detect(rds.file, "_\\+$"), "+",
        ifelse(str_detect(rds.file, "_-$"), "-", "*")
      )

      cat("Plot sashimi: ", rds.file, "is for junctions on strand", rds.data$junction.strand, "\n")
      rds.data$species <- file.name.parts[1]
      rds.data$tissue <- file.name.parts[2]
      rds.data$timepoint <- file.name.parts[3]
      rds.data$gene.id <- file.name.parts[4]
      rds.data$gene.name <- GENE.LOCATIONS |>
        dplyr::filter(GeneId == rds.data$gene.id) |>
        dplyr::select(Gene) |>
        dplyr::pull()

      rds.data$filename <- basename(rds.file)

      rds.data$density_list <- rds.data$density_list[[1]]
      rds.data$junction_list <- rds.data$junction_list[[1]]

      rds.data$canonical.transcript.id <- GENE.LOCATIONS |>
        dplyr::filter(GeneId == rds.data$gene.id) |>
        dplyr::select(CanonicalTranscriptId) |>
        dplyr::pull()

      rds.data$gtf.data <- GTF.DATA[[rds.data$species]]

      # Is the gene on the forward or reverse strand? Note - this is the gene, not the junctions or reads
      rds.data$transcript.is.reverse.strand <- transcript.is.reverse.strand(rds.data$canonical.transcript.id, rds.data$gtf.data)
      rds.data$transcript.strand <- ifelse(rds.data$transcript.is.reverse.strand, "-", "+")

      rds.data$junction_list$matches.known.exons <- mapply(junction.is.in.GTF,
        junction.start = rds.data$junction_list$x,
        junction.end = rds.data$junction_list$xend,
        MoreArgs = list(
          transcript.id = rds.data$canonical.transcript.id,
          gtf.data = rds.data$gtf.data
        )
      )


      return(rds.data)
    },
    error = \(e) {
      cat("Plot sashimi: ", "Error reading Rds data from", rds.file, "\n", paste(e))
      e
    }
  )
}

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
  # Reduce any overlapping exons if we have multiple transcripts
  exon.data <- exon.data[exon.data$strand == strand, ]
  intron.data <- intron.data[intron.data$strand == strand]

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
  intron.ranges$overlappingExons <- GenomicRanges::countOverlaps(intron.ranges, exon.ranges, minoverlap = 10)
  intron.ranges <- intron.ranges[intron.ranges$overlappingExons == 0, ]

  # Some introns may be missing. Fill in gaps from min start to max end that are
  # not covered by intron or exons
  missing.introns <- GenomicRanges::gaps(GenomicRanges::reduce(c(intron.ranges, exon.ranges)),
    start = min(exons$start)
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
      new.end = min(start) + cumsum(new.length),
      new.start = new.end - new.length # how much offset to apply
    )

  # Create a function that uses the above tables to convert a coordinate to the new ranges
  function(coordinate) {
    # If a coordinate is out of bounds, don't adjust it
    if (coordinate < min(full.ranges$start)) {
      return(coordinate)
    }

    if (coordinate > max(full.ranges$end)) {
      old.dist <- max(full.ranges$end) - coordinate
      return(max(full.ranges$new.end) + old.dist)
    }

    # Take the first matching row
    feature <- full.ranges %>%
      dplyr::filter(start <= coordinate & end >= coordinate) %>%
      dplyr::slice_head(n = 1)

    # If we find nothing, do not adjust
    if (nrow(feature) == 0) {
      return(coordinate)
    }

    # How far along the feature are we?
    f.feature <- (coordinate - feature$start) / feature$original.length

    # Nearest integer to the same fraction of the new coordinate space
    result <- unique(round(f.feature * feature$new.length + feature$new.start))

    # If this was NA becuase we are outside the bounds of the main transcript,
    # return the original coordinate
    # if (is.na(result)) {
    #   return(coordinate)
    # }
    return(result)
  }
}


collapse.introns <- function(sashimi.data, exon.data, intron.data) {
  # Calculate offsets to make all introns at most 500bp
  intron.collapser <- create.intron.collapser(exon.data, intron.data, sashimi.data$transcript.strand)

  # Copy the existing coordinates
  junctions <- sashimi.data$junction_list
  anns <- sashimi.data$ann_list

  # Apply offsets to coordinates
  anns$exons <- as.data.frame(anns$exons)
  anns$introns <- as.data.frame(anns$introns)

  anns$exons$old.start <- anns$exons$start
  anns$exons$old.end <- anns$exons$end
  anns$exons$start <- mapply(intron.collapser, coordinate = anns$exons$old.start, SIMPLIFY = TRUE)
  anns$exons$end <- mapply(intron.collapser, coordinate = anns$exons$old.end, SIMPLIFY = TRUE)

  anns$introns$old.start <- anns$introns$start
  anns$introns$old.end <- anns$introns$end
  anns$introns$start <- mapply(intron.collapser, coordinate = anns$introns$old.start, SIMPLIFY = TRUE)
  anns$introns$end <- mapply(intron.collapser, coordinate = anns$introns$old.end, SIMPLIFY = TRUE)

  junctions$old.x <- junctions$x
  junctions$old.xend <- junctions$xend
  junctions$x <- sapply(junctions$old.x, intron.collapser)
  junctions$xend <- sapply(junctions$old.xend, intron.collapser)

  list(
    anns = anns,
    junctions = junctions
  )
}


#### Functions to build chart components ####

#
#' Create a ggplot with transcripts covering the region in ggsashimi data.
#' The reference transcript is highlighted.
#'
#' @param sashimi.data data from ggsashimi Rds output
#'
#' @returns a ggplot showing all transcripts
#' @export
#'
#' @examples
make.gene.track <- function(sashimi.data) {
  data <- sashimi.data$density_list
  junctions <- sashimi.data$junction_list
  anns <- sashimi.data$ann_list

  # Set coordinates for the x axis, reversing if on reverse strand
  xmin <- ifelse(sashimi.data$transcript.is.reverse.strand, max(data$x), min(data$x))
  xmax <- ifelse(sashimi.data$transcript.is.reverse.strand, min(data$x), max(data$x))

  # Create the arrows to draw on introns
  make.tx.arrows <- function() {
    # Create data table for strand arrows
    txarrows <- data.table()
    introns <- sashimi.data$ann_list[["introns"]]
    # Add right-pointing arrows for plus strand
    if ("+" %in% introns$strand && nrow(introns[strand == "+" & end - start > 5, ]) > 0) {
      txarrows <- rbind(
        txarrows,
        introns[strand == "+" & end - start > 5, list(
          seq(start + 4, end, by = 453) - 1,
          seq(start + 4, end, by = 453)
        ), by = .(tx, start, end)]
      )
    }
    # Add left-pointing arrows for minus strand
    if ("-" %in% introns$strand && nrow(introns[strand == "-" & end - start > 5, ]) > 0) {
      txarrows <- rbind(
        txarrows,
        introns[strand == "-" & end - start > 5, list(
          seq(start, max(start + 1, end - 4), by = 453),
          seq(start, max(start + 1, end - 4), by = 453) - 1
        ), by = .(tx, start, end)]
      )
    }
    txarrows
  }

  # Make the gene track
  sashimi.plot <- ggplot() +
    # Introns and arrows
    geom_segment(data = sashimi.data$ann_list$introns, aes(x = start, xend = end, y = tx, yend = tx), linewidth = 0.3) +
    # geom_segment(data=make.tx.arrows(), aes(x=V1,xend=V2,y=tx,yend=tx), arrow=arrow(length=unit(0.02,"npc")))+

    # Exons
    geom_segment(
      data = sashimi.data$ann_list$exons, aes(
        x = start, xend = end,
        y = tx, yend = tx,
        col = tx == sashimi.data$canonical.transcript.id
      ),
      linewidth = 5, alpha = 1
    ) +
    scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "grey")) +
    scale_y_discrete(expand = c(0, 0.5)) +
    coord_cartesian(xlim = c(xmin, xmax)) +
    scale_x_continuous(expand = c(0, 0.25), labels = comma) +
    theme_minimal() +
    theme(
      axis.line.y = element_blank(),
      axis.line.x = element_line(),
      axis.ticks.x = element_line(),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )

  if (sashimi.data$transcript.is.reverse.strand) {
    sashimi.plot <- sashimi.plot +
      scale_x_reverse()
  }
  sashimi.plot
}

# Create a sashimi panel for the given ggsashimi data
# min.spanning.reads - the minimum number of junction-spanning reads needed to display on the chart
# y.label - the data value to be used for the y axis label
make.sashimi.panel <- function(sashimi.data, min.spanning.reads = 5, label = "tissue") {
  # cat("Making sashimi panel\n")
  data <- sashimi.data$density_list
  junctions <- sashimi.data$junction_list
  anns <- sashimi.data$ann_list

  # Set coordinates for the x axis
  is.minus.strand <- sashimi.data$transcript.is.reverse.strand
  xmin <- min(data$x)
  xmax <- max(data$x)
  ymax <- max(data$y)

  # Create the coverage plot
  splot <- ggplot() +
    geom_bar(data = data, aes(x, y), width = 1, position = "identity", stat = "identity", alpha = 0.5) +
    coord_cartesian(
      ylim = c(0 - (ymax), ymax * 1.5),
      xlim = c(xmin, xmax)
    ) +
    scale_x_continuous(expand = c(0, 0.25)) +
    labs(y = label)

  if (sashimi.data$transcript.is.reverse.strand) {
    splot <- splot +
      scale_x_reverse()
  }


  add.junction <- function(splot, xmin, xmax, ymin, ymax, is.even, count, matches.known.exons) {
    # Define the spline shapes that make the junction lines
    spline.color <- ifelse(matches.known.exons, "blue", "black")
    spline.size <- ifelse(matches.known.exons, 1, 2)

    l.spline.btm <- grid::xsplineGrob(x = c(0, 0, 1, 1), y = c(1, 0, 0, 0), shape = 1, gp = gpar(lwd = spline.size, col = spline.color))
    l.spline.top <- grid::xsplineGrob(x = c(0, 0, 1, 1), y = c(0, 1, 1, 1), shape = 1, gp = gpar(lwd = spline.size, col = spline.color))
    r.spline.btm <- grid::xsplineGrob(x = c(1, 1, 0, 0), y = c(1, 0, 0, 0), shape = 1, gp = gpar(lwd = spline.size, col = spline.color))
    r.spline.top <- grid::xsplineGrob(x = c(1, 1, 0, 0), y = c(0, 1, 1, 1), shape = 1, gp = gpar(lwd = spline.size, col = spline.color))

    # Determine which splines to use for minus strand versus plus strand transcripts
    l.grob.btm <- l.spline.btm
    if (is.minus.strand) l.grob.btm <- r.spline.btm

    l.grob.top <- l.spline.top
    if (is.minus.strand) l.grob.top <- r.spline.top

    r.grob.btm <- r.spline.btm
    if (is.minus.strand) r.grob.btm <- l.spline.btm

    r.grob.top <- r.spline.top
    if (is.minus.strand) r.grob.top <- l.spline.top

    xmid <- (xmin + xmax) / 2
    ymid <- ((ymin + ymax) / 2) * 1.25

    # Left arc
    if (is.even) { # Junctions below zero
      splot <- splot + annotation_custom(
        grob = l.grob.btm,
        xmin = xmin,
        xmax = xmid,
        ymin = -ymid * 0.5,
        ymax = 0
      )
    } else {
      splot <- splot + annotation_custom(
        grob = l.grob.top,
        xmin = xmin,
        xmax = xmid,
        ymax = ymid,
        ymin = 0
      )
    }

    # Right arc
    if (is.even) { # Junctions below zero
      splot <- splot + annotation_custom(
        grob = r.grob.btm,
        xmin = xmid,
        xmax = xmax,
        ymin = -ymid * 0.5,
        ymax = 0
      )
    } else {
      splot <- splot + annotation_custom(
        grob = r.grob.top,
        xmin = xmid,
        xmax = xmax,
        ymax = ymid,
        ymin = 0
      )
    }

    splot <- splot + annotate("label",
      x = xmid,
      y = ifelse(is.even, -ymid * 1.1, ymid * 1.2),
      label = as.character(count),
      size = 2, col = spline.color, fill = NA, label.size = NA
    )

    return(splot)
  }

  if (nrow(junctions) > 0) {
    junctions <- junctions %>%
      dplyr::filter(count >= min.spanning.reads)
  }

  # Add the junctions, adjusting for plus vs minus strand
  if (nrow(junctions) > 0) {
    # Calculate charting coordinates
    junctions$total.count <- sum(junctions$count)
    junctions$isEven <- sapply(1:nrow(junctions), \(x) x %% 2 == 0)

    for (i in 1:nrow(junctions)) {
      splot <- add.junction(
        splot, junctions[i, ]$x, junctions[i, ]$xend, junctions[i, ]$y,
        junctions[i, ]$yend, junctions[i, ]$isEven, junctions[i, ]$count,
        junctions[i, ]$matches.known.exons
      )
    }
  }

  # Format the final plot
  splot <- splot +
    theme_minimal() +
    theme(
      axis.line = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  list(plot = splot, junctions = junctions)
}


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
  junctions <- sashimi.data$junction_list
  anns <- sashimi.data$ann_list

  if (is.collapse.introns) {
    collapse.data <- collapse.introns(sashimi.data, sashimi.data$ann_list$exons, sashimi.data$ann_list$introns)
    junctions <- collapse.data$junctions
    anns <- collapse.data$anns
  }



  if (any(!is.numeric(anns$exons$end)) | any(!is.numeric(anns$exons$start))) {
    cat("Plot sashimi: Error in annotations: at least one start or end is NA\n")
    print(anns$exons)
    str(anns$exons)
    # stop("Error in annotations: at least one start or end is NA")
  }

  # Set coordinate range for the x axis
  xmin <- min(anns$exons$start, na.rm = T) - 500
  xmax <- max(anns$exons$end, na.rm = T) + 500

  # Only plot exons from the reference transcript
  canonical.exons <- anns$exons |>
    dplyr::filter(
      tx == sashimi.data$canonical.transcript.id,
      strand == sashimi.data$transcript.strand
    )

  canonical.introns <- anns$introns |> dplyr::filter(
    start > min(canonical.exons$start),
    end < max(canonical.exons$start),
    strand == sashimi.data$transcript.strand
  )

  # Make the gene track
  splot <- ggplot() +
    # Introns
    geom_segment(data = canonical.introns, aes(x = start, xend = end, y = 0, yend = 0), linewidth = 0.3) +

    # Reference transcript exons
    geom_rect(data = canonical.exons, aes(xmin = start, xmax = end, ymin = -0.5, ymax = 0.5), fill = "blue", alpha = 1) +
    scale_y_discrete(expand = c(0.5, 0.5)) +
    coord_cartesian(xlim = c(xmin, xmax)) +
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

  if (sashimi.data$transcript.is.reverse.strand) {
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
  add.junction <- function(splot, xmin, xmax, ymin, ymax, is.even, count, matches.known.exons, f.length) {
    # Define the spline shapes that make the junction lines
    spline.color <- ifelse(matches.known.exons, "blue", "black")
    spline.size <- ifelse(matches.known.exons, 1, 2)
    spline.alpha <- ifelse(matches.known.exons, 0.5, 1)

    l.spline.btm <- grid::xsplineGrob(x = c(0, 0, 1, 1), y = c(1, 0, 0, 0), shape = 1, gp = gpar(lwd = spline.size, col = spline.color, alpha = spline.alpha))
    l.spline.top <- grid::xsplineGrob(x = c(0, 0, 1, 1), y = c(0, 1, 1, 1), shape = 1, gp = gpar(lwd = spline.size, col = spline.color, alpha = spline.alpha))
    r.spline.btm <- grid::xsplineGrob(x = c(1, 1, 0, 0), y = c(1, 0, 0, 0), shape = 1, gp = gpar(lwd = spline.size, col = spline.color, alpha = spline.alpha))
    r.spline.top <- grid::xsplineGrob(x = c(1, 1, 0, 0), y = c(0, 1, 1, 1), shape = 1, gp = gpar(lwd = spline.size, col = spline.color, alpha = spline.alpha))

    spline.y.offset <- 0.5 # separation between exon and spline
    label.y.offset <- 1.5 # separation between spline and label

    # Determine which splines to use for minus strand versus plus strand transcripts
    l.grob.btm <- l.spline.btm
    if (sashimi.data$transcript.is.reverse.strand) l.grob.btm <- r.spline.btm

    l.grob.top <- l.spline.top
    if (sashimi.data$transcript.is.reverse.strand) l.grob.top <- r.spline.top

    r.grob.btm <- r.spline.btm
    if (sashimi.data$transcript.is.reverse.strand) r.grob.btm <- l.spline.btm

    r.grob.top <- r.spline.top
    if (sashimi.data$transcript.is.reverse.strand) r.grob.top <- l.spline.top

    xmid <- (xmin + xmax) / 2

    # Left arc
    if (is.even) { # Junctions below zero
      splot <- splot + annotation_custom(
        grob = l.grob.btm,
        xmin = xmin,
        xmax = xmid,
        ymin = -(1 + f.length),
        ymax = -spline.y.offset
      )
    } else {
      splot <- splot + annotation_custom(
        grob = l.grob.top,
        xmin = xmin,
        xmax = xmid,
        ymax = 1 + f.length,
        ymin = spline.y.offset
      )
    }

    # Right arc
    if (is.even) { # Junctions below zero
      splot <- splot + annotation_custom(
        grob = r.grob.btm,
        xmin = xmid,
        xmax = xmax,
        ymin = -(1 + f.length),
        ymax = -spline.y.offset
      )
    } else {
      splot <- splot + annotation_custom(
        grob = r.grob.top,
        xmin = xmid,
        xmax = xmax,
        ymax = 1 + f.length,
        ymin = spline.y.offset
      )
    }

    splot <- splot + annotate("label",
      x = xmid,
      y = ifelse(is.even, -label.y.offset - f.length, label.y.offset + f.length),
      label = as.character(count),
      size = 2, col = spline.color, fill = NA, label.size = NA
    )

    return(splot)
  }

  # Check junctions meet plot criteria
  if (nrow(junctions) > 0) {
    junctions$x <- as.vector(as.numeric(junctions$x))

    junctions <- junctions %>%
      dplyr::filter(count >= min.spanning.reads) %>%
      na.omit() %>%
      dplyr::arrange(x, xend)


    if (any(!is.numeric(junctions$xend)) | any(!is.numeric(junctions$x))) {
      print(junctions)
      str(junctions)
      stop("Plot sashimi: Error in junctions: at least one x or xend is NA")
    }

    junctions <- junctions %>%
      dplyr::mutate(length = abs(xend - x))
  }

  # Add the junctions, adjusting for plus vs minus strand
  if (nrow(junctions) > 0) {
    # Calculate charting coordinates
    junctions$isEven <- sapply(1:nrow(junctions), \(x) x %% 2 == 0)
    junctions$y.offset <- rep(seq(0, 1, 0.25), length.out = nrow(junctions)) # give each junction a separate y offset

    for (i in 1:nrow(junctions)) {
      splot <- add.junction(
        splot, junctions[i, ]$x, junctions[i, ]$xend,
        junctions[i, ]$y, junctions[i, ]$yend,
        junctions[i, ]$isEven, junctions[i, ]$count,
        junctions[i, ]$matches.known.exons, junctions[i, ]$y.offset
      )
    }
  }

  list(plot = splot, junctions = junctions)
}



#### Functions to create multi-panel charts ####

# Create a sashimi panel plot for all tissues of the given species and timepoint
make.species.panel <- function(species, timepoint, gene.id) {
  data.files <- list.files(path = "data/merged", pattern = paste0(species, ".*\\.", timepoint, "\\.", gene.id, "\\..*Rds_*"), full.names = TRUE)
  if (length(data.files) == 0) {
    return()
  }
  data <- lapply(data.files, read.rds.file)

  gene.name <- GENE.LOCATIONS[GENE.LOCATIONS$GeneId == gene.id, ]$Gene

  out.png.file <- paste0("report/species/", species, ".", timepoint, ".", gene.id, ".", gene.name, ".png")

  plots <- lapply(data, \(x)  make.sashimi.panel(x, label = paste0(x$tissue, " ", x$junction.strand))$plot)
  track <- make.gene.track(data[[1]]) # only one gene, only need one track
  plots[[length(plots) + 1]] <- track

  patchwork::wrap_plots(plots, nrow = length(plots))
  save.double.width(filename = out.png.file, plot = last_plot(), height = 300)
}

# Create a sashimi panel plot for all species of the given tissue and timepoint
make.tissue.panel <- function(tissue, timepoint) {
  cat("Plot sashimi: Making", tissue, "at", timepoint, "\n")
  data.files <- list.files(path = "data/merged", pattern = paste0(".*\\.", tissue, "\\.", timepoint, "\\..*Rds_*"), full.names = TRUE)
  if (length(data.files) == 0) {
    return()
  }
  data <- lapply(data.files, read.rds.file)

  out.png.file <- paste0("report/tissues/", tissue, ".", timepoint, ".png")

  plots <- lapply(data, \(x) make.sashimi.panel(x, label = paste0(x$species, "\n", x$gene.name, " ", x$junction.strand))$plot)
  tracks <- lapply(data, make.gene.track)
  plots <- c(rbind(plots, tracks))

  patchwork::wrap_plots(plots, nrow = length(plots))
  save.double.width(filename = out.png.file, plot = last_plot(), height = 300)
}

# Create a sashimi panel plot for all timepoint of the given tissue and species
make.timepoint.panel <- function(species, tissue, gene.id) {
  data.files <- list.files(path = "data/merged", pattern = paste0(species, ".*\\.", tissue, "\\..*", gene.id, "\\..*Rds_*"), full.names = TRUE)
  if (length(data.files) == 0) {
    return()
  }

  # Ensure files are plotted in time order
  ordered.files <- list()
  for (t in TIME.ORDER) {
    for (f in data.files) {
      if (str_detect(f, t)) {
        ordered.files <- c(ordered.files, f)
      }
    }
  }

  data <- lapply(ordered.files, read.rds.file)

  gene.name <- GENE.LOCATIONS[GENE.LOCATIONS$GeneId == gene.id, ]$Gene

  out.png.file <- paste0("report/timepoints/", species, ".", tissue, ".", gene.id, ".", gene.name, ".png")

  plots <- lapply(data, \(x) make.sashimi.panel(x, label = paste0(x$species, " ", gene.name, "\n", x$timepoint, " ", x$junction.strand))$plot)
  track <- make.gene.track(data[[1]])
  plots[[length(plots) + 1]] <- track

  patchwork::wrap_plots(plots, nrow = length(plots))
  save.double.width(filename = out.png.file, plot = last_plot(), height = 300)
}

#### Select groups for plotting ####

# find all species combinations for plotting

all.samples <- merge(make.sample.groups(), GENE.LOCATIONS, by = "CommonName")
cat("Plot sashimi: Selected samples for plotting\n")

species.groups <- all.samples %>%
  dplyr::group_by(CommonName, Timepoint, GeneId) %>%
  dplyr::summarise(SampleCount = n(), .groups = "drop_last")

tissue.groups <- all.samples %>%
  dplyr::group_by(Organism_part, Timepoint) %>%
  dplyr::summarise(SampleCount = n(), .groups = "drop_last")

timepoint.groups <- all.samples %>%
  dplyr::group_by(CommonName, Organism_part, GeneId) %>%
  dplyr::summarise(SampleCount = n(), .groups = "drop_last")

#### Functions to condense introns for neater plotting ####

# Make species plot using combined panels
make.condensed.species.panels <- function(species, timepoint, gene.id) {
  data.files <- list.files(path = "data/merged", pattern = paste0(species, ".*\\.", timepoint, "\\.", gene.id, "\\..*Rds_*"), full.names = TRUE)
  if (length(data.files) == 0) {
    return()
  }
  data <- lapply(data.files, read.rds.file)

  gene.name <- GENE.LOCATIONS[GENE.LOCATIONS$GeneId == gene.id, ]$Gene

  out.png.file <- paste0("report/species/", species, ".", timepoint, ".", gene.id, ".", gene.name, ".condensed.png")
  cat("Plot sashimi: Making", out.png.file, "\n")

  plots <- lapply(data, \(x)  make.gene.track.sashimi.panel(x, label = paste0(species, " ", gene.name, "\n", x$tissue, " ", x$junction.strand), show.x.axis = FALSE, is.collapse.introns = TRUE)$plot)

  patchwork::wrap_plots(plots, nrow = length(plots))

  save.double.width(filename = out.png.file, plot = last_plot(), height = 300)
}

# Make tissue plot using condensed panels
make.condensed.tissue.panels <- function(tissue, timepoint) {
  cat("Plot sashimi: Making", tissue, "at", timepoint, "\n")
  data.files <- list.files(path = "data/merged", pattern = paste0(".*\\.", tissue, "\\.", timepoint, "\\..*Rds_*"), full.names = TRUE)
  if (length(data.files) == 0) {
    return()
  }
  data <- lapply(data.files, read.rds.file)

  out.png.file <- paste0("report/tissues/", tissue, ".", timepoint, ".condensed.png")

  plots <- lapply(data, \(x) make.gene.track.sashimi.panel(x,
    label = paste0(x$species, "\n", x$gene.name, " ", x$junction.strand),
    show.x.axis = FALSE, is.collapse.introns = TRUE
  )$plot)

  patchwork::wrap_plots(plots, nrow = length(plots))
  save.double.width(filename = out.png.file, plot = last_plot(), height = 300)
}

# Make species plot using combined panels
make.condensed.timepoint.panels <- function(species, tissue, gene.id) {
  cat("Plot sashimi: Making", species, tissue, "for", gene.id, "\n")
  data.files <- list.files(path = "data/merged", pattern = paste0(species, ".*\\.", tissue, "\\..*", gene.id, "\\..*Rds_*"), full.names = TRUE)
  if (length(data.files) == 0) {
    return()
  }

  # Ensure files are plotted in time order
  ordered.files <- list()
  for (t in TIME.ORDER) {
    for (f in data.files) {
      if (str_detect(f, t)) {
        ordered.files <- c(ordered.files, f)
      }
    }
  }

  data <- lapply(ordered.files, read.rds.file)

  gene.name <- GENE.LOCATIONS[GENE.LOCATIONS$GeneId == gene.id, ]$Gene

  out.png.file <- paste0("report/timepoints/", species, ".", tissue, ".", gene.id, ".", gene.name, ".condensed.png")

  plots <- lapply(data, \(x) make.gene.track.sashimi.panel(x,
    label = paste0(x$species, " ", gene.name, "\n", x$timepoint, " ", x$junction.strand),
    show.x.axis = FALSE, is.collapse.introns = TRUE
  )$plot)

  patchwork::wrap_plots(plots, nrow = length(plots))
  save.double.width(filename = out.png.file, plot = last_plot(), height = 300)
}

#### Make condensed figures ####
cat("Plot sashimi: Making figures\n")
mapply(make.condensed.species.panels, species.groups$CommonName, species.groups$Timepoint, species.groups$GeneId)

# e.g. make.condensed.tissue.panels("forebrain", "adult")
mapply(make.condensed.tissue.panels, tissue.groups$Organism_part, tissue.groups$Timepoint)

# e.g. make.condensed.timepoint.panels("platypus", "testis", "ENSOANG00000046710")
make.condensed.timepoint.panels("anole", "testis", "ENSACAG00000007227")
mapply(make.condensed.timepoint.panels, timepoint.groups$CommonName, timepoint.groups$Organism_part, timepoint.groups$GeneId)


#### Make non-condensed figures ####

# Make species plots showing variation over tissues as a specific timepoint
mapply(make.species.panel, species.groups$CommonName, species.groups$Timepoint, species.groups$GeneId)

# Make tissue plots showing variation over species as a specific timepoint
mapply(make.tissue.panel, tissue.groups$Organism_part, tissue.groups$Timepoint)

# Make timepoint plots showing variation over times in a specific tissue
# e.g. make.timepoint.panel("chicken", "testis", "ENSGALG00010003052")
mapply(make.timepoint.panel, timepoint.groups$CommonName, timepoint.groups$Organism_part, timepoint.groups$GeneId)
cat("Plot sashimi: Done!\n")
