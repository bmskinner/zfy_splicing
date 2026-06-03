# Make sashimi plot based on ggsashimi.py modified to write data objects to Rds
# when run This custom sashimi plot ensures the transcripts are always left to
# right irrespective of strand.
source("src/functions.R")
source("src/ggsashimi.R")

cat("Plot sashimi: running shashimi plotting\n")

#### Custom plotting functions ####

#' Plot splice junctions on a gene exon track
#'
#' @param sashimi.data sashimi data as produced by read_sashimi_data
#' @param min.spanning.reads  the minimum number of spanning reads to include
#' @param label the label for the gene track
#' @param show.x.axis if true, display the x axis
#' @param is.collapse_introns if true, make each intron at most 500bp wide
#'
#' @returns
#' @export
#'
#' @examples
make_sashimi_coverage_plot <- function(sashimi.data, min.spanning.reads = 5, label = "label",
                                       show.x.axis = TRUE, is.collapse.introns = FALSE) {
  if (is.collapse.introns) {
    sashimi.data <- collapse_introns(
      sashimi.data,
      sashimi.data$reference.transcript.boundaries$exons,
      sashimi.data$reference.transcript.boundaries$introns
    )
  }

  junctions <- sashimi.data$junctions
  anns <- sashimi.data$reference.transcript.boundaries
  is.x.reverse <- sashimi.data$reference.transcript.strand == "-"

  if (any(!is.numeric(anns$exons$end)) | any(!is.numeric(anns$exons$start))) {
    cat("Error in annotations: at least one start or end is NA\n")
    print(anns$exons)
    str(anns$exons)
  }

  # Plot exons from the reference transcript
  reference.exons <- anns$exons |>
    dplyr::filter(
      transcript_id == sashimi.data$reference.transcript.id
    )

  non.reference.exons <- anns$exons |>
    dplyr::filter(
      transcript_id != sashimi.data$reference.transcript.id
    )

  reference.introns <- anns$introns |>
    dplyr::filter(
      start > min(reference.exons$start),
      end < max(reference.exons$start),
      transcript_id == sashimi.data$reference.transcript.id
    )

  non.reference.introns <- anns$introns |>
    dplyr::filter(
      start > min(reference.exons$start),
      end < max(reference.exons$start),
      transcript_id != sashimi.data$reference.transcript.id
    )

  # Set coordinate range for the x axis
  xmin <- min(c(reference.exons$start, reference.exons$end), na.rm = T) - 500
  xmax <- max(c(reference.exons$start, reference.exons$end), na.rm = T) + 500

  max.coverage <- max(c(sashimi.data$coverage$positive.strand, sashimi.data$coverage$negative.strand))

  intron.y <- 1.5
  exon.ymin <- 1.1
  exon.ymax <- 2

  # Make the gene track
  splot <- ggplot() +
    # Add read coverage, scaled to 0-1
    geom_area(data = sashimi.data$coverage, aes(x = position, y = positive.strand / max.coverage), fill = "darkgrey", col = "darkgrey") +
    geom_area(data = sashimi.data$coverage, aes(x = position, y = -negative.strand / max.coverage), fill = "darkgrey", col = "darkgrey") +

    # Non reference transcript introns
    geom_segment(
      data = non.reference.introns, aes(
        x = start,
        xend = end,
        y = ifelse(strand == sashimi.data$reference.transcript.strand, intron.y, -intron.y),
        yend = ifelse(strand == sashimi.data$reference.transcript.strand, intron.y, -intron.y)
      ),
      linewidth = 0.3, col = "grey"
    ) +
    # Reference transcript introns
    geom_segment(
      data = reference.introns, aes(
        x = start,
        xend = end,
        y = ifelse(strand == sashimi.data$reference.transcript.strand, intron.y, -intron.y),
        yend = ifelse(strand == sashimi.data$reference.transcript.strand, intron.y, -intron.y)
      ),
      linewidth = 0.3, col = "blue"
    ) +

    # Add upper strand label
    geom_hline(yintercept = 0) +
    annotate("text",
      x = ifelse(is.x.reverse, Inf, -Inf),
      y = intron.y + 1,
      label = sashimi.data$reference.transcript.strand,
      hjust = 0, vjust = 0.5
    ) +

    # Add lower strand label
    annotate("text",
      x = ifelse(is.x.reverse, Inf, -Inf),
      y = -intron.y - 1,
      label = ifelse(sashimi.data$reference.transcript.strand == "+", "-", "+"),
      hjust = 0, vjust = 0.5
    ) +

    # Non-reference transcript exons
    geom_rect(
      data = non.reference.exons, aes(
        xmin = start,
        xmax = end,
        ymin = ifelse(strand == sashimi.data$reference.transcript.strand, exon.ymin, -exon.ymin),
        ymax = ifelse(strand == sashimi.data$reference.transcript.strand, exon.ymax, -exon.ymax)
      ),
      fill = "grey", alpha = 1
    ) +

    # Reference transcript exons
    geom_rect(
      data = reference.exons, aes(
        xmin = start,
        xmax = end,
        ymin = ifelse(strand == sashimi.data$reference.transcript.strand, exon.ymin, -exon.ymin),
        ymax = ifelse(strand == sashimi.data$reference.transcript.strand, exon.ymax, -exon.ymax)
      ),
      fill = "blue", alpha = 1
    ) +
    scale_y_discrete(expand = c(0.35, 0.35)) +
    coord_cartesian(xlim = c(xmin, xmax), ylim = c(-3, 3)) +
    scale_x_continuous(expand = c(0.05, 0.05)) +
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
      scale_x_reverse(expand = c(0.05, 0.05))
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

    label.y.offset <- 0.4 # separation between spline and label

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
        ymax = -ymax
      )
    } else {
      splot <- splot + annotation_custom(
        grob = l.grob.top,
        xmin = xmin,
        xmax = xmid,
        ymax = ymin,
        ymin = ymax
      )
    }

    # Right arc
    if (is.above) { # Junctions below zero
      splot <- splot + annotation_custom(
        grob = r.grob.btm,
        xmin = xmid,
        xmax = xmax,
        ymin = -ymin,
        ymax = -ymax
      )
    } else {
      splot <- splot + annotation_custom(
        grob = r.grob.top,
        xmin = xmid,
        xmax = xmax,
        ymax = ymin,
        ymin = ymax
      )
    }

    splot <- splot + annotate("label",
      x = xmid,
      y = ifelse(is.above, 0 - ymax - label.y.offset, ymax + label.y.offset),
      label = as.character(count),
      size = 2, col = spline.color, fill = "white"
    )

    return(splot)
  }

  # Check junctions meet plot criteria
  if (nrow(junctions) > 0) {
    junctions <- junctions %>%
      dplyr::filter(count >= min.spanning.reads) %>%
      na.omit() %>%
      dplyr::arrange(start, end)


    if (any(!is.numeric(junctions$end)) | any(!is.numeric(junctions$start))) {
      print(junctions)
      str(junctions)
      stop("Error in junctions: at least one x or xend is NA")
    }

    junctions <- junctions %>%
      dplyr::mutate(length = end - start + 1)
  }

  # Add the junctions, adjusting for plus vs minus strand
  if (nrow(junctions) > 0) {
    # Calculate charting coordinates
    junctions$y.offset <- rep(seq(0, 1.5, 0.5), length.out = nrow(junctions)) # give each junction a separate y offset

    for (i in 1:nrow(junctions)) {
      jrow <- junctions[i, ]

      splot <- add.junction(
        splot = splot,
        xmin = jrow$start, xmax = jrow$end,
        is.above = jrow$strand == sashimi.data$reference.transcript.strand,
        ymin = 2, ymax = 2.5 + jrow$y.offset,
        count = jrow$count,
        matches.known.exons = jrow$matches.known.exon
      )
    }
  }

  list(plot = splot, sashimi.data = sashimi.data)
}

#### Main script ####

# Ensure output dirs exist
fs::dir_create(c("report/species", "report/timepoints", "report/tissues", "report/raw_sashimi"))

# Read all GTF files once, since we have multiple genes/tissues per species
GTF.DATA <- read_gtf_data(GENOME.DATA$GTF_FILE, GENOME.DATA$CommonName)

cat("Plot sashimi: Making figures\n")

bam.files <- data.frame(path = list.files(path = "data/stringtie", pattern = ".*.bam", full.names = TRUE)) |>
  dplyr::mutate(file = basename(path)) |>
  tidyr::separate_wider_delim(file,
    delim = ".", names = c("species", "tissue", "timepoint", "gene_id", "gene_name", "ext"),
    too_few = "debug", too_many = "debug"
  )

for (i in 1:nrow(bam.files)) {
  bam.row <- bam.files[i, ]
  species <- bam.row$species
  tissue <- bam.row$tissue
  timepoint <- bam.row$timepoint
  gene_id <- bam.row$gene_id
  gene_name <- bam.row$gene_name

  cat("Detecting splice junctions for", species, tissue, timepoint, gene_id, "\n")

  gene.data <- GENE.LOCATIONS[GENE.LOCATIONS$GeneId == gene_id, ]
  coords <- parse_coordinates(gene.data$Location)

  tryCatch(
    {
      sashimi.data <- read_sashimi_data(
        bam.file = bam.row$path,
        gtf.data = GTF.DATA[[species]],
        chr = coords$coord.chr, start = coords$coord.start, end = coords$coord.end,
        reference.gene.id = gene.data$GeneId,
        reference.transcript.id = gene.data$CanonicalTranscriptId
      )

      # Create plot with collapsed introns
      sashimi.plot.collapsed <- make_sashimi_coverage_plot(sashimi.data,
        is.collapse.introns = TRUE, show.x.axis = FALSE,
        min.spanning.reads = 5, label = paste0(species, "\n", tissue, "\n", timepoint, "\n", gene_name)
      )

      save.double.width(paste0("report/raw_sashimi/", species, ".", tissue, ".", timepoint, ".", gene_id, ".", gene_name, ".condensed.png"),
        sashimi.plot.collapsed$plot,
        height = 50
      )
      # Create plot with expanded introns
      sashimi.plot.expanded <- make_sashimi_coverage_plot(sashimi.data,
        is.collapse.introns = FALSE,
        min.spanning.reads = 5, label = paste0(species, "\n", tissue, "\n", timepoint, "\n", gene_name)
      )

      save.double.width(paste0("report/raw_sashimi/", species, ".", tissue, ".", timepoint, ".", gene_id, ".", gene_name, ".expanded.png"),
        sashimi.plot.expanded$plot,
        height = 50
      )

      # Create a plot with coverage info. Only expanded available
      # coverage.plot <- make_coverage_plot(sashimi.data, label = paste0(species, "\n", tissue, "\n", timepoint, "\n", gene_name))
      # save.double.width(paste0("report/raw_sashimi/", species, ".", tissue, ".", timepoint, ".", gene_id, ".", gene_name, ".coverage.png"),
      #   coverage.plot,
      #   height = 80
      # )
    },
    error = \(e) {
      cat("Error plotting sashimi data from", bam.row$path, "\n", paste(e))
      e
    }
  )
}

cat("Plot sashimi: Done!\n")
