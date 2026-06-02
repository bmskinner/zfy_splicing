# R implementation of ggsashimi.py from
# https://github.com/guigolab/ggsashimi

SAM.FLAG.READ.UNMAPPED <- 0x4
SAM.FLAG.MATE.UNMAPPED <- 0x8
SAM.FLAG.READ.REVERSE.STRAND <- 0x10
SAM.FLAG.FIRST.IN.PAIR <- 0x40
SAM.FLAG.SECOND.IN.PAIR <- 0x80

#### Objects for storing junction information ####

#' Create a coverage map object
#'
#' @param start the starting base of the region of interest in the genome
#' @param end the end base of the region of interest in the genome
#'
#' @returns
#'
#' @examples
make_coverage_map <- function(start, end) {
  coverage.array <- new.env()

  n.indexes <- end - start + 1

  coverage.array$"+" <- rep(0, n.indexes)
  names(coverage.array$`+`) <- start:end

  coverage.array$"-" <- rep(0, n.indexes)
  names(coverage.array$`-`) <- start:end

  coverage.array$add <- function(pos, strand, value) {
    coverage.array[[strand]][as.character(pos)] <- coverage.array[[strand]][as.character(pos)] + value
  }

  coverage.array$increment <- function(pos, strand) {
    coverage.array$add(pos, strand, 1)
  }

  coverage.array$incrementRange <- function(start, end, strand) {
    # coverage.array$add(pos, strand, 1)
    base.range <- as.character(start:end)
    coverage.array[[strand]][base.range] <- coverage.array[[strand]][base.range] + 1
  }

  coverage.array$get <- function(strand) {
    pos <- data.frame(
      position = names(coverage.array[["+"]]),
      positive.strand = coverage.array[["+"]]
    )

    neg <- data.frame(
      position = names(coverage.array[["-"]]),
      negative.strand = coverage.array[["-"]]
    )

    merge(pos, neg, by = "position", all = TRUE) |>
      dplyr::mutate(coverage = positive.strand + negative.strand)
  }
  coverage.array
}

#' Create a map for splice junctions
#'
#' @returns
#'
#' @examples
make_junction_map <- function() {
  junction.map <- new.env()
  junction.map$add <- function(start, end, strand, value) {
    if (rlang::env_has(junction.map, "values")) {
      existing <- rlang::env_get(junction.map, "values")
      if (any(existing$start == start & existing$end == end & existing$strand == strand)) {
        row.id <- which(existing$start == start & existing$end == end & existing$strand == strand)
        existing[row.id, "count"] <- existing[row.id, "count"] + value
        rlang::env_poke(junction.map, "values", existing)
      } else {
        df <- rbind(existing, data.frame(start = start, end = end, strand = strand, count = value))
        rlang::env_poke(junction.map, "values", df)
      }
    } else { # first iteration, create the df
      rlang::env_poke(
        junction.map, "values",
        data.frame(start = start, end = end, strand = strand, count = value)
      )
    }
  }
  junction.map
}

#### Core functions ####

#' Check for a SAM flag.
#'
#' SAM flags are an integer, with each bit representing a property. This
#' function uses a bitwise AND to check for the existence of the property.
#'
#' @param sam.flag the read flag to check
#' @param property the SAM flag property to look for
#'
#' @returns true if the SAM property is present within sam.flag, false otherwise
#'
#' @examples
has_sam_flag <- function(sam.flag, property) {
  as.integer(sam.flag) %&% property == property
}


#' Parse a genomic coordinate string. The string should have the format
#' "chr:start-end"
#'
#' @param coordinate.string the coordinate string
#'
#' @returns a list with chromosome, start and end elements
#'
#' @examples
#' parse_coordinates("Y:1234-5678")
parse_coordinates <- function(coordinate.string) {
  coordinate.string <- stringr::str_remove_all(coordinate.string, ",")

  coordinate.elements <- stringr::str_split(coordinate.string, ":")

  coord.chr <- coordinate.elements[[1]][1]

  coord.locations <- stringr::str_split(coordinate.elements[[1]][2], "-")
  coord.start <- as.integer(coord.locations[[1]][1])
  coord.end <- as.integer(coord.locations[[1]][2])

  list(
    coord.chr = coord.chr,
    coord.start = coord.start,
    coord.end = coord.end
  )
}

#' Determine if a read is sense or antisense from SAM flag and strand string.
#'
#' @param strand.string the strand of the read
#' @param samflag integer SAM flag
#'
#' @returns 0 if the read is sense, 1 if the read is antisense
#' @export
#'
#' @examples
# flip_read <- function(strand.string, samflag) {
#   if (strand.string == "NONE" | strand.string == "SENSE") {
#     return(0)
#   }
#   if (strand.string == "ANTISENSE") {
#     return(1)
#   }
#   if (strand.string == "MATE1_SENSE") {
#     # 64 = first in pair, 128 = second in pair
#     if (has_sam_flag(samflag, SAM.FLAG.FIRST.IN.PAIR)) {
#       return(0)
#     }
#     if (has_sam_flag(samflag, SAM.FLAG.SECOND.IN.PAIR)) {
#       return(1)
#     }
#   }
#   if (strand.string == "MATE2_SENSE") {
#     if (has_sam_flag(samflag, SAM.FLAG.FIRST.IN.PAIR)) {
#       return(1)
#     }
#     if (has_sam_flag(samflag, SAM.FLAG.SECOND.IN.PAIR)) {
#       return(0)
#     }
#   }
# }

#' Determine if a read is sense or antisense from SAM flag relative to a
#' reference strand.
#'
#' @param strand.string the reference strand to test against
#' @param samflag integer SAM flag
#'
#' @returns + if the read is sense relative to the reference strand, - otherwise
#'
#' @examples
#' find_read_strand("SENSE", 42)
find_read_strand <- function(strand.string, samflag) {
  if (strand.string == "NONE") {
    return("+")
  }

  is.reverse.strand <- has_sam_flag(samflag, SAM.FLAG.READ.REVERSE.STRAND)
  if (strand.string == "SENSE") {
    return(ifelse(is.reverse.strand, "-", "+"))
  }

  if (strand.string == "ANTISENSE") {
    return(ifelse(is.reverse.strand, "+", "-"))
  }
  if (strand.string == "MATE1_SENSE") {
    # 64 = first in pair, 128 = second in pair
    if (has_sam_flag(samflag, SAM.FLAG.FIRST.IN.PAIR)) {
      return(ifelse(is.reverse.strand, "-", "+"))
    }
    if (has_sam_flag(samflag, SAM.FLAG.SECOND.IN.PAIR)) {
      return(ifelse(is.reverse.strand, "+", "-"))
    }
  }
  if (strand.string == "MATE2_SENSE") {
    if (has_sam_flag(samflag, SAM.FLAG.FIRST.IN.PAIR)) {
      return(ifelse(is.reverse.strand, "+", "-"))
    }
    if (has_sam_flag(samflag, SAM.FLAG.SECOND.IN.PAIR)) {
      return(ifelse(is.reverse.strand, "-", "+"))
    }
  }
}


#' Read a bam file. Count per-base coverage and detect splice junctions.
#'
#' @param bam.file the bam file to read
#' @param coordinate.string the location in which to count e.g. 1:1234-5678
#' @param strand.string the strand to test relative to. One of NONE, SENSE, ANTISENSE,
#'   MATE1_SENSE, MATE2_SENSE
#'
#' @returns a list with a read depth coverage dataframe and a junctions
#'   dataframe for the region
#' @export
#'
#' @examples
#' read_bam("/path/to/reads.bam", "1:118296484-118319811", "NONE")
read_bam <- function(bam.file, coordinate.string, strand.string) {
  coordinates <- parse_coordinates(coordinate.string)
  cat("Reading bam file:", bam.file, "at", coordinates$coord.chr, ":", coordinates$coord.start, "-", coordinates$coord.end, "\n")
  # Initialize empty coverage array and junction maps
  coverage.array <- make_coverage_map(coordinates$coord.start, coordinates$coord.end)
  junction.map <- make_junction_map()

  bam.conn <- Rsamtools::BamFile(bam.file)
  bam.data <- scanBam(bam.conn)[[1]]

  # Go read by read
  for (i in 1:length(bam.data$qname)) {
    if (i %% 500 == 0) cat(sprintf("Processed %.2f%% of reads\n", i / length(bam.data$qname) * 100))
    read.data <- lapply(bam.data, function(xx) xx[i])
    # Skip if read is unmapped
    if (has_sam_flag(read.data$flag, SAM.FLAG.READ.UNMAPPED) |
      has_sam_flag(read.data$flag, SAM.FLAG.MATE.UNMAPPED)) {
      next
    }

    # Ignore reads with more exotic CIGAR operators
    # print(read.data$cigar)
    if (any(stringr::str_detect(read.data$cigar, c("H", "P", "X", "=")))) next

    # Determine the strand of the read
    read.strand <- find_read_strand(strand.string, read.data$flag)

    # Parse the cigar string
    # Will be e.g. 10M3I21D would parse into 10, 3, 21 and M, I, D
    cigar.lengths <- stringr::str_extract_all(read.data$cigar, "\\d+")[[1]]
    cigar.ops <- stringr::str_extract_all(read.data$cigar, "[MIDNS]")[[1]]

    current.position <- read.data$pos



    #' Count coverage and splice junctions from a cigar string
    #'
    #' @param op the cigar op code
    #' @param len the cigar length
    #' @param pos the position of the pointer in the genome
    #' @param strand the read strand
    #'
    #' @returns
    #' @export
    #'
    #' @examples
    count.operator <- function(op, len, pos, strand) {
      # Match - increase coverage across the range
      if (op == "M") {
        pos.range <- pos:(pos + len - 1)
        pos.range <- pos.range[pos.range >= coordinates$coord.start & pos.range < coordinates$coord.end]
        if (length(pos.range) > 0) {
          coverage.array$incrementRange(pos.range[1], pos.range[length(pos.range)], strand)
        }
      }

      # Insertion or Soft-clip
      # we do not change position
      if (op == "I" | op == "S") {
        return(pos)
      }

      # Deletion
      if (op == "D") {
        # no action
      }

      # Junction
      if (op == "N") {
        don <- pos # splice donor
        acc <- pos + len # splice acceptor - somewhere downstream
        if (don >= coordinates$coord.start & acc <= coordinates$coord.end) {
          junction.map$add(don, acc, strand, 1)
        }
      }

      pos <- pos + len

      return(pos)
    }

    # Check all cigar ops for this read.
    for (i in 1:length(cigar.ops)) {
      curr.cigar.length <- as.integer(cigar.lengths[i])
      curr.cigar.op <- cigar.ops[i]
      current.position <- count.operator(
        curr.cigar.op, curr.cigar.length, current.position, read.strand
      )
    }
  }

  coverage <- coverage.array$get() |>
    dplyr::mutate(position = as.integer(position))

  return(list(
    coverage = coverage,
    junctions = junction.map$values,
    reference.strand = strand.string,
    reference.location = coordinate.string
  ))
}



#' Test if the given transcript is on the forward or reverse strand
#'
#' @param transcript.id the transcript to test
#' @param gtf.data the GTF data for the genome as GenomicRanges e.g. as read by rtracklayer
#'
#' @returns true if any exons of the given transcript are on the reverse strand, false otherwise
#' @export
#'
#' @examples
transcript_is_reverse_strand <- function(transcript.id, gtf.data) {
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
get_transcript_junctions <- function(transcript.id, gtf.data) {
  # Find the exons for this transcript
  exons <- gtf.data[gtf.data$type == "exon" & gtf.data$transcript_id == transcript.id, ]

  if (transcript_is_reverse_strand(transcript.id, gtf.data)) {
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
#' junction_is_in_GTF("ENSGALT00010007119", chicken.gtf, 118299992, 118301364)
junction_is_in_GTF <- function(transcript.id, gtf.data, junction.start, junction.end) {
  junctions <- get_transcript_junctions(transcript.id, gtf.data)
  any(junctions$j1 == junction.start & junctions$j2 == junction.end)
}

#' Get intron and exon boundaries for all transcripts within the given region of
#' a GTF
#'
#' @param gtf.data the GTF data as read by e.g. rtracklayer
#' @param chr the chromosome
#' @param start the start of the window
#' @param end the end of the window
#'
#' @returns a data frame of exon and intron locations within the window
#' @export
#'
#' @examples
get_exon_boundaries <- function(gtf.data, chr, loc.start, loc.end) {
  cat("Getting exon boundaries\n")

  region <- gtf.data |>
    dplyr::filter(
      seqnames == as.character(chr),
      start >= loc.start - 5000 &
        end <= loc.end + 5000
    )

  cat("Region ", chr, ":", loc.start, "-", loc.end, "contains", nrow(region), "data rows\n")

  region.introns <- region |>
    as.data.frame() |>
    dplyr::filter(type == "exon") |>
    dplyr::arrange(start, end) |>
    dplyr::group_by(transcript_id) |>
    dplyr::mutate(
      intron.start = end + 1,
      intron.end = dplyr::lead(start - 1),
      length = intron.end - intron.start + 1,
      type = "intron"
    ) |>
    dplyr::filter(!is.na(intron.start), !is.na(intron.end)) |>
    dplyr::select(seqnames, gene_id, transcript_id, gene_name, transcript_name, type,
      start = intron.start, end = intron.end, strand, length
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(transcript_id, start, end)

  cat("Region ", chr, ":", loc.start, "-", loc.end, "contains", nrow(region.introns), "intron rows\n")

  region.exons <- region |>
    as.data.frame() |>
    dplyr::filter(type == "exon") |>
    dplyr::select(
      seqnames, gene_id, transcript_id, gene_name, transcript_name, type,
      start, end, strand
    ) |>
    dplyr::mutate(length = end - start + 1) |>
    dplyr::arrange(transcript_id, start, end)

  cat("Region ", chr, ":", loc.start, "-", loc.end, "contains", nrow(region.exons), "exon rows\n")

  list(
    exons = region.exons,
    introns = region.introns
  )
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
read_sashimi_data <- function(bam.file, gtf.data, chr, start, end,
                              reference.gene.id, reference.transcript.id = NA) {
  if (is.null(bam.file)) stop("No bam file specified")
  if (!file.exists(bam.file)) stop(paste("Cannot find the bam file:", bam.file))
  if (is.null(gtf.data)) stop("No gtf data given")
  if (is.null(chr)) stop("No chromosome given")
  if (is.null(start)) stop("No start coordinate given")
  if (is.null(end)) stop("No end coordinate given")
  if (is.null(reference.gene.id)) stop("No reference gene id given")

  tryCatch(
    {
      cat("Reading sashimi data from bam file in region '", paste0(chr, ":", start, "-", end), "'\n")
      sashimi.data <- list()
      sashimi.data$input.file <- bam.file
      bam.data <- read_bam(bam.file, paste0(chr, ":", start, "-", end), "SENSE")

      sashimi.data$reference.gene.id <- reference.gene.id
      # Get the longest transcript in the gene if none specified
      if (is.na(reference.transcript.id)) {
        cat("No reference transcript given, selecting longest for gene id", reference.gene.id, "\n")
        reference.transcript.id <- gtf.data |>
          dplyr::filter(gene_id == reference.gene.id, type == "transcript") |>
          dplyr::mutate(length = end - start + 1) |>
          dplyr::arrange(length) |>
          dplyr::slice_tail(n = 1) |>
          dplyr::select(transcript_id) |>
          dplyr::pull()
      }
      sashimi.data$reference.transcript.id <- reference.transcript.id
      cat("Reference transcript is", reference.transcript.id, "\n")

      if (nrow(gtf.data[gtf.data$gene_id == reference.gene.id, ]) == 0) {
        stop("Unable to detect reference gene id in genome GTF\n")
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

      cat("Reference transcript is on strand '", paste(sashimi.data$reference.transcript.strand, collapse = ","), "'\n")

      sashimi.data$reference.transcript.boundaries <- get_exon_boundaries(gtf.data, chr, start, end)

      cat("Reference exon/intron bounds detected:\n")

      sashimi.data$reference.gene.name <- gtf.data |>
        dplyr::filter(gene_id == reference.gene.id) |>
        dplyr::select(gene_name) |>
        dplyr::distinct() |>
        dplyr::pull()

      sashimi.data$coverage <- bam.data$coverage
      sashimi.data$junctions <- bam.data$junctions
      sashimi.data$junctions.reference.strand <- bam.data$reference.strand

      sashimi.data$junctions$matches.known.exon <- mapply(junction_is_in_GTF,
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
      cat("Error reading sashimi data from", bam.file, "\n", paste(e))
      e
    }
  )
}



#' Make a conversion table that can rescale coordinates to collapse introns to a
#' maximum length. This allows plots to show exon junctions more clearly in
#' genes with long introns.
#'
#' @param exon.data data frame of exon coordinates with columns start, end and
#'   strand
#' @param intron.data data frame of intron coordinates with columns start, end
#'   and strand
#' @param strand the strand with the transcript to convert against
#' @param max.intron.length the maximum length for an intron
#'
#' @returns a function that will convert coordinates collapsing introns
#' @export
#'
#' @examples
create_intron_collapser <- function(exon.data, intron.data, strand, max.intron.length = 500) {
  if (is.null(exon.data)) stop("No exon data provided")

  cat("Creating intron collapser\n")
  # Reduce any overlapping exons if we have multiple transcripts
  exon.data <- exon.data[exon.data$strand == strand, ]
  intron.data <- intron.data[intron.data$strand == strand, ]

  cat("Detected", nrow(exon.data), "exons and", nrow(intron.data), "introns\n")

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
      warning("No range table entry covering", coordinate)
      return(coordinate)
    }

    # How far along the feature are we?
    fractional.distance <- (coordinate - feature$start) / feature$original.length

    # Nearest integer to the same fraction of the new coordinate space
    result <- (fractional.distance * feature$new.length) + feature$new.start
    return(result)
  }
  cat("Created intron collapser\n")
  list(
    calculate = calculate,
    full.ranges = full.ranges
  )
}


#' Collapse introns to a maximum length for pretty plotting of sashimi splice
#' data.
#'
#' Note that the locations of the introns and exons are no longer accurate
#' against the reference genome. Chart axes will not be meaningful.
#'
#' @param sashimi.data the sashimi data as created by read_sashimi_data
#' @param exon.data the exon annotation in the region
#' @param intron.data the intron annotations in the region
#'
#' @returns the sashimi data object with collapsed intron locations added.
#' @export
#'
#' @examples
collapse_introns <- function(sashimi.data, exon.data, intron.data) {
  # Calculate offsets to make all introns at most 500bp
  intron.collapser <- create_intron_collapser(exon.data, intron.data, sashimi.data$reference.transcript.strand)

  cat("Collapsing introns\n")
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

#### Charting functions ####

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
make_sashimi_plot <- function(sashimi.data, min.spanning.reads = 5, label = "tissue",
                              show.x.axis = TRUE, is.collapse_introns = FALSE) {
  if (is.collapse_introns) {
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
    annotate("text",
      x = ifelse(is.x.reverse, xmax - 250, xmin + 250),
      y = 0.5,
      label = sashimi.data$reference.transcript.strand,
      size = 5, col = "black"
    ) +
    annotate("text",
      x = ifelse(is.x.reverse, xmax - 250, xmin + 250),
      y = -0.5,
      label = ifelse(sashimi.data$reference.transcript.strand == "+", "-", "+"),
      size = 5, col = "black"
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

    label.y.offset <- 0.2 # separation between spline and label

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
      size = 2, col = spline.color, fill = "white", label.size = NA
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
    junctions$y.offset <- rep(seq(0, 1, 0.5), length.out = nrow(junctions)) # give each junction a separate y offset

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
