# R implementation of ggsashimi.py from
# https://github.com/guigolab/ggsashimi

SAM.FLAG.READ.UNMAPPED <- 0x4
SAM.FLAG.MATE.UNMAPPED <- 0x8
SAM.FLAG.READ.REVERSE.STRAND <- 0x10
SAM.FLAG.FIRST.IN.PAIR <- 0x40
SAM.FLAG.SECOND.IN.PAIR <- 0x80

#' Create a coverage map object
#'
#' @param start the starting base of the region of interest in the genome
#' @param end the end base of the region of interest in the genome
#'
#' @returns
#'
#' @examples
make.coverage.map <- function(start, end) {
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
make.junction.map <- function() {
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
has.flag <- function(sam.flag, property) {
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
#' parse.coordinates("Y:1234-5678")
parse.coordinates <- function(coordinate.string) {
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
#     if (has.flag(samflag, SAM.FLAG.FIRST.IN.PAIR)) {
#       return(0)
#     }
#     if (has.flag(samflag, SAM.FLAG.SECOND.IN.PAIR)) {
#       return(1)
#     }
#   }
#   if (strand.string == "MATE2_SENSE") {
#     if (has.flag(samflag, SAM.FLAG.FIRST.IN.PAIR)) {
#       return(1)
#     }
#     if (has.flag(samflag, SAM.FLAG.SECOND.IN.PAIR)) {
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

  is.reverse.strand <- has.flag(samflag, SAM.FLAG.READ.REVERSE.STRAND)
  if (strand.string == "SENSE") {
    return(ifelse(is.reverse.strand, "-", "+"))
  }

  if (strand.string == "ANTISENSE") {
    return(ifelse(is.reverse.strand, "+", "-"))
  }
  if (strand.string == "MATE1_SENSE") {
    # 64 = first in pair, 128 = second in pair
    if (has.flag(samflag, SAM.FLAG.FIRST.IN.PAIR)) {
      return(ifelse(is.reverse.strand, "-", "+"))
    }
    if (has.flag(samflag, SAM.FLAG.SECOND.IN.PAIR)) {
      return(ifelse(is.reverse.strand, "+", "-"))
    }
  }
  if (strand.string == "MATE2_SENSE") {
    if (has.flag(samflag, SAM.FLAG.FIRST.IN.PAIR)) {
      return(ifelse(is.reverse.strand, "+", "-"))
    }
    if (has.flag(samflag, SAM.FLAG.SECOND.IN.PAIR)) {
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
  cat("Reading bam file:", bam.file, "\n")
  coordinates <- parse.coordinates(coordinate.string)

  # Initialize empty coverage array and junction maps
  coverage.array <- make.coverage.map(coordinates$coord.start, coordinates$coord.end)
  junction.map <- make.junction.map()

  bam.conn <- Rsamtools::BamFile(bam.file)
  # bam.data <- Rsamtools::BamFile(bam.file, index = paste0(bam.file, ".csi"))
  bam.data <- scanBam(bam.conn)[[1]]

  # Go read by read
  for (i in 1:length(bam.data$qname)) {
    if (i %% 500 == 0) cat(sprintf("Processed %.2f%% of reads\n", i / length(bam.data$qname) * 100))
    read.data <- lapply(bam.data, function(xx) xx[i])
    # Skip if read is unmapped
    if (has.flag(read.data$flag, SAM.FLAG.READ.UNMAPPED) |
      has.flag(read.data$flag, SAM.FLAG.MATE.UNMAPPED)) {
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
get_exon_boundaries <- function(gtf.data, chr, start, end) {
  region <- gtf.data[gtf.data$seqid == chr & gtf.data$start >= start - 5000 & gtf.data$end <= end + 5000, ]
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
    dplyr::select(gene_id, transcript_id, gene_name, transcript_name, type,
      start = intron.start, end = intron.end, strand, length
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(transcript_id, start, end)

  region.exons <- region |>
    as.data.frame() |>
    dplyr::filter(type == "exon") |>
    dplyr::select(
      gene_id, transcript_id, gene_name, transcript_name, type,
      start, end, strand
    ) |>
    dplyr::mutate(length = end - start + 1) |>
    dplyr::arrange(transcript_id, start, end)

  list(
    exons = region.exons,
    introns = region.introns
  )
}


#' Translate junction count formats to the original ggsashimi objects for plotting.
#' Temp function - migrate the plot function to use the original data structures
#'
#' @param bam.file
#' @param gtf.data
#' @param chr
#' @param start
#' @param end
#'
#' @returns
#' @export
#'
#' @examples
make_ggasashimi_output <- function(bam.file, gtf.data, chr, start, end) {
  result <- read_bam(bam.file, paste0(chr, ":", start, "-", end), "SENSE")
  boundaries <- get_exon_boundaries(gtf.data, chr, start, end)

  ann_list <- list()
  ann_list$exons <- boundaries$exons |>
    dplyr::select(tx = transcript_id, start, end, strand)
  ann_list$introns <- boundaries$introns |>
    dplyr::select(tx = transcript_id, start, end, strand)

  density_list <- result$coverage |>
    dplyr::mutate(position = as.integer(position)) |>
    dplyr::select(x = position, y = coverage)

  junction_list <- result$junctions |>
    dplyr::select(
      x = start, xend = end, y = value, yend = value, count = value
    ) |>
    dplyr::mutate(matches.known.exons = FALSE)

  list(
    ann_list = ann_list,
    density_list = density_list,
    junction_list = junction_list
  )
}
