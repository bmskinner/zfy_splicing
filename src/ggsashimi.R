# R implementation of ggsashimi.py from
# https://github.com/guigolab/ggsashimi

packages <- c(
  "parallel", "tidyverse", "GenomicRanges", "grid",
  "fs", "data.table", "rtracklayer", "Rsamtools", "bitops", "rlang",
  "data.table"
)

suppressPackageStartupMessages({
  is.installed <- sapply(packages, require, character.only = TRUE)
  if (!all(is.installed)) stop("The following packages are required:", paste(packages[!is.installed], collapse = ", "))
})


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

  coverage.array$"FORWARD" <- rep(0, n.indexes)
  names(coverage.array$FORWARD) <- start:end

  coverage.array$"REVERSE" <- rep(0, n.indexes)
  names(coverage.array$REVERSE) <- start:end

  coverage.array$add <- function(pos, strand, value) {
    coverage.array[[strand]][as.character(pos)] <- coverage.array[[strand]][as.character(pos)] + value
  }

  coverage.array$increment <- function(pos, strand) {
    coverage.array$add(pos, strand, 1)
  }

  coverage.array$incrementRange <- function(start, end, strand) {
    base.range <- as.character(start:end)
    coverage.array[[strand]][base.range] <- coverage.array[[strand]][base.range] + 1
  }

  coverage.array$get <- function(strand) {
    pos <- data.frame(
      position = names(coverage.array[["FORWARD"]]),
      forward.strand = coverage.array[["FORWARD"]]
    )

    neg <- data.frame(
      position = names(coverage.array[["REVERSE"]]),
      reverse.strand = coverage.array[["REVERSE"]]
    )

    merge(pos, neg, by = "position", all = TRUE) |>
      dplyr::mutate(coverage = forward.strand + reverse.strand)
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
  # An empty data frame that will be overwritten if there are any reads
  rlang::env_poke(
    junction.map, "values",
    data.frame(start = c(), end = c(), strand = c(), count = c())
  )

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
    return("FORWARD")
  }

  is.reverse.strand <- has_sam_flag(samflag, SAM.FLAG.READ.REVERSE.STRAND)
  if (strand.string == "SENSE") {
    return(ifelse(is.reverse.strand, "REVERSE", "FORWARD"))
  }

  if (strand.string == "ANTISENSE") {
    return(ifelse(is.reverse.strand, "FORWARD", "REVERSE"))
  }
  if (strand.string == "MATE1_SENSE") {
    # 64 = first in pair, 128 = second in pair
    if (has_sam_flag(samflag, SAM.FLAG.FIRST.IN.PAIR)) {
      return(ifelse(is.reverse.strand, "REVERSE", "FORWARD"))
    }
    if (has_sam_flag(samflag, SAM.FLAG.SECOND.IN.PAIR)) {
      return(ifelse(is.reverse.strand, "FORWARD", "REVERSE"))
    }
  }
  if (strand.string == "MATE2_SENSE") {
    if (has_sam_flag(samflag, SAM.FLAG.FIRST.IN.PAIR)) {
      return(ifelse(is.reverse.strand, "FORWARD", "REVERSE"))
    }
    if (has_sam_flag(samflag, SAM.FLAG.SECOND.IN.PAIR)) {
      return(ifelse(is.reverse.strand, "REVERSE", "FORWARD"))
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

  # Use an index if available - otherwise read directly but slower
  bai.index.file <- paste0(bam.file, ".bai")
  csi.index.file <- paste0(bam.file, ".csi")
  if (file.exists(bai.index.file)) {
    bam.conn <- Rsamtools::BamFile(bam.file, bai.index.file)
  } else if (file.exists(csi.index.file)) {
    bam.conn <- Rsamtools::BamFile(bam.file, csi.index.file)
  } else {
    bam.conn <- Rsamtools::BamFile(bam.file)
  }
  bam.data <- scanBam(bam.conn)[[1]]

  if (length(bam.data$qname) == 0) {
    stop("There are no reads in bam file", bam.file)
  } else {
    cat("There are", length(bam.data$qname), "reads in the bam file\n")
  }

  # Go read by read
  for (i in 1:length(bam.data$qname)) {
    if (i %% 500 == 0) cat(sprintf("Processed %.2f%% of %s reads\n", i / length(bam.data$qname) * 100, length(bam.data$qname)))
    read.data <- lapply(bam.data, function(xx) xx[i])
    # Skip if read is unmapped
    if (has_sam_flag(read.data$flag, SAM.FLAG.READ.UNMAPPED) |
      has_sam_flag(read.data$flag, SAM.FLAG.MATE.UNMAPPED)) {
      next
    }

    # Ignore reads with more exotic CIGAR operators
    # print(read.data$cigar)
    if (any(stringr::str_detect(read.data$cigar, c("H", "P", "X", "=")))) next

    # Determine the strand of this read in the genome
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
        don <- pos - 1 # splice donor TODO
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

#' Given a transcript id, find the splice junctions within that transcript
#'
#' @param transcript.id the transcript to test
#' @param gtf.data the GTF data for the genome as GenomicRanges e.g. as read by rtracklayer
#'
#' @returns
#' @export
#'
#' @examples
get_transcript_junctions <- function(transcript.id, gtf.data) {
  # Junction coordinates should be final base of exon to first base of next exon
  gtf.data |>
    dplyr::filter(type == "exon", transcript_id == transcript.id) |>
    dplyr::arrange(start, end) |>
    dplyr::mutate(j1 = end, j2 = lead(start)) |>
    dplyr::select(gene_id, transcript_id, j1, j2, strand) |>
    na.omit() # final exon of the transcript has no next junction
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
  any((junctions$j1 == junction.start & junctions$j2 == junction.end) |
    (junctions$j1 == junction.end & junctions$j2 == junction.start))
}

#' Classify a splice junction relative to a reference gene and transcript.
#'
#' @param reference.gene.id the reference gene id in the GTF
#' @param reference.transcript.id the reference transcript variant in the GTF
#' @param reference.junctions the splice junctions expected from the GTF
#' @param junction.start the start base of the junction
#' @param junction.end the end base of the junction
#' @param junction.strand the strand of the junction, one of FORWARD or REVERSE
#' @param count the number of instance of this junction observed
#'
#' @returns
#' @export
#'
#' @examples
classify_junction <- function(reference.gene.id, reference.transcript.id, reference.junctions,
                              junction.start, junction.end, junction.strand, count) {
  if (is.na(junction.start) | is.na(junction.end)) {
    return(data.frame(
      start = junction.start,
      end = junction.end,
      strand = ifelse(junction.strand == "FORWARD", "+", "-"),
      count = count,
      junction.strand = junction.strand,
      type = "novel"
    ))
  }

  # Check for matching transcripts
  full.matches <- reference.junctions |>
    dplyr::filter(j1 == junction.start, j2 == junction.end)

  if (nrow(full.matches) > 0) {
    transcript.strand <- unique(full.matches$strand)
    if (reference.transcript.id %in% full.matches$transcript_id) {
      return(data.frame(
        start = junction.start,
        end = junction.end,
        strand = transcript.strand,
        count = count,
        junction.strand = junction.strand,
        type = "reference_transcript_full_junction"
      ))
    }

    if (reference.gene.id %in% full.matches$gene_id) {
      return(data.frame(
        start = junction.start,
        end = junction.end,
        strand = transcript.strand,
        count = count,
        junction.strand = junction.strand,
        type = "reference_gene_alternative_transcript_full_junction"
      ))
    }
    return(data.frame(
      start = junction.start,
      end = junction.end,
      strand = ifelse(junction.strand == "FORWARD", "+", "-"),
      count = count,
      junction.strand = junction.strand,
      type = "non-reference_gene_full_junction"
    ))
  }


  partial.matches <- reference.junctions |>
    dplyr::filter(j1 == junction.start | j2 == junction.end)

  if (nrow(partial.matches) > 0) {
    transcript.strand <- unique(partial.matches$strand)
    if (reference.transcript.id %in% partial.matches$transcript_id) {
      return(data.frame(
        start = junction.start,
        end = junction.end,
        strand = transcript.strand,
        count = count,
        junction.strand = junction.strand,
        type = "reference_transcript_one_junction"
      ))
    }

    if (reference.transcript.id %in% partial.matches$gene_id) {
      return(data.frame(
        start = junction.start,
        end = junction.end,
        strand = transcript.strand, count = count,
        junction.strand = junction.strand,
        type = "reference_gene_alternative_transcript_one_junction"
      ))
    }
    return(data.frame(
      start = junction.start,
      end = junction.end,
      strand = ifelse(junction.strand == "FORWARD", "+", "-"),
      count = count,
      junction.strand = junction.strand,
      type = "non-reference_gene_one_junction"
    ))
  }

  return(data.frame(
    start = junction.start,
    end = junction.end, strand = ifelse(junction.strand == "FORWARD", "+", "-"),
    junction.strand = junction.strand, count = count, type = "novel"
  ))
}

#' Import GTF data from a vector of file paths and store in a named list
#'
#' The GTF data is converted to a data frame.
#'
#' @param gtf.files the file paths to read
#' @param gtf.names the names to give the GTF data
#'
#' @returns a list in which [[gtf.names[1]]] contains the data from gtf.files[1]
#' @export
#'
#' @examples
#' gtf.data <- read_gtf_data(
#'   c("genomes/GRCh38.gtf", "genomes/Sscrofa11.1.gtf"),
#'   c("human", "pig")
#' )
#' gtf.data[[human]]
read_gtf_data <- function(gtf.files, gtf.names) {
  cat("Plot sashimi: Reading full genome GTF files\n")
  gtf.data <- mclapply(gtf.files, \(x){
    df <- as.data.frame(rtracklayer::import(x)) %>%
      # Create gene_name column if missing. Ensembl vs NCBI GTF format
      # We need to use magrittr pipe for this to access input data frame as .
      dplyr::mutate(
        gene_name = if ("gene_name" %in% colnames(.)) gene_name else if ("gene" %in% colnames(.)) gene else gene_id,
        transcript_name = if ("transcript_name" %in% colnames(.)) transcript_name else if ("transcript" %in% colnames(.)) transcript else transcript_id,
      )
    cat("Plot sashimi: Read ", x, "\n")
    df
  },
  mc.cores = ifelse(installr::is.windows(), 1, 6)
  )
  names(gtf.data) <- gtf.names
  cat("Plot sashimi: Read full genome GTF files\n")
  gtf.data
}


#' Get intron and exon boundaries for all transcripts within the given region of
#' a GTF
#'
#' @param gtf.data the GTF data as read by e.g. rtracklayer
#' @param chr the chromosome
#' @param loc.start the start of the window
#' @param loc.end the end of the window
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

  # cat("Region ", chr, ":", loc.start, "-", loc.end, "contains", nrow(region), "data rows\n")

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
    dplyr::select(seqnames, gene_id, gene_name, transcript_id, transcript_name, type,
      start = intron.start, end = intron.end, strand, length
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(transcript_id, start, end)

  # cat("Region ", chr, ":", loc.start, "-", loc.end, "contains", nrow(region.introns), "intron rows\n")

  region.exons <- region |>
    as.data.frame() |>
    dplyr::filter(type == "exon") |>
    dplyr::select(
      seqnames, gene_id, gene_name, transcript_id, transcript_name, type,
      start, end, strand
    ) |>
    dplyr::mutate(length = end - start + 1) |>
    dplyr::arrange(transcript_id, start, end)

  # cat("Region ", chr, ":", loc.start, "-", loc.end, "contains", nrow(region.exons), "exon rows\n")

  region.junctions <- do.call(rbind, lapply(unique(region.exons$transcript_id), get_transcript_junctions, gtf.data = gtf.data))

  list(
    exons = region.exons,
    introns = region.introns,
    junctions = region.junctions
  )
}


#' Read splice junctions from a bam file and integrate with GTF data
#'
#' Read junction and coverage data from the bam. Build exon and intron
#' annotations from the GTF data and count known or novel splice junctions.
#'
#' @param bam.file the bam file to read
#' @param gtf.data the full genome annotation
#' @param chr the chromosome of the window to report
#' @param start the start coordinate of the window to report
#' @param end the end coordinate of the window to report
#' @param reference.gene.id the gene to use as the strand reference
#' @param reference.transcript.id the transcript to use as the strand reference.
#'   If not provided, the longest transcript of the gene will be used
#'
#' @returns a list with junction, coverage and annotation data elements
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

  cat("Reading sashimi data from bam file in region '", paste0(chr, ":", start, "-", end), "'\n")
  sashimi.data <- list()
  sashimi.data$input.file <- bam.file
  bam.data <- read_bam(bam.file, paste0(chr, ":", start, "-", end), "SENSE")

  sashimi.data$reference.gtf.region <- gtf.data[gtf.data$seqnames == chr &
    gtf.data$start >= start &
    gtf.data$end <= end, ]

  sashimi.data$reference.gene.id <- reference.gene.id
  cat("Reference gene id is '", reference.gene.id, "'\n")

  # Get the longest transcript in the gene if none specified
  if (is.na(reference.transcript.id)) {
    cat("No reference transcript given, selecting longest for gene id", sashimi.data$reference.gene.id, "\n")
    reference.transcript.id <- sashimi.data$reference.gtf.region |>
      dplyr::filter(gene_id == sashimi.data$reference.gene.id, type == "transcript") |>
      dplyr::mutate(length = end - start + 1) |>
      dplyr::arrange(length) |>
      dplyr::slice_tail(n = 1) |>
      dplyr::select(transcript_id) |>
      dplyr::pull()
  }
  sashimi.data$reference.transcript.id <- reference.transcript.id
  cat("Reference transcript is", reference.transcript.id, "\n")

  if (nrow(sashimi.data$reference.gtf.region[sashimi.data$reference.gtf.region$gene_id == sashimi.data$reference.gene.id, ]) == 0) {
    stop("Unable to detect reference gene id in genome region GTF\n")
  }

  # Is the gene on the forward or reverse strand? Note - this is the gene, not the junctions or reads
  sashimi.data$reference.transcript.strand <- sashimi.data$reference.gtf.region |>
    dplyr::filter(
      transcript_id == sashimi.data$reference.transcript.id,
      type == "exon"
    ) |>
    dplyr::select(strand) |>
    dplyr::distinct() |>
    dplyr::pull(strand) |>
    as.character()

  sashimi.data$non.reference.transcript.strand <- ifelse(sashimi.data$reference.transcript.strand == "+", "-", "+")

  cat("Reference transcript is on strand '", paste(sashimi.data$reference.transcript.strand, collapse = ","), "'\n")

  sashimi.data$reference.transcript.boundaries <- get_exon_boundaries(sashimi.data$reference.gtf.region, chr, start, end)

  sashimi.data$reference.gene.name <- sashimi.data$reference.gtf.region |>
    dplyr::filter(gene_id == reference.gene.id) |>
    dplyr::select(gene_name) |>
    dplyr::distinct() |>
    dplyr::pull()

  cat("Detected", nrow(bam.data$junctions), "splice junctions\n")
  if (nrow(bam.data$junctions) > 0) {
    sashimi.data$junctions.all <- do.call(rbind, mapply(classify_junction,
      junction.start = bam.data$junctions$start,
      junction.end = bam.data$junctions$end,
      junction.strand = bam.data$junctions$strand,
      count = bam.data$junctions$count,
      MoreArgs = list(
        reference.transcript.id = reference.transcript.id,
        reference.gene.id = reference.gene.id,
        reference.junctions = sashimi.data$reference.transcript.boundaries$junctions
      ),
      SIMPLIFY = FALSE
    ))
  } else {
    # No junctions detected, return the coverage
    sashimi.data$junctions <- data.frame(
      start = c(),
      end = c(),
      strand = c(),
      count = c(),
      junction.strand = c(),
      type = c()
    )

    # With no junctions, just use coverage. The strand with the max coverage
    # is probably the transcript stand - but we have very low coverage anyway
    sashimi.data$reference.read.strand <- ifelse(sum(bam.data$coverage$forward.strand) > sum(bam.data$coverage$reverse.strand),
      "FORWARD", "REVERSE"
    )

    sashimi.data$coverage <- bam.data$coverage |>
      dplyr::rowwise() |>
      dplyr::mutate(
        reference.strand = ifelse(sashimi.data$reference.read.strand == "FORWARD", forward.strand, reverse.strand),
        non.reference.strand = ifelse(sashimi.data$reference.read.strand == "FORWARD", reverse.strand, forward.strand),
      )

    return(sashimi.data)
  }

  # The mapping may be reversed e.g. if the wrong strand was set in mapping.
  # Check the strand the reference transcript is on, and swap forward and reverse
  # strands of coverage and junctions if needed
  sashimi.data$reference.read.strand <- sashimi.data$junctions.all |>
    dplyr::filter(type == "reference_transcript_full_junction") |>
    dplyr::group_by(junction.strand) |>
    dplyr::summarise(count = sum(count)) |>
    dplyr::arrange(count) |>
    dplyr::slice_tail(n = 1) |>
    dplyr::select(junction.strand) |>
    dplyr::pull()

  # Now we know the reference strand, we can map the forward and reverse reads here
  # to the genome + or - strand confidently
  sashimi.data$junctions.all <- sashimi.data$junctions.all |>
    dplyr::mutate(corrected.strand = ifelse(junction.strand == sashimi.data$reference.read.strand,
      sashimi.data$reference.transcript.strand,
      sashimi.data$non.reference.transcript.strand
    )) |>
    dplyr::rowwise() |>
    dplyr::mutate(corrected.strand = case_when(type == "reference_transcript_full_junction" ~ sashimi.data$reference.transcript.strand,
      type == "reference_transcript_one_junction" ~ sashimi.data$reference.transcript.strand,
      type == "reference_gene_alternative_transcript_full_junction" ~ sashimi.data$reference.transcript.strand,
      type == "reference_gene_alternative_transcript_one_junction" ~ sashimi.data$reference.transcript.strand,
      .default = corrected.strand
    ))

  # Combine junctions on the same strand and location
  sashimi.data$junctions <- sashimi.data$junctions.all |>
    dplyr::group_by(start, end, type, corrected.strand) |>
    dplyr::summarise(
      count = sum(count),
      .groups = "drop_last"
    ) |>
    dplyr::rename(strand = corrected.strand)

  # Do we need to reverse the coverage values? Check if the reference transcript reads
  # are on the expected strand.
  sashimi.data$coverage <- bam.data$coverage |>
    dplyr::rowwise() |>
    dplyr::mutate(
      reference.strand = ifelse(sashimi.data$reference.read.strand == "FORWARD", forward.strand, reverse.strand),
      non.reference.strand = ifelse(sashimi.data$reference.read.strand == "FORWARD", reverse.strand, forward.strand),
    )

  # sashimi.data$junctions <- bam.data$junctions
  # sashimi.data$junctions.reference.strand <- bam.data$reference.strand


  # junctions <- get_transcript_junctions(transcript.id, gtf.data)
  # any(junctions$j1 == junction.start & junctions$j2 == junction.end)
  # sashimi.data$junctions$matches.reference.splice.junction <- mapply(junction_is_in_GTF,
  #   junction.start = bam.data$junctions$start,
  #   junction.end = bam.data$junctions$end,
  #   MoreArgs = list(
  #     transcript.id = reference.transcript.id,
  #     gtf.data = sashimi.data$reference.gtf.region
  #   )
  # )
  # If we have junctions that match know exon boundaries, but are on the
  # wrong strand, this is usually because the library was not stranded.
  # Combine onto the correct strand.
  #
  #   junctions.grouped <- sashimi.data$junctions |>
  #     dplyr::group_by(start, end) |>
  #     dplyr::filter(matches.reference.splice.junction) |> # reference transcript only
  #     dplyr::summarise(
  #       strand = "SENSE",
  #       count = sum(count),
  #       matches.reference.splice.junction = TRUE,
  #       .groups = "drop_last"
  #     )
  #
  #   sashimi.data$junctions <- sashimi.data$junctions |>
  #     dplyr::filter(!matches.reference.splice.junction) |>
  #     rbind(junctions.grouped)


  return(sashimi.data)
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

  exon.data <- exon.data[exon.data$strand == strand, ]
  intron.data <- intron.data[intron.data$strand == strand, ]
  cat("Creating intron collapser; detected", nrow(exon.data), "exons and", nrow(intron.data), "introns\n")

  # Reduce any overlapping exons if we have multiple transcripts
  exon.ranges <- GenomicRanges::reduce(GenomicRanges::GRanges(
    seqnames = rep("test", nrow(exon.data)),
    ranges = IRanges::IRanges(
      start = exon.data$start,
      end = exon.data$end
    ),
    strand = strand
  ))

  # Break overlapping introns apart, since they can be part of an exon for a
  # different transcript
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

  # Some introns or intergenic sequence may be missing. Fill in gaps from min start to max end that are
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

  # Create the new start and end coordinates with the desired scaling
  full.ranges <- full.ranges |>
    dplyr::arrange(start, end) |>
    dplyr::distinct() |>
    dplyr::mutate(
      original.length = end - start + 1,
      new.length = ifelse(original.length > max.intron.length & Type == "intron", max.intron.length, original.length),
      new.end = min(start) + cumsum(new.length) - 1,
      new.start = new.end - new.length + 1, # how much offset to apply
      is.ordered = new.start < new.end,
      is.contiguous = new.start == dplyr::lag(new.end) + 1,
      is.full.coverage = start == dplyr::lag(end) + 1,
      step.size = new.length / original.length
    )

  range.min <- min(full.ranges$start)
  range.max <- max(full.ranges$end)
  new.range.max <- max(full.ranges$new.end)

  # Create a mapping of old to new coordinate
  lookup.table <- do.call(rbind, mapply(\(new.start, new.end, old.start, old.end, original.length, new.length){
    data.frame(
      old.position = old.start:old.end,
      new.position = ((as.double(0:(old.end - old.start))) / original.length * new.length) + new.start
    )
  }, full.ranges$new.start, full.ranges$new.end, full.ranges$start, full.ranges$end, full.ranges$original.length, full.ranges$new.length, SIMPLIFY = FALSE))

  # rownames(lookup.table) <- lookup.table$old.position

  # Create a function that uses the above tables to convert a coordinate vector
  # to the new ranges.
  calculate <- function(coordinate) {
    pre <- coordinate[coordinate < range.min]
    post <- coordinate[coordinate > range.max]

    middle <- coordinate[coordinate >= range.min & coordinate <= range.max]
    indexes <- sapply(middle, \(x)which(lookup.table$old.position == x))

    middle <- lookup.table[indexes, "new.position"]

    post <- (post - range.max) + new.range.max

    c(pre, middle, post)
  }

  list(
    calculate = calculate,
    full.ranges = full.ranges,
    lookup.table = lookup.table
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
  # Calculate offsets
  intron.collapser <- create_intron_collapser(exon.data, intron.data, sashimi.data$reference.transcript.strand)

  # Exons
  sashimi.data$reference.transcript.boundaries$exons$old.start <- sashimi.data$reference.transcript.boundaries$exons$start
  sashimi.data$reference.transcript.boundaries$exons$old.end <- sashimi.data$reference.transcript.boundaries$exons$end
  sashimi.data$reference.transcript.boundaries$exons$start <- intron.collapser$calculate(sashimi.data$reference.transcript.boundaries$exons$start)
  sashimi.data$reference.transcript.boundaries$exons$end <- intron.collapser$calculate(sashimi.data$reference.transcript.boundaries$exons$old.end)

  # Introns
  sashimi.data$reference.transcript.boundaries$introns$old.start <- sashimi.data$reference.transcript.boundaries$introns$start
  sashimi.data$reference.transcript.boundaries$introns$old.end <- sashimi.data$reference.transcript.boundaries$introns$end
  sashimi.data$reference.transcript.boundaries$introns$start <- intron.collapser$calculate(sashimi.data$reference.transcript.boundaries$introns$start)
  sashimi.data$reference.transcript.boundaries$introns$end <- intron.collapser$calculate(sashimi.data$reference.transcript.boundaries$introns$old.end)

  # Splice junctions
  if (nrow(sashimi.data$junctions) > 0) {
    sashimi.data$junctions$old.start <- sashimi.data$junctions$start
    sashimi.data$junctions$old.end <- sashimi.data$junctions$end
    sashimi.data$junctions$start <- intron.collapser$calculate(sashimi.data$junctions$old.start)
    sashimi.data$junctions$end <- intron.collapser$calculate(sashimi.data$junctions$old.end)
  }

  # Coverage
  sashimi.data$coverage$old.position <- sashimi.data$coverage$position
  sashimi.data$coverage$position <- intron.collapser$calculate(sashimi.data$coverage$old.position)

  sashimi.data
}

#### Charting functions ####
#' Plot splice junctions on a gene exon track
#'
#' @param sashimi.data sashimi data as produced by read_sashimi_data
#' @param min.spanning.reads  the minimum number of spanning reads for a splice junction to be plotted
#' @param label the label for the gene track
#' @param show.x.axis if true, display the x axis with genome coordinates
#' @param is.collapse_introns if true, make each intron at most 500bp wide. If this option is used, the x-axis will no longer match the genome coordinates
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
  xmin <- min(sashimi.data$coverage$position)
  xmax <- max(sashimi.data$coverage$position)

  max.coverage <- max(sashimi.data$coverage$coverage)

  intron.y <- 1.5
  exon.ymin <- 1.1
  exon.ymax <- 2

  # Make the gene track
  splot <- ggplot() +
    # Add read coverage, scaled to 0-1
    geom_area(data = sashimi.data$coverage, aes(x = position, y = reference.strand / max.coverage), fill = "darkgrey", col = "darkgrey") +
    geom_area(data = sashimi.data$coverage, aes(x = position, y = -non.reference.strand / max.coverage), fill = "darkgrey", col = "darkgrey") +

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
      label = sashimi.data$non.reference.transcript.strand,
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
  add.junction <- function(splot, xmin, xmax, is.lower, ymin, ymax, count, junction.type) {
    # Define the spline line styles that make the junction lines
    SPLINE.COLOURS <- c(
      "reference_transcript_full_junction" = "blue",
      "reference_gene_alternative_transcript_full_junction" = "black",
      "reference_transcript_one_junction" = "black",
      "reference_gene_alternative_transcript_one_junction" = "black",
      "non-reference_gene_full_junction" = "grey",
      "non-reference_gene_one_junction" = "grey",
      "novel" = "grey"
    )

    spline.color <- SPLINE.COLOURS[junction.type]
    spline.size <- 1
    spline.alpha <- 1

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
    if (is.lower) { # Junctions below zero
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
    if (is.lower) { # Junctions below zero
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
      y = ifelse(is.lower, 0 - ymax - label.y.offset, ymax + label.y.offset),
      label = as.character(count),
      size = 2, col = spline.color, fill = NA
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
        is.lower = jrow$strand != sashimi.data$reference.transcript.strand,
        ymin = 2, ymax = 2.5 + jrow$y.offset,
        count = jrow$count,
        junction.type = jrow$type
      )
    }
  }

  list(plot = splot, sashimi.data = sashimi.data)
}
# TODO - chicken has an issue - the antisense transcript is on the wrong strand.
# Check for other anomalies before finalising the junction merging.
