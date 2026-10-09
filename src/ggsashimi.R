# R implementation of ggsashimi.py from
# https://github.com/guigolab/ggsashimi

packages <- c(
  "parallel", "tidyverse", "GenomicRanges", "grid",
  "fs", "data.table", "rtracklayer", "Rsamtools", "bitops", "rlang",
  "data.table", "installr"
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

# Number of parallel cores to use in mclapply - only works on Unix
DEFAULT.MC.CORES <- ifelse(installr::is.windows(), 1, 6)

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
  # Fastest hash access in R is via environment name lookup.
  # Make an environment for read counts

  coverage.array <- new.env()
  coverage.array$forward <- new.env()
  coverage.array$reverse <- new.env()

  n.indexes <- end - start + 1

  for (i in as.character(start:end)) {
    env_poke(coverage.array$forward, i, 0)
    env_poke(coverage.array$reverse, i, 0)
  }

  coverage.array$add <- function(pos, strand, value) {
    pos <- as.character(pos)
    if (strand == "FORWARD") {
      curr <- env_get(coverage.array$forward, pos)
      env_poke(coverage.array$forward, pos, curr + value)
    } else {
      curr <- env_get(coverage.array$reverse, pos)
      env_poke(coverage.array$reverse, pos, curr + value)
    }
  }

  coverage.array$increment <- function(pos, strand) {
    coverage.array$add(pos, strand, 1)
  }

  coverage.array$incrementRange <- function(start, end, strand) {
    base.range <- as.character(start:end)

    for (i in base.range) {
      coverage.array$increment(i, strand)
    }
  }

  coverage.array$get <- function() {
    fwd <- env_get_list(coverage.array$forward, nms = as.character(start:end))
    rev <- env_get_list(coverage.array$reverse, nms = as.character(start:end))

    merge(data.frame(position = names(fwd), forward.strand = unlist(fwd)),
      data.frame(position = names(rev), reverse.strand = unlist(rev)),
      by = "position", all = TRUE
    ) |>
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
#' parse_coordinates("NC1234.5:-1234-5678")
parse_coordinates <- function(coordinate.string) {
  coordinate.string <- stringr::str_remove_all(coordinate.string, ",")
  # Ensure we can also handle a coordinate that is negative e.g. treeshrew ZFX
  parsed <- stringr::str_extract(coordinate.string, "^([\\._A-Za-z\\d]+):(-?\\d+)-(-?\\d+)$", group = 1:3)

  list(
    coord.chr = parsed[1],
    coord.start = as.integer(parsed[2]),
    coord.end = as.integer(parsed[3])
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

  # scanBam will automatically use a bai index if available
  bam.data <- Rsamtools::scanBam(bam.file)[[1]]

  bam.data <- data.frame(
    flag = bam.data$flag,
    cigar = bam.data$cigar,
    pos = bam.data$pos
  ) |>
    dplyr::mutate(
      # Determine the strand of each read in the genome
      read.strand = sapply(flag, find_read_strand, strand.string = strand.string),
      # Skip if read or mate is unmapped
      is.unmapped = sapply(flag, has_sam_flag, property = SAM.FLAG.READ.UNMAPPED),
      is.mate.unmapped = sapply(flag, has_sam_flag, property = SAM.FLAG.MATE.UNMAPPED),
      # Ignore reads with more exotic CIGAR operators
      is.nonstandard.read = sapply(cigar, \(x) any(stringr::str_detect(x, c("H", "P", "X", "="))))
    ) |>
    dplyr::filter(!is.unmapped, !is.mate.unmapped, !is.nonstandard.read)

  total.reads <- nrow(bam.data)

  if (total.reads == 0) {
    cat("There are no reads in bam file", bam.file, "\n")
    return(list(
      total.reads = 0,
      coverage = data.frame(),
      junctions = data.frame(),
      reference.strand = strand.string,
      reference.location = coordinate.string
    ))
  } else {
    cat("There are", total.reads, "valid reads in the bam file\n")
  }

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
    # Insertion or Soft-clip
    # we do not change position
    if (op == "I" | op == "S") {
      return(pos)
    }

    # Match - increase coverage across the range
    if (op == "M") {
      pos.range <- pos:(pos + len - 1)
      pos.range <- pos.range[pos.range >= coordinates$coord.start & pos.range < coordinates$coord.end]
      if (length(pos.range) > 0) {
        coverage.array$incrementRange(pos.range[1], pos.range[length(pos.range)], strand)
      }
    }

    # Deletion
    # if (op == "D") {
    #   # no action
    # }

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


  # Go read by read
  for (i in 1:total.reads) {
    if (i %% 500 == 0) cat(sprintf("Processed %.2f%% of %s reads\n", i / total.reads * 100, total.reads))
    flag <- bam.data$flag[i]
    cigar <- bam.data$cigar[i]
    pos <- bam.data$pos[i]
    read.strand <- bam.data$read.strand[i]

    # Parse the cigar string
    # e.g. 10M3I21D would parse into 10, 3, 21 and M, I, D
    cigar.lengths <- as.integer(stringr::str_extract_all(cigar, "\\d+")[[1]])
    cigar.ops <- stringr::str_extract_all(cigar, "[MIDNS]")[[1]]

    current.position <- pos

    # Check all cigar ops for this read.
    for (i in 1:length(cigar.ops)) {
      curr.cigar.length <- cigar.lengths[i]
      curr.cigar.op <- cigar.ops[i]
      current.position <- count.operator(
        curr.cigar.op, curr.cigar.length, current.position, read.strand
      )
    }
  }

  coverage <- coverage.array$get() |>
    dplyr::mutate(position = as.integer(position))

  return(list(
    total.reads = total.reads,
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
#' The GTF data is converted to a data frame. Allows GTF data to be loaded for
#' several species and reused in scripts.
#'
#' @param gtf.files a vector of GTF file paths to read
#' @param gtf.names the names to give the GTF data. Should have the same number
#'   of elements as gtf.files
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
    cat("Plot sashimi: Read", x, "\n")
    df
  },
  mc.cores = DEFAULT.MC.CORES
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
  region <- gtf.data |>
    dplyr::filter(
      seqnames == as.character(chr),
      start >= loc.start - 5000 &
        end <= loc.end + 5000
    )

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

  region.exons <- region |>
    as.data.frame() |>
    dplyr::filter(type == "exon") |>
    dplyr::select(
      seqnames, gene_id, gene_name, transcript_id, transcript_name, type,
      start, end, strand
    ) |>
    dplyr::mutate(length = end - start + 1) |>
    dplyr::arrange(transcript_id, start, end)

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
#' @param gtf.data the full genome annotation. See read_gtf_data()
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

  sashimi.data <- list()
  sashimi.data$input.file <- bam.file
  bam.data <- read_bam(bam.file, paste0(chr, ":", start, "-", end), "SENSE")
  sashimi.data$total.reads <- bam.data$total.reads
  if (bam.data$total.reads == 0) {
    return(sashimi.data)
  }

  sashimi.data$reference.gtf.region <- gtf.data[gtf.data$seqnames == chr &
    gtf.data$start >= start &
    gtf.data$end <= end, ]

  sashimi.data$reference.gene.id <- reference.gene.id

  # Get the longest transcript in the gene if none specified
  if (is.na(reference.transcript.id)) {
    cat("No reference transcript given, selecting longest for gene id '", sashimi.data$reference.gene.id, "'\n")
    reference.transcript.id <- sashimi.data$reference.gtf.region |>
      dplyr::filter(gene_id == sashimi.data$reference.gene.id, type == "transcript") |>
      dplyr::mutate(length = end - start + 1) |>
      dplyr::arrange(length) |>
      dplyr::slice_tail(n = 1) |>
      dplyr::select(transcript_id) |>
      dplyr::pull()
  }
  sashimi.data$reference.transcript.id <- reference.transcript.id

  if (nrow(sashimi.data$reference.gtf.region[sashimi.data$reference.gtf.region$transcript_id == sashimi.data$reference.transcript.id, ]) == 0) {
    stop("Cannot detect a reference transcript with id '", reference.transcript.id, "' in genome region GTF\n")
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

  cat("Reference gene is '",reference.gene.id, "', transcript is '", reference.transcript.id, "' on strand '", paste(sashimi.data$reference.transcript.strand, collapse = ","), "'\n")
  if (is.null(sashimi.data$reference.transcript.strand)) stop("Unable to find reference strand for", sashimi.data$reference.transcript.id)

  sashimi.data$reference.transcript.boundaries <- get_exon_boundaries(sashimi.data$reference.gtf.region, chr, start, end)

  sashimi.data$reference.gene.name <- sashimi.data$reference.gtf.region |>
    dplyr::filter(gene_id == reference.gene.id) |>
    dplyr::select(gene_name) |>
    dplyr::distinct() |>
    dplyr::pull()


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
    cat("No splice junctions detected\n")
    sashimi.data$junctions <- data.frame(
      start = c(),
      end = c(),
      strand = c(),
      count = c(),
      junction.strand = c(),
      type = c(),
      Pct_junction_reads = c(),
      Pct_total_reads = c()
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
    dplyr::filter(type == "reference_transcript_full_junction" | type == "reference_gene_alternative_transcript_full_junction") |>
    dplyr::group_by(junction.strand) |>
    dplyr::summarise(count = sum(count)) |>
    dplyr::arrange(count) |>
    dplyr::slice_tail(n = 1) |>
    dplyr::select(junction.strand) |>
    dplyr::pull()

  # What if there were no reference transcript junctions? Not going to use the
  # plot, so just default to forward strand.
  if (is.empty(sashimi.data$reference.read.strand)) {
    print(sashimi.data$junctions.all)
    sashimi.data$reference.read.strand <- "+"
  }

  # Now we know the reference strand, we can map the forward and reverse reads here
  # to the genome + or - strand confidently.
  sashimi.data$junctions.all <- sashimi.data$junctions.all |>
    dplyr::ungroup() |>
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
    dplyr::rename(strand = corrected.strand) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      Pct_junction_reads = count / sum(count) * 100,
      Pct_total_reads = count / bam.data$total.reads * 100
    )

  # Do we need to reverse the coverage values? Check if the reference transcript reads
  # are on the expected strand.
  sashimi.data$coverage <- bam.data$coverage |>
    dplyr::ungroup() |>
    dplyr::rowwise() |>
    dplyr::mutate(
      reference.strand = ifelse(sashimi.data$reference.read.strand == "FORWARD", forward.strand, reverse.strand),
      non.reference.strand = ifelse(sashimi.data$reference.read.strand == "FORWARD", reverse.strand, forward.strand),
    )

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

  tryCatch({
    exon.ranges <- GenomicRanges::GRanges(
      seqnames = rep("test", nrow(exon.data)),
      ranges = IRanges::IRanges(
        start = exon.data$start,
        end = exon.data$end
      ),
      strand = strand
    )

    # Reduce any overlapping exons if we have multiple transcripts
    # cat("Reducing exons\n")
    exon.ranges <- GenomicRanges::reduce(exon.ranges)
    
    # cat("Creating introns\n")
    intron.ranges <- GenomicRanges::GRanges(
      seqnames = rep("test", nrow(intron.data)),
      ranges = IRanges::IRanges(
        start = intron.data$start,
        end = intron.data$end
      ),
      strand = strand
    )
    

    
    # Remove introns that overlap an exon due to multiple transcripts. Standard
    # GenomicRanges::disjoin fails with <simpleError in
    # .Call2("C_find_overlaps_NCList", start(query), end(query), start(subject),
    # end(subject), nclist, nclist_is_q, maxgap, minoverlap, type,     select,
    # circle.length, PACKAGE = "IRanges"): build_NCList: memory allocation failed>
    # on our HPC. Replace GenomicRanges::countOverlaps with custom function here
    # only. Standard code: intron.ranges <- GenomicRanges::disjoin(intron.ranges).
    # Added a check to only run this if there are overlapping introns
    has.overlaps <- function(intron.ranges, exon.ranges){
      
      sapply(1:length(intron.ranges), \(i){
        i.start <- start(intron.ranges[i])
        i.end <- end(intron.ranges[i])
        for(j in 1:length(exon.ranges)){
          j.start <- start(exon.ranges[j])
          j.end <- end(exon.ranges[j])
          
          if( (i.start <= j.start & i.end >= j.start) |
              (i.start <= j.end & i.end >= j.end)  |
              (i.start >= j.start & i.end <=j.end )  
          ){
            return(1) # don't need an absolute count, just a numeric non-zero
          }
        }
        return(0)
      })
    }
    
    
    # Break overlapping introns apart, since they can be part of an exon for a
    # different transcript
    intron.ranges$overlappingIntrons <- has.overlaps(intron.ranges, intron.ranges)
    
    if(any(intron.ranges$overlappingIntrons>1)){
      cat("Disjoining introns\n")
      intron.ranges <- GenomicRanges::disjoin(intron.ranges)
    }

    # cat("Checking introns for overlaps with exons\n")
    
    # intron.ranges$overlappingExons <- GenomicRanges::countOverlaps(intron.ranges, exon.ranges, minoverlap = 1)
    intron.ranges$overlappingExons <- has.overlaps(intron.ranges, exon.ranges)
    intron.ranges <- intron.ranges[intron.ranges$overlappingExons == 0, ]
    
    # cat("Filling intron gaps\n")
    # Some introns or intergenic sequence may be missing. Fill in gaps from min start to max end that are
    # not covered by intron or exons
    missing.introns <- GenomicRanges::gaps(GenomicRanges::reduce(c(intron.ranges, exon.ranges)),
                                           start = min(exon.data$start)
    )
    
    # cat("Found intron and exon set\n")
    
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
    
    # cat("Created lookup table\n")
  }, error = \(e){
    cat("Error creating intron collapser, not collapsing\n")
    print(e)
    return(list(
      calculate = \(x) x,
      full.ranges = NA,
      lookup.table = NA
    ))
  })
  
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
                                       show.x.axis = TRUE, is.collapse.introns = FALSE,
                                       is.show.junction.percent = FALSE) {
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

  max.coverage <- max(sashimi.data$coverage$forward.strand, sashimi.data$coverage$reverse.strand)

  intron.y <- 1.5
  exon.ymin <- 1.1
  exon.ymax <- 2

  # Make the gene track
  splot <- ggplot() +
    # Add read coverage, scaled to 0-1
    geom_area(data = sashimi.data$coverage, aes(x = position, y = reference.strand / max.coverage), fill = "darkgrey", col = "darkgrey") +
    geom_area(data = sashimi.data$coverage, aes(x = position, y = -non.reference.strand / max.coverage), fill = "darkgrey", col = "darkgrey") +

    # Add label for the max coverage value on the coverage chart at max x value
    annotate("text",
             x = ifelse(is.x.reverse, -Inf, Inf),
             y = 1,
             label = max.coverage,
             hjust = 1.2, vjust = 0.5,
             size = 2
    ) +
    annotate("text",
             x = ifelse(is.x.reverse, -Inf, Inf),
             y = -1,
             label = max.coverage,
             hjust = 1.1, vjust = 0.5,
             size = 2
    ) +
    # Axis for coverage region
   annotate("segment",
            x = ifelse(is.x.reverse, -Inf, Inf),
            xend = ifelse(is.x.reverse, -Inf, Inf),
            y=1, yend=-1,
            linewidth=0.3, col="black"
            ) +
    
    
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

    # Add centre line
    geom_hline(yintercept = 0, linewidth=0.3, col="black") +
    
    # Add upper strand label
    annotate("text",
      x = ifelse(is.x.reverse, -Inf, Inf),
      y = intron.y + 1,
      label = sashimi.data$reference.transcript.strand,
      hjust = 1, vjust = 0.5
    ) +

    # Add lower strand label
    annotate("text",
      x = ifelse(is.x.reverse, -Inf, Inf),
      y = -intron.y - 1,
      label = sashimi.data$non.reference.transcript.strand,
      hjust = 1, vjust = 0.5
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
  add.junction <- function(splot, xmin, xmax, is.lower, ymin, ymax, count, pct.junction.reads, junction.type) {
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

    label.y.offset <- 0.0 # separation between spline and label

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

    # We may want to show the percentage of junction spanning reads
    # If so, the label offset will change.
    junction.label <- ifelse(is.show.junction.percent,
      sprintf("%s\n%.1f%%", as.character(count), pct.junction.reads),
      sprintf("%s", as.character(count))
    )

    junction.label.y.offset <- ifelse(is.show.junction.percent, 0.0, 0.4) # separation between spline and label

    # Junction count label
    splot <- splot + annotate("label",
      x = xmid,
      y = ifelse(is.lower, 0 - ymax - junction.label.y.offset, ymax + junction.label.y.offset),
      label = junction.label,
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
        pct.junction.reads = jrow$Pct_junction_reads,
        junction.type = jrow$type
      )
    }
  }

  list(plot = splot, sashimi.data = sashimi.data)
}

