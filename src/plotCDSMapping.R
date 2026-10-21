# Read a complete bam file from CDS mapping and plot the reads. We want
# to be able to detect missing e2 if present

source("src/functions.R")

fs::dir_create("cds/report")

cds.samples <- read.csv(CDS.MAPPING.FILE) |>
  dplyr::mutate(path = paste0("cds/data/", CommonName, "/", Run, ".bam"))|>
  dplyr::group_by(CommonName, Genome, Sex, Tissue, Timepoint) |>
  dplyr::summarise(bams = list(path))|>
  dplyr::mutate(name = paste(CommonName, Sex, Tissue, Timepoint, sep="_"))

# Given an output file name, and bam files to aggregate, plot the presence
# of spliced reads in the file
plot.bam.reads <- function(name, bams){
  
  flog.info("Making %s", name)
  
  if(!all(file.exists(bams))) {
    flog.warn("Missing %s", paste(bams[!file.exists(bams)], collapse = ", "))
    return()
  }
 
  bam.data <- do.call(rbind, lapply(bams, \(x){ as.data.frame(Rsamtools::scanBam(x)[[1]])}))
  
  # Split the cigar strings to find read positions on the reference
  bam.table <- bam.data |>
    tidyr::separate_wider_regex(cigar, patterns = c(match1="^\\d+[A-Z]", match2="\\d+[A-Z]", 
                                                    match3="\\d+[A-Z]", match4="\\d+[A-Z]",
                                                    match5="\\d+[A-Z]"), too_few = "debug") |>
    tidyr::pivot_longer(match1:match5, names_to = "Match", values_to = "cigar_match") |>
    dplyr::filter(!is.na(cigar_match)) |>
    dplyr::mutate(
      cigar_type = stringr::str_extract(cigar_match, "[A-Z]"),
      cigar_number = as.numeric(stringr::str_extract(cigar_match, "\\d+"))) |>
    dplyr::group_by(qname, strand) |>
    dplyr::arrange(qname, strand, Match) |>
    dplyr::mutate(
      has_junction = any(cigar_type=="N"),
      cigar_start = pos) |>
    dplyr::mutate(
      cigar_start = case_when(Match == "match1" ~ cigar_start,
                              Match=="match2"  ~ cigar_start + lag(cigar_number),
                              Match=="match3"  ~ cigar_start + lag(cigar_number) +lag(cigar_number, 2),
                              Match=="match3"  ~ cigar_start + lag(cigar_number) +lag(cigar_number, 2)+lag(cigar_number, 3),
                              Match=="match4"  ~ cigar_start + lag(cigar_number) +lag(cigar_number, 2)+lag(cigar_number, 3),
                              .default = NA),
      cigar_stop = cigar_start + cigar_number-1,
      
    ) |>
    dplyr::filter(cigar_type!="S")
  
  
  # Plot the reads. We want to highlight reads with a potential splice junction. 
  # If any samples look interesting, we can then match up the exon coordinates.
  bam.plot <- ggplot(bam.table)+
    geom_segment(aes(x = cigar_start, xend = cigar_stop, 
                     # y = yPos, yend = yPos,
                     y = reorder(qname, cigar_start), yend=reorder(qname, cigar_start),
                     col = has_junction,
                     linetype=cigar_type))+
    labs(x = "Position", y = "Read", title = name)+
    scale_color_manual(values = c(`TRUE`="red", `FALSE`="lightgrey"))+
    scale_linetype_manual(values = c("M"=1, "N"=3, "S"=1, "I"=1, "D"=1))+
    facet_wrap(~rname, scales="free_y")+
    theme_bw()+
    theme(axis.text.y = element_blank(),
          panel.grid = element_blank())
  
  save.double.width(paste0("cds/report/",name),  bam.plot)
  
}

invisible(mcmapply(plot.bam.reads, cds.samples$name, cds.samples$bams, 
                   mc.cores = DEFAULT.MC.CORES))


# # Read all the bam files
# bam.files <- data.frame(path = list.files(path = "cds/data", pattern = "[SDE]RR\\d+.*\\.bam$", 
#                                           full.names = TRUE, recursive=TRUE)) |>
#   dplyr::mutate(file = basename(path),
#                 mtime = file.mtime(path)) |>
#   tidyr::separate_wider_delim(file,
#                               delim = ".", names = c("Run", "ext")
#   )
# 
# bam.data <- do.call(rbind, lapply(bam.files$path, \(x){ as.data.frame(Rsamtools::scanBam(x)[[1]])}))
# 
# # Split the cigar strings to find read positions on the reference
# bam.table <- bam.data |>
#   tidyr::separate_wider_regex(cigar, patterns = c(match1="^\\d+[A-Z]", match2="\\d+[A-Z]", 
#                                                   match3="\\d+[A-Z]", match4="\\d+[A-Z]",
#                                                   match5="\\d+[A-Z]"), too_few = "debug") |>
#   tidyr::pivot_longer(match1:match5, names_to = "Match", values_to = "cigar_match") |>
#   dplyr::filter(!is.na(cigar_match)) |>
#   dplyr::mutate(
#     cigar_type = stringr::str_extract(cigar_match, "[A-Z]"),
#     cigar_number = as.numeric(stringr::str_extract(cigar_match, "\\d+"))) |>
#   dplyr::group_by(qname, strand) |>
#   dplyr::arrange(qname, strand, Match) |>
#   dplyr::mutate(
#     has_junction = any(cigar_type=="N"),
#     cigar_start = pos) |>
#   dplyr::mutate(
#     cigar_start = case_when(Match == "match1" ~ cigar_start,
#                             Match=="match2"  ~ cigar_start + lag(cigar_number),
#                             Match=="match3"  ~ cigar_start + lag(cigar_number) +lag(cigar_number, 2),
#                             Match=="match3"  ~ cigar_start + lag(cigar_number) +lag(cigar_number, 2)+lag(cigar_number, 3),
#                             Match=="match4"  ~ cigar_start + lag(cigar_number) +lag(cigar_number, 2)+lag(cigar_number, 3),
#                             .default = NA),
#     cigar_stop = cigar_start + cigar_number-1,
#     
#  ) |>
#   dplyr::filter(cigar_type!="S")
# 
# 
# # Plot the reads. We want to highlight reads with a potential splice junction. 
# # If any samples look interesting, we can then match up the exon coordinates.
# ggplot(bam.table)+
#   geom_segment(aes(x = cigar_start, xend = cigar_stop, 
#                    # y = yPos, yend = yPos,
#                    y = reorder(qname, cigar_start), yend=reorder(qname, cigar_start),
#                    col = has_junction,
#                    linetype=cigar_type))+
#   labs(x = "Position", y = "Read")+
#   scale_color_manual(values = c(`TRUE`="red", `FALSE`="lightgrey"))+
#   scale_linetype_manual(values = c("M"=1, "N"=3, "S"=1, "I"=1, "D"=1))+
#   facet_wrap(~rname, scales="free_y")+
#   theme_bw()+
#   theme(axis.text.y = element_blank(),
#         panel.grid = element_blank())

