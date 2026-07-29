# Find paired end reads in bed files with mate1, mate2 coordinates

#### Read bed files ####

bed.files <- list.files(path = "data", pattern = "[SDE]RR.*.bed$", 
                                          full.names = TRUE, recursive=TRUE)

# Read the bed file. Keep only reads with pairs on opposite strands and on the
# same chromosome.
read.bed <- function(file){
  if(file.size(file)==0) return()
  read_tsv(file, col_names = c("Seqname", "Mate1Start", "Mate1End", "Seqname2", 
                               "Mate2Start", "Mate2End", "ReadName", "Score", 
                               "Mate1Strand", "Mate2Strand"),
           col_types = c("ciiciicicc")) |>
    dplyr::filter(Mate1Strand!=Mate2Strand, 
                  Seqname==Seqname2,
                  Mate1Start!=Mate2Start,
                  Mate1End != Mate2End) |>
    dplyr::mutate(File = file,
                  Run = str_extract(File, "([SDE]RR\\d+)", group = 1),
                  GeneId = str_extract(File, "[SDE]RR\\d+\\.([\\w\\d]+)\\.bam.bed$", group = 1)
                  ) |>
    dplyr::select(-Seqname2, -File, -Score)
}


bed.data <- do.call(dplyr::bind_rows, lapply(bed.files, read.bed) ) |>
  dplyr::rowwise() |>
  dplyr::mutate(  
    Mate1Size = Mate1End - Mate1Start + 1,
    Mate2Size = Mate2End - Mate2Start + 1,
    InsertStart = min(Mate1Start, Mate2Start),
    InsertEnd = max(Mate1End, Mate2End),
    InsertSize = InsertEnd - InsertStart + 1,
    JunctionStart = min(Mate1End, Mate2End),
    JunctionEnd = max(Mate1Start, Mate2Start)
    ) |>
  merge(SELECTED.SAMPLES, by = c("Run")) |>
  merge(GENE.LOCATIONS, by = c("CommonName", "GeneId", "GTF_FILE"))


#### Match reads to exon 2 boundaries ####

# junction coordinates are set in plotSashimi
splice.data <- bed.data |>
  dplyr::filter(Group %in% c("ZFX", "ZFY"),
                ) |>
  merge(junction.coordinates, by = c("CommonName", "GeneId")) |>
  dplyr::mutate(SpansExon2 = JunctionStart <= start & JunctionEnd >= end) |>
  dplyr::group_by(Run, CommonName, GeneId, Gene, sex, Organism_part, Timepoint) |>
  dplyr::mutate(RunMedianInsertSize = median(InsertSize))

splice.summ <- splice.data |> 
  dplyr::group_by(CommonName, GeneId, Gene, sex, Organism_part, Timepoint, SpansExon2) |>
  dplyr::summarise(GroupInsertSize = median(InsertSize))|>
  tidyr::pivot_wider(names_from = SpansExon2, names_prefix = "SpansExon2", values_from = c(GroupInsertSize)) |>
  dplyr::filter(!is.na(SpansExon2TRUE)) |>
  as.data.frame()

plot.spanning <- function(i){

# for(i in 1:nrow(splice.summ)){
  species <- splice.summ[i, "CommonName"]
  tissue <- splice.summ[i, "Organism_part"]
  timepoint <- splice.summ[i, "Timepoint"]
  gene <- splice.summ[i, "Gene"]
  s <- splice.summ[i, "sex"]
  
  subset.data <- splice.data |>
    dplyr::filter(CommonName==species, Timepoint==timepoint, Organism_part==tissue, Gene==gene,
                  sex==s,
                  SpansExon2) |>
    dplyr::arrange(Mate1Start, Mate2Start)

  if(nrow(subset.data)==0) next
  if(!any(subset.data$SpansExon2)) next

  subset.data$id <- 1:nrow(subset.data)
  
  location <- parse_coordinates(unique(subset.data$Location))
  
  anns <- get_exon_boundaries(GTF.DATA[[species]], location$coord.chr, location$coord.start, location$coord.end)
  
  # Plot exons from the reference transcript
  reference.exons <- anns$exons |>
    dplyr::filter(
      transcript_id == unique(subset.data$CanonicalTranscriptId)
    )
  
  reference.introns <- anns$introns |>
    dplyr::filter(
      start > min(reference.exons$start),
      end < max(reference.exons$start),
      transcript_id == unique(subset.data$CanonicalTranscriptId)
    )
  
  spanning.plot <- ggplot(subset.data) +
    geom_rect(
      data = reference.exons, aes(
        xmin = start,
        xmax = end,
        ymin = -Inf,
        ymax = Inf
      ),
      fill = "grey", alpha = 1
    ) +
    geom_segment(aes(x = Mate1Start, xend = Mate1End, y = id, yend = id), linewidth = 2, col="blue", , alpha=0.5)+
    geom_segment(aes(x = Mate2Start, xend = Mate2End, y = id, yend = id), linewidth = 2, col="orange", alpha=0.5)+
    geom_segment(aes(x = Mate1End, xend = Mate2Start, y = id, yend = id), col="black")+
    labs(x = "Position", y = "Read pair", title = paste(species, s, tissue, timepoint, gene))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank())
  
  save.double.width(paste0("report/spanning_pairs/", species, ".",s, ".", tissue,
                           ".", timepoint, ".", gene, ".png"),spanning.plot, width = 170, height = min(1000, 50 + nrow(subset.data)*0.5))
  spanning.plot
}

spanning.plots <- lapply(1:nrow(splice.summ), plot.spanning)

save.double.width("report/spanning_pairs/spanning_pairs.png", patchwork::wrap_plots(spanning.plots, ncol = 1), height = 250)
