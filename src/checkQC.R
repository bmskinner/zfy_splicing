# Check sample QC

library(tidyverse)
library(ggbeeswarm)
library(patchwork)
library(fs)
source("src/functions.R")

fs::dir_create("report/_qc")

cat("Reading samples\n")
filtered.samples <- read.filtered.samples()

#### Trimming report ####

cat("Checking trimming\n")

trimming.summary.files <- list.files(path = "data", pattern = "fastq.gz_trimming_report.txt$", full.names = T, recursive = T)

extract.trimming.info <- function(f){
  trimming.report <- read_file(f)
  reads.processed <- str_extract(trimming.report, "Total reads processed:\\s+[\\d|,]+")
  reads.processed <- str_replace(reads.processed, "Total reads processed:\\s+", "")
  
  reads.with.adapters <- str_extract(trimming.report, "Reads with adapters:\\s+[\\d|,|(|)|%| |\\.]+")
  reads.with.adapters <- str_replace(reads.with.adapters, "Reads with adapters:\\s+", "")
  
  reads.passing.filters <- str_extract(trimming.report, "Reads written \\(passing filters\\):\\s+[\\d|,|(|)|%| |\\.]+")
  reads.passing.filters <- str_replace(reads.passing.filters, "Reads written \\(passing filters\\):\\s+", "")
  
  file.name <- basename(f)
  file.name<- str_replace(file.name, "_trimming_report.txt", "")
  data.frame("File"=file.name, 
             "Total_Reads"=reads.processed, 
             "Reads_with_Adapters"=reads.with.adapters, 
             "Reads_Passing"=reads.passing.filters)
}
trimming.summary <- do.call(rbind, lapply(trimming.summary.files, extract.trimming.info))

create.xlsx(trimming.summary, file.name = "report/_qc/trimming_report.xlsx")


#### FASTQC report ####
cat("Checking FastQC\n")
fastqc.summary.files <- list.files(path = "report/FASTQC", pattern = "summary.txt$", full.names = T, recursive = T)
fastqc.data <- do.call(rbind, lapply(fastqc.summary.files, read.table, sep="\t"))
colnames(fastqc.data) <- c("Outcome", "Measure", "Sample")
fastqc.check <- fastqc.data %>%
  dplyr::filter(Measure %in% c("Basic Statistics", "Adapter Content", "Per base sequence quality" )) %>%
  tidyr::pivot_wider(id_cols = Sample, names_from = Measure, values_from = Outcome)

create.xlsx(fastqc.check, file.name = "report/_qc/FASTQC_report.xlsx")


#### Mapping efficiencies ####
cat("Checking mapping\n")
extract.pct <- function(x){
  x <- stringr::str_extract(x, "\\(.*\\)")
  x <- stringr::str_replace(x, "\\(", "")
  x <- stringr::str_replace(x, "%\\)", "")
  as.numeric(x)
}

extract.val <- function(x){
  as.numeric(stringr::str_replace(x, " \\(.*\\)" , ""))
}

# Extract the mapping summary from stdout files
# cat bash.o* | grep -w -e 'mapping' -e 'Aligned' -e 'rate' | tr -d '\t' > report/mapping.txt
system2("cat", "logs/bash.o* | grep -w -e 'mapping' -e 'Aligned' -e 'rate' | tr -d '\t' > report/_qc/mapping.txt")

# Spread to separate columns
map.data <- read.table("report/_qc/mapping.txt", sep="$") %>% #  sep char does not exist, force single column
  tidyr::extract(V1, c("Run"), "([S|E]RR\\d+)", remove = FALSE) %>% # column from regex
  tidyr::fill(Run, .direction="down") %>% # fill missing run values
  tidyr::extract(V1, c("Measure", "Reads", "Pct"), "(.*): ?(\\d+)? \\(?(\\d+\\.\\d+)%\\)?", remove = FALSE, convert = TRUE) %>% # find columns from regex
  dplyr::mutate(Measure = str_replace_all(Measure, " ", "_")) %>% # ensure colnames will not have spaces
  dplyr::filter(! str_detect(V1, "mapping")) %>% # remove rows with just 'SRRxxxx mapping'
  dplyr::group_by(Run, Measure) %>%
  dplyr::slice_tail(n=1) %>% # if a sample has been mapped more than once, take only the most recent
  tidyr::pivot_wider(id_cols = Run, names_from = Measure, values_from = c(Reads, Pct)) %>% # make new columns
  dplyr::mutate(Reads_Overall_alignment_rate = rowSums(across(dplyr::starts_with("Reads_")), na.rm=TRUE)) %>%
  merge(., filtered.samples, by="Run")   # Merge in the sample info

# Plot the mapping efficiencies

plot.single.end.mapping <- function(map.data){
  
  p2 <- ggplot(map.data, aes(x = CommonName, y = Pct_Aligned_0_time, col=Timepoint))+
    geom_beeswarm(size=1)+
    labs(y="Percentage of reads (%)", col="Timepoint", title="Unmapped")+
    coord_cartesian(ylim = c(0, 100))+
    facet_wrap(~Organism_part)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_blank())
  
  p3 <- ggplot(map.data, aes(x = CommonName, y = Pct_Aligned_1_time, col=Timepoint))+
    geom_beeswarm(size=1)+
    labs(y="Percentage of reads (%)", col="Timepoint", title="Uniquely mapped")+
    coord_cartesian(ylim = c(0, 100))+
    facet_wrap(~Organism_part)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_blank())
  
  p4 <- ggplot(map.data, aes(x = CommonName, y = `Pct_Aligned_>1_times`, col=Timepoint))+
    geom_beeswarm(size=1)+
    labs(y="Percentage of reads (%)", col="Timepoint", title="Multiple mapped")+
    coord_cartesian(ylim = c(0, 100))+
    facet_wrap(~Organism_part)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_blank())
  
  p5 <- ggplot(map.data, aes(x = CommonName, y = Reads_Overall_alignment_rate, col=Timepoint))+
    geom_beeswarm(size=1)+
    labs(y="Number of reads", col="Timepoint", title="Total reads")+
    facet_wrap(~Organism_part)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_blank())
  
  
  p1 <- ggplot(map.data, aes(x = CommonName, y = Pct_Overall_alignment_rate, col=Timepoint))+
    geom_beeswarm(size=1)+
    labs(y="Overall mapping (%)", col="Timepoint", title="Overall mapping")+
    theme_bw()+
    facet_wrap(~Organism_part)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_blank())
  
  
  p2 + p3 +p4 + patchwork::plot_layout(guides = "collect", axes = "collect") & theme(legend.position = "bottom")
  
  ggsave(plot=last_plot(), filename = "report/_qc/mapping.qc.pct.se.png", dpi=300, units="mm",
         width=200, height = 170)
  
  
  p1 +p5 + patchwork::plot_layout(guides = "collect")  & theme(legend.position = "bottom")
  ggsave(plot=last_plot(), filename = "report/_qc/mapping.qc.total.se.png", dpi=300, units="mm",
         width=200, height = 170)
}

plot.paired.end.mapping <- function(map.data){
  
  p2 <- ggplot(map.data, aes(x = CommonName, y = Pct_Aligned_concordantly_or_discordantly_0_time, col=Timepoint))+
    geom_beeswarm(size=1)+
    labs(y="Percentage of reads (%)", col="Timepoint", title="Unmapped")+
    coord_cartesian(ylim = c(0, 100))+
    facet_wrap(~Organism_part)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_blank())
  
  p3 <- ggplot(map.data, aes(x = CommonName, y = Pct_Aligned_concordantly_1_time, col=Timepoint))+
    geom_beeswarm(size=1)+
    labs(y="Percentage of reads (%)", col="Timepoint", title="Uniquely mapped")+
    coord_cartesian(ylim = c(0, 100))+
    facet_wrap(~Organism_part)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_blank())
  
  p4 <- ggplot(map.data, aes(x = CommonName, y = `Pct_Aligned_concordantly_>1_times`, col=Timepoint))+
    geom_beeswarm(size=1)+
    labs(y="Percentage of reads (%)", col="Timepoint", title="Multiple mapped")+
    coord_cartesian(ylim = c(0, 100))+
    facet_wrap(~Organism_part)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_blank())
  
  p5 <- ggplot(map.data, aes(x = CommonName, y = Reads_Overall_alignment_rate, col=Timepoint))+
    geom_beeswarm(size=1)+
    labs(y="Number of reads", col="Timepoint", title="Total reads")+
    facet_wrap(~Organism_part)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_blank())
  
  
  p1 <- ggplot(map.data, aes(x = CommonName, y = Pct_Overall_alignment_rate, col=Timepoint))+
    geom_beeswarm(size=1)+
    labs(y="Overall mapping (%)", col="Timepoint", title="Overall mapping")+
    theme_bw()+
    facet_wrap(~Organism_part)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_blank())
  
  
  p2 + p3 +p4 + patchwork::plot_layout(guides = "collect", axes = "collect") & theme(legend.position = "bottom")
  
  ggsave(plot=last_plot(), filename = "report/_qc/mapping.qc.pct.pe.png", dpi=300, units="mm",
         width=200, height = 170)
  
  
  p1 + p5 + patchwork::plot_layout(guides = "collect")  & theme(legend.position = "bottom")
  ggsave(plot=last_plot(), filename = "report/_qc/mapping.qc.total.pe.png", dpi=300, units="mm",
         width=200, height = 170)
}

plot.single.end.mapping(map.data[map.data$LibraryLayout=="SINGLE",])
plot.paired.end.mapping(map.data[map.data$LibraryLayout=="PAIRED",])
#### ####

