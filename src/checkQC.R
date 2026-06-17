# Check sample QC
cat("QC check: Beginning\n")
source("src/functions.R")


fs::dir_create("report/_qc")

#### Trimming report ####

cat("QC check: Checking trimming\n")

trimming.summary.files <- list.files(path = "data", pattern = "fastq.gz_trimming_report.txt$", full.names = T, recursive = T)

extract.trimming.info <- function(f) {
  trimming.report <- read_file(f)
  reads.processed <- str_extract(trimming.report, "Total reads processed:\\s+[\\d|,]+")
  reads.processed <- str_replace(reads.processed, "Total reads processed:\\s+", "")

  reads.with.adapters <- str_extract(trimming.report, "Reads with adapters:\\s+[\\d|,|(|)|%| |\\.]+")
  reads.with.adapters <- str_replace(reads.with.adapters, "Reads with adapters:\\s+", "")

  reads.passing.filters <- str_extract(trimming.report, "Reads written \\(passing filters\\):\\s+[\\d|,|(|)|%| |\\.]+")
  reads.passing.filters <- str_replace(reads.passing.filters, "Reads written \\(passing filters\\):\\s+", "")

  file.name <- basename(f)
  file.name <- str_replace(file.name, "_trimming_report.txt", "")
  data.frame(
    "File" = file.name,
    "Total_Reads" = reads.processed,
    "Reads_with_Adapters" = reads.with.adapters,
    "Reads_Passing" = reads.passing.filters
  )
}
trimming.summary <- do.call(rbind, lapply(trimming.summary.files, extract.trimming.info))

create.xlsx(trimming.summary, file.name = "report/_qc/trimming_report.xlsx")


#### FASTQC report ####
cat("QC check: Checking FastQC\n")
fastqc.summary.files <- list.files(path = "report/FASTQC", pattern = "summary.txt$", full.names = T, recursive = T)
fastqc.data <- do.call(rbind, lapply(fastqc.summary.files, read.table, sep = "\t"))
colnames(fastqc.data) <- c("Outcome", "Measure", "Sample")
fastqc.check <- fastqc.data %>%
  dplyr::filter(Measure %in% c("Basic Statistics", "Adapter Content", "Per base sequence quality")) %>%
  tidyr::pivot_wider(id_cols = Sample, names_from = Measure, values_from = Outcome) |>
  dplyr::mutate(Run = stringr::str_extract(Sample, "^([A-Z\\d]+)_", group = 1)) |>
  merge(SELECTED.SAMPLES, by = "Run")

create.xlsx(fastqc.check, file.name = "report/_qc/FASTQC_report.xlsx")


#### Mapping efficiencies ####
cat("QC check: Checking mapping reports\n")
extract.pct <- function(x) {
  x <- stringr::str_extract(x, "\\(.*\\)")
  x <- stringr::str_replace(x, "\\(", "")
  x <- stringr::str_replace(x, "%\\)", "")
  as.numeric(x)
}

extract.val <- function(x) {
  as.numeric(stringr::str_replace(x, " \\(.*\\)", ""))
}

fs::file_delete("report/_qc/mapping.txt")

# Extract the mapping summary from stdout files directed to logs
# if (length(list.files(path = "logs", pattern = "bash.o.*")) > 0) {
#   tryCatch(
#     {
#       system2("cat", "logs/bash.o* | grep -w -e 'mapping' -e 'Aligned' -e 'rate' | tr -d '\t' >> report/_qc/mapping.txt")
#     },
#     error = function(e) warning(e)
#   )
# }
# if (length(list.files(path = "logs", pattern = "mapSamples.sh.o.*")) > 0) {
#   tryCatch(
#     {
#       system2("cat", "logs/mapSamples.sh.o* | grep -w -e 'mapping' -e 'Aligned' -e 'rate' | tr -d '\t' >> report/_qc/mapping.txt")
#     },
#     error = function(e) warning(e)
#   )
# }

if (length(list.files(path = "logs", pattern = "*.mapping.log")) > 0) {
  tryCatch(
    {
      system2("cat", "logs/*.mapping.log | grep -w -e 'mapping' -e 'Aligned' -e 'rate' | tr -d '\t' > report/_qc/mapping.txt")
    },
    error = function(e) warning(e)
  )
}

# Spread to separate columns
map.data <- read.table("report/_qc/mapping.txt", sep = "$") %>% #  sep char does not exist, force single column
  tidyr::extract(V1, c("Run"), "([A-Z]RR\\d+)", remove = FALSE) %>% # column from regex
  tidyr::fill(Run, .direction = "down") %>% # fill missing run values
  tidyr::extract(V1, c("Measure", "Reads", "Pct"), "(.*): ?(\\d+)? \\(?(\\d+\\.\\d+)%\\)?", remove = FALSE, convert = TRUE) %>% # find columns from regex
  dplyr::mutate(Measure = str_replace_all(Measure, " ", "_")) %>% # ensure colnames will not have spaces
  dplyr::filter(!str_detect(V1, "mapping")) %>% # remove rows with just 'SRRxxxx mapping'
  dplyr::group_by(Run, Measure) %>%
  dplyr::slice_tail(n = 1) %>% # if a sample has been mapped more than once, take only the most recent
  tidyr::pivot_wider(id_cols = Run, names_from = Measure, values_from = c(Reads, Pct)) %>% # make new columns
  dplyr::mutate(
    Reads_Overall_alignment_rate = rowSums(across(dplyr::starts_with("Reads_")), na.rm = TRUE),
    Single_mapped_pct = sum(Pct_Aligned_1_time, Pct_Aligned_concordantly_1_time),
    Multi_mapped_pct = sum(`Pct_Aligned_>1_times`, `Pct_Aligned_concordantly_>1_times`),
    Unmapped_pct = sum(Pct_Aligned_0_time, `Pct_Aligned_concordantly_or_discordantly_0_time`)
  ) %>%
  merge(., SELECTED.SAMPLES, by = "Run") # Merge in the sample info

create.xlsx(map.data, "report/_qc/mapping.xlsx")

if (any(map.data$Pct_Overall_alignment_rate < 80)) cat("QC check: Some samples have poor mapping rates\n")

# Plot the mapping efficiencies

plot.mapping.rates <- function(map.data) {
  ggplot(map.data, aes(x = Run)) +
    geom_col(aes(y = Unmapped_pct), fill = "salmon", position = "stack") +
    geom_col(aes(y = Single_mapped_pct), fill = "lightgreen", position = "stack") +
    geom_col(aes(y = Multi_mapped_pct), fill = "orange", position = "stack") +
    labs(y = "Percentage of reads (%)", title = "Mapping groups: Single mapped, multimapped, unmapped") +
    coord_cartesian(ylim = c(0, 100)) +
    facet_wrap(Organism_part ~ CommonName, scales = "free_x") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.title.x = element_blank()
    )

  ggsave(
    plot = last_plot(), filename = "report/_qc/mapping.qc.pct.png", dpi = 300, units = "mm",
    width = 200, height = 230
  )

  p1 <- ggplot(map.data, aes(x = CommonName, y = Pct_Overall_alignment_rate, col = Pct_Overall_alignment_rate > 80)) +
    geom_beeswarm(size = 1) +
    labs(y = "Overall mapping (%)", col = "OK", title = "Overall mapping") +
    scale_color_manual(values = c(`FALSE` = "salmon", `TRUE` = "lightgreen")) +
    theme_bw() +
    facet_wrap(~Organism_part) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.title.x = element_blank()
    )

  ggsave(
    plot = last_plot(), filename = "report/_qc/mapping.qc.total.png", dpi = 300, units = "mm",
    width = 200, height = 170
  )
}

tryCatch(
  {
    cat("QC check: Making mapping plots\n")
    plot.mapping.rates(map.data)
  },
  error = \(e) warning(e)
)

#### ####
cat("QC check: Done!\n")
