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
  merge(SELECTED.SAMPLES, by = "Run", all.y = TRUE)

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

file.remove("report/_qc/mapping.txt")

# Extract the mapping summary from stdout files directed to logs
if (length(list.files(path = "logs", pattern = "*.mapping.log")) > 0) {
  tryCatch(
    {
      # Combine all mapping outputs into one file, with one line per sample
      system("for f in logs/*.mapping.log; do  grep -w -e 'Aligned' -e 'rate' $f | echo $f `tr --delete '\t'`  >> report/_qc/mapping.txt; done")
    },
    error = function(e) warning(e)
  )
}

# Spread to separate columns
map.data <- read.table("report/_qc/mapping.txt", sep = "$") |> #  sep char does not exist, force single column
  tidyr::extract(V1, "Run", "([A-Z]RR\\d+)", remove = FALSE) |> # columns from regex
  tidyr::extract(V1,
    into = c("Aligned_concordantly_or_discordantly_0_time", "Pct_aligned_concordantly_or_discordantly_0_time"),
    regex = "Aligned concordantly or discordantly 0 time: (\\d+) \\(([\\d+\\.]+)%\\)",
    remove = FALSE
  ) |>
  tidyr::extract(V1,
    into = c("Aligned_concordantly_1_time", "Pct_aligned_concordantly_1_time"),
    regex = "Aligned concordantly 1 time: (\\d+) \\(([\\d+\\.]+)%\\)",
    remove = FALSE
  ) |>
  tidyr::extract(V1,
    into = c("Aligned_concordantly_greater_1_time", "Pct_aligned_concordantly_greater_1_time"),
    regex = "Aligned concordantly >1 times: (\\d+) \\(([\\d+\\.]+)%\\)",
    remove = FALSE
  ) |>
  tidyr::extract(V1,
    into = c("Aligned_discordantly_1_time", "Pct_aligned_discordantly_1_time"),
    regex = "Aligned discordantly 1 time: (\\d+) \\(([\\d+\\.]+)%\\)",
    remove = FALSE
  ) |>
  tidyr::extract(V1,
    into = c("Pct_overall_alignment_rate"),
    regex = "Overall alignment rate: ([\\d+\\.]+)%",
    remove = FALSE
  ) |>
  tidyr::extract(V1,
    into = c("Aligned_1_time", "Pct_aligned_1_time"),
    regex = "Aligned 1 time: (\\d+) \\(([\\d+\\.]+)%\\)",
    remove = FALSE
  ) |>
  tidyr::extract(V1,
    into = c("Aligned_0_time", "Pct_aligned_0_time"),
    regex = "Aligned 0 time: (\\d+) \\(([\\d+\\.]+)%\\)",
    remove = FALSE
  ) |>
  tidyr::extract(V1,
    into = c("Aligned_greater_1_time", "Pct_aligned_greater_1_time"),
    regex = "Aligned >1 times: (\\d+) \\(([\\d+\\.]+)%\\)",
    remove = FALSE
  ) |>
  dplyr::mutate(across(contains("time"), as.numeric),
    Pct_overall_alignment_rate = as.numeric(Pct_overall_alignment_rate)
  ) |>
  dplyr::mutate(
    Single_mapped = ifelse(is.na(Aligned_concordantly_1_time),
      Aligned_1_time,
      Aligned_concordantly_1_time
    ),
    Single_mapped_pct = ifelse(is.na(Pct_aligned_concordantly_1_time),
      Pct_aligned_1_time,
      Pct_aligned_concordantly_1_time
    ),
    Multi_mapped_pct = ifelse(is.na(Pct_aligned_concordantly_greater_1_time),
      Pct_aligned_greater_1_time,
      Pct_aligned_concordantly_greater_1_time
    ),
    Unmapped_pct = ifelse(is.na(Pct_aligned_concordantly_or_discordantly_0_time),
      Pct_aligned_0_time,
      Pct_aligned_concordantly_or_discordantly_0_time
    )
  ) |>
  dplyr::select(-V1) |>
  merge(SELECTED.SAMPLES, by = "Run", all.y = TRUE) # Merge in the sample info

create.xlsx(map.data, "report/_qc/mapping.xlsx")

if (any(map.data$Pct_overall_alignment_rate < 80)) cat("QC check: Some samples have poor mapping rates\n")

# Plot the mapping efficiencies

plot.mapping.rates <- function(map.data) {
  ggplot(map.data, aes(x = Run)) +
    geom_col(aes(y = Unmapped_pct + Multi_mapped_pct + Single_mapped_pct), fill = "lightgreen", position = "stack") +
    geom_col(aes(y = Unmapped_pct + Multi_mapped_pct), fill = "orange", position = "stack") +
    geom_col(aes(y = Unmapped_pct), fill = "salmon", position = "stack") +
    labs(y = "Percentage of reads (%)", title = "Mapping groups: Single mapped, multimapped, unmapped") +
    coord_cartesian(ylim = c(0, 100)) +
    facet_wrap(CommonName ~ Tissue, scales = "free_x") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.title.x = element_blank()
    )

  ggsave(
    plot = last_plot(), filename = "report/_qc/mapping.qc.pct.png", dpi = 300, units = "mm",
    width = 300, height = 400
  )

  ggplot(map.data, aes(x = CommonName, y = Pct_overall_alignment_rate, col = Pct_overall_alignment_rate > 80)) +
    geom_beeswarm(size = 1) +
    labs(y = "Overall mapping (%)", col = "OK", title = "Overall mapping") +
    scale_color_manual(values = c(`FALSE` = "salmon", `TRUE` = "lightgreen")) +
    theme_bw() +
    facet_wrap(~Tissue) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.title.x = element_blank()
    )

  ggsave(
    plot = last_plot(), filename = "report/_qc/mapping.qc.total.png", dpi = 300, units = "mm",
    width = 200, height = 170
  )

  total.mapped.bases <- map.data |>
    dplyr::filter(Timepoint %in% c("birth", "mid-meiosis", "adult")) |>
    dplyr::group_by(CommonName, Tissue, Timepoint, Sex) |>
    dplyr::summarise(
      TotalMappedReads = sum(Single_mapped),
      .groups = "drop_last"
    ) |>
    dplyr::mutate(BaseSizeGroup = case_when(TotalMappedReads < 1e8 ~ "Low",
      TotalMappedReads < 4e8 ~ "Mid",
      .default = "High"
    ))

  ggplot(
    total.mapped.bases,
    aes(
      x = CommonName, y = TotalMappedReads / 1e6,
      fill = BaseSizeGroup
    )
  ) +
    geom_hline(yintercept = 100, col = "lightgreen") +
    geom_hline(yintercept = 400, col = "darkgreen") +
    scale_fill_manual(values = c("Low" = "salmon", "Mid" = "lightgreen", "High" = "darkgreen")) +
    geom_col() +
    labs(y = "Total mapped reads (Millions)") +
    theme_bw() +
    facet_grid(Tissue ~ Timepoint) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.title.x = element_blank(),
      legend.position = "none"
    )
  ggsave(
    plot = last_plot(), filename = "report/_qc/mapping.qc.mapped_reads.png", dpi = 300, units = "mm",
    width = 170, height = 170
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
