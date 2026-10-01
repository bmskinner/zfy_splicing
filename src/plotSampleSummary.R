#!/bin/Rscript
source("src/functions.R")

#### Make summary tables ####

cat("Making sample summary tables\n")

# What are the timepoints, tissues and species we can look at?
sample.groups <- SELECTED.SAMPLES %>%
  dplyr::rename(
    OriginalTimepoint = DevStage,
    MappedTimepoint = Timepoint
  ) %>%
  dplyr::group_by(Organism, CommonName, MappedTimepoint, Tissue, Sex) %>%
  dplyr::summarise(
    count = n(), TotalBases = sum(Bases),
    .groups = "drop_last"
  ) %>%
  dplyr::mutate(BaseSizeGroup = case_when(TotalBases < 1e10 ~ "Poor",
    TotalBases < 5e10 ~ "OK",
    .default = "Good"
  )) %>%
  dplyr::arrange(CommonName) %>%
  dplyr::ungroup()

# Export summary tables
fs::dir_create("report")
create.xlsx(SELECTED.SAMPLES, "report/analysed.samples.xlsx")
create.xlsx(sample.groups, "report/sample.groups.xlsx")

# Make summary plot of total bases
sample.plot <- ggplot(
  sample.groups %>% dplyr::filter(MappedTimepoint %in% c("adult", "mid-meiosis", "birth")),
  aes(x = interaction(CommonName, Sex), y = TotalBases / 1e9, fill = BaseSizeGroup)
) +
  geom_hline(yintercept = 10, col = "lightgreen") +
  geom_hline(yintercept = 50, col = "darkgreen") +
  geom_col() +
  scale_y_log10() +
  scale_size_manual(values = c(1, 3), guide = "none") +
  scale_fill_manual(values = c("Poor" = "salmon", "OK" = "lightgreen", "Good" = "darkgreen")) +
  labs(y = "Total bases (Gb)") +
  facet_grid(Tissue ~ MappedTimepoint) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    legend.position = "none"
  )
save.double.width("report/read.depths.png", sample.plot, height = 230)

# And the mouse specific timepoints
mouse.samples <- SELECTED.SAMPLES |>
  dplyr::filter(CommonName == "mouse" & str_starts(Timepoint, "Day")) |>
  dplyr::group_by(CommonName, Timepoint, Tissue) %>%
  dplyr::summarise(
    count = n(), TotalBases = sum(Bases),
    .groups = "drop_last"
  ) %>%
  dplyr::mutate(BaseSizeGroup = case_when(TotalBases < 1e10 ~ "Poor",
    TotalBases < 5e10 ~ "OK",
    .default = "Good"
  )) %>%
  dplyr::arrange(CommonName) %>%
  dplyr::ungroup()

mouse.plot <- ggplot(
  mouse.samples,
  aes(x = Timepoint, y = TotalBases / 1e9, fill = BaseSizeGroup)
) +
  geom_hline(yintercept = 10, col = "lightgreen") +
  geom_hline(yintercept = 50, col = "darkgreen") +
  geom_col() +
  scale_y_log10() +
  scale_size_manual(values = c(1, 3), guide = "none") +
  scale_fill_manual(values = c("Poor" = "salmon", "OK" = "lightgreen", "Good" = "darkgreen")) +
  labs(y = "Total bases (Gb)") +
  facet_wrap(~Tissue) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    legend.position = "none"
  )
save.double.width("report/read.depths.mouse.png", mouse.plot, height = 230)

cat("Sample selection: Done!\n")
