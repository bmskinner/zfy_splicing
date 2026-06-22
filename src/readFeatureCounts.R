# Read featureCounts outputs, combine per species and calculate TPM for genes
source("src/functions.R")

sample.groups <- SELECTED.SAMPLES %>%
  dplyr::group_by(Organism, Organism_part, Timepoint, CommonName) |>
  dplyr::summarise(Count = n(), .groups = "drop_last")


feature.data <- SELECTED.SAMPLES |>
  dplyr::mutate(featureCountsFile = paste0("data/", CommonName, "/", Run, ".counts.txt"))

read.feature.count.file <- function(file) {
  cat("Reading feature counts from", file, "\n")
  read_tsv(file,
    skip = 2, show_col_types = FALSE,
    col_names = c("GeneId", "Chr", "Start", "End", "Strand", "Length", "Reads"),
    col_types = "cccccii"
  ) |>
    dplyr::mutate(
      Run = str_remove(file, ".counts.txt"),
      Run = str_remove(Run, ".*/"),
      LengthKb = Length / 1000,
      RPK = Reads / LengthKb,
      ScaleFactor = sum(RPK) / 1e6,
      TPM = RPK / ScaleFactor
    ) |>
    dplyr::filter(GeneId %in% GENE.LOCATIONS$GeneId) |>
    dplyr::select(Run, GeneId, Length, Reads, TPM)
}


feature.files <- list.files(path = "data", pattern = ".*.counts.txt$", full.names = TRUE, recursive = TRUE)

feature.values <- do.call(rbind, parallel::mclapply(feature.files,
  read.feature.count.file,
  mc.cores = ifelse(installr::is.windows(), 1, 6)
))

write_csv(feature.values, "report/tpm.csv", quote = "needed")
create.xlsx(feature.values, "report/tpm.xlsx")

feature.values <- read_csv("report/tpm.csv", show_col_types = FALSE)

feature.values <- feature.values |>
  merge(SELECTED.SAMPLES, by = "Run") |>
  merge(GENE.LOCATIONS, by = c("GeneId", "CommonName", "GTF_FILE")) |>
  merge(GENOME.DATA, by = c("CommonName", "Genome", "GTF_FILE")) |>
  dplyr::mutate(
    Clade = fct_relevel(Clade, "Outgroup", "Birds", "Monotremes", "Marsupials", "Artiodactyls", "Primates", "Rodents"),
    Timepoint = fct_relevel(as.factor(Timepoint), "birth", "mid-meiosis", "adult", "Day_00-06", "Day_07-13", "Day_14-20", "Day_21-27")
  ) |>
  dplyr::mutate(Group = case_when(str_detect(Gene, "Z[Ff][Xx]") ~ "ZFX",
    str_detect(Gene, "Z[Ff][Yy]") ~ "ZFY",
    str_detect(Gene, "R[Bb][Mm][Yy]") ~ "RBMY",
    .default = "NA"
  )) |>
  dplyr::group_by(Organism_part, Timepoint, CommonName, Group) |>
  dplyr::mutate(MedianTPM = median(TPM), nSamples = n())


# Standard timepoints
plt <- ggplot(feature.values |> dplyr::filter(Timepoint %in% c("adult", "mid-meiosis", "birth"), nSamples > 1), aes(x = TPM, y = interaction(Gene, CommonName), , fill = Group)) +
  geom_boxplot() +
  facet_wrap(Organism_part ~ Timepoint, scales = "free_y") +
  theme_bw()


plt <- ggplot(feature.values, aes(x = TPM, y = interaction(CommonName, Gene), , fill = Group)) +
  geom_boxplot() +
  facet_wrap(~Timepoint, scales = "free_y") +
  theme_bw()


save.double.width("report/tpm.png", plt)

#### Look at ZFX / ZFY ####
plt <- ggplot(feature.values |> dplyr::filter(Group != "RBMY"), aes(x = Organism_part, y = interaction(Timepoint, CommonName, Gene), fill = MedianTPM)) +
  geom_tile() +
  labs(x = "Tissue", y = "Gene", fill = "Median TPM") +
  facet_wrap(~Clade, scales = "free_y", ncol = 2) +
  scale_fill_viridis_c() +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

save.plot("report/tpm_clade.png", plt, width = 170, height = 250)

#### Look at RBMY ####

plt <- ggplot(feature.values |> dplyr::filter(Group == "RBMY"), aes(x = Organism_part, y = interaction(Timepoint, CommonName, Gene), fill = MedianTPM)) +
  geom_tile() +
  labs(x = "Tissue", y = "Gene", fill = "Median TPM") +
  facet_wrap(~Clade, scales = "free_y", ncol = 2) +
  scale_fill_viridis_c() +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

save.plot("report/tpm_rbmy_clade.png", plt, width = 170, height = 250)
