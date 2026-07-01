# Read featureCounts outputs, combine per species and calculate TPM for genes
source("src/functions.R")

sample.groups <- SELECTED.SAMPLES %>%
  dplyr::group_by(Organism, Organism_part, Timepoint, CommonName, sex) |>
  dplyr::summarise(Count = n(), .groups = "drop_last")


#### Extract TPM values from featureCounts output ####

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
      LengthKb = Length / 1e3,
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
  mc.cores = DEFAULT.MC.CORES
))

write_csv(feature.values, "report/tpm.csv", quote = "needed")
create.xlsx(feature.values, "report/tpm.xlsx")


#### Read the TPM values ####

feature.values <- read_csv("report/tpm.csv", show_col_types = FALSE)

feature.values <- feature.values |>
  merge(SELECTED.SAMPLES, by = "Run") |>
  merge(GENE.LOCATIONS, by = c("GeneId", "CommonName", "GTF_FILE")) |>
  merge(GENOME.DATA, by = c("CommonName", "Genome", "GTF_FILE")) |>
  dplyr::mutate(
    Clade = fct_relevel(Clade, "Outgroup", "Birds", "Monotremes", "Marsupials", "Artiodactyls", "Primates", "Rodents"),
    Timepoint = fct_relevel(as.factor(Timepoint), "birth", "mid-meiosis", "adult", "Day_00-06", "Day_07-13", "Day_14-20", "Day_21-27")
  ) |>
  dplyr::group_by(Organism_part, Timepoint, CommonName, Group, sex, GeneId) |>
  dplyr::mutate(
    MedianTPM = median(TPM, na.rm = TRUE),
    MeanTPM = mean(TPM, na.rm = TRUE), nSamples = n()
  )


##### Standard timepoints #####

###### Tissue by tissue ######
for (tissue in unique(feature.values$Organism_part)) {
  tissue.data <- feature.values[feature.values$Organism_part == tissue, ]

  zfxy.data <- tissue.data |>
    dplyr::filter(
      Timepoint %in% c("adult", "mid-meiosis", "birth"),
      Group %in% c("ZFX", "ZFY")
    )
  zfxy.mouse.data <- tissue.data |>
    dplyr::filter(
      Timepoint %in% c("Day_00-06", "Day_07-13", "Day_14-20", "Day_21-27"),
      Group %in% c("ZFX", "ZFY")
    )

  rbmy.data <- tissue.data |>
    dplyr::filter(
      Timepoint %in% c("adult", "mid-meiosis", "birth"),
      Group %in% c("RBMX", "RBMY")
    )
  rbmy.mouse.data <- tissue.data |>
    dplyr::filter(
      Timepoint %in% c("Day_00-06", "Day_07-13", "Day_14-20", "Day_21-27"),
      Group %in% c("RBMX", "RBMY")
    )

  # ZFX/Y plot by clade
  if (nrow(zfxy.data) > 0) {
    n.rows <- length(unique(interaction(zfxy.data$sex, zfxy.data$CommonName, zfxy.data$Gene)))

    plt <- ggplot(
      zfxy.data,
      aes(x = TPM, y = interaction(Gene, sex, CommonName), fill = Group, col = Group)
    ) +
      geom_boxplot() +
      scale_fill_manual(values = c("lightgreen", "lightblue")) +
      scale_color_manual(values = c("darkgreen", "darkblue")) +
      facet_grid(interaction(Clade, sep = "\n") ~ Timepoint, scales = "free_y", space = "free_y") +
      theme_bw() +
      theme(
        axis.title.y = element_blank(),
        legend.position = "top"
      )
    save.double.width(paste0("report/tpm/", tissue, ".zfxy.tpm.png"), plt, height = 50 + n.rows * 5)
  }

  # Mouse specific timepoints
  if (nrow(zfxy.mouse.data) > 0) {
    n.rows <- length(unique(interaction(zfxy.mouse.data$sex, zfxy.mouse.data$CommonName, zfxy.mouse.data$Gene)))
    plt <- ggplot(
      zfxy.mouse.data,
      aes(x = TPM, y = interaction(Gene, sex, CommonName), fill = Group, col = Group)
    ) +
      geom_boxplot() +
      scale_fill_manual(values = c("lightgreen", "lightblue")) +
      scale_color_manual(values = c("darkgreen", "darkblue")) +
      facet_grid(interaction(Clade, sep = "\n") ~ Timepoint, scales = "free_y", space = "free_y") +
      theme_bw() +
      theme(
        axis.title.y = element_blank(),
        legend.position = "top"
      )
    save.double.width(paste0("report/tpm/", tissue, ".zfxy.mouse.tpm.png"), plt, height = 50 + n.rows * 5)
  }

  # RBMX/Y plot by clade
  if (nrow(rbmy.data) > 0) {
    n.rows <- length(unique(interaction(rbmy.data$sex, rbmy.data$CommonName, rbmy.data$Gene)))
    plt <- ggplot(
      rbmy.data,
      aes(x = TPM, y = interaction(Gene, sex, CommonName), fill = Group, col = Group)
    ) +
      geom_boxplot() +
      scale_fill_manual(values = c("lightgreen", "lightblue")) +
      scale_color_manual(values = c("darkgreen", "darkblue")) +
      facet_grid(interaction(Clade, sep = "\n") ~ Timepoint, scales = "free_y", space = "free_y") +
      theme_bw() +
      theme(
        axis.title.y = element_blank(),
        legend.position = "top"
      )
    save.double.width(paste0("report/tpm/", tissue, ".rbmxy.tpm.png"), plt, height = 50 + n.rows * 5)
  }

  # Mouse specific timepoints
  if (nrow(rbmy.mouse.data) > 0) {
    n.rows <- length(unique(interaction(rbmy.mouse.data$sex, rbmy.mouse.data$CommonName, rbmy.mouse.data$Gene)))
    plt <- ggplot(
      rbmy.mouse.data,
      aes(x = TPM, y = interaction(Gene, sex, CommonName), fill = Group, col = Group)
    ) +
      geom_boxplot() +
      scale_fill_manual(values = c("lightgreen", "lightblue")) +
      scale_color_manual(values = c("darkgreen", "darkblue")) +
      facet_grid(interaction(Clade, sep = "\n") ~ Timepoint, scales = "free_y", space = "free_y") +
      theme_bw() +
      theme(
        axis.title.y = element_blank(),
        legend.position = "top"
      )
    save.double.width(paste0("report/tpm/", tissue, ".rbmy.mouse.tpm.png"), plt, height = 50 + n.rows * 5)
  }
}


#### Look at ZFX / ZFY ####
plt <- ggplot(
  feature.values |> dplyr::filter(Group != "RBMY"),
  aes(
    x = Organism_part,
    y = interaction(Timepoint, CommonName, Gene, sex),
    fill = MedianTPM
  )
) +
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

plt <- ggplot(
  feature.values |> dplyr::filter(Group == "RBMY"),
  aes(
    x = Organism_part, y = interaction(Timepoint, CommonName, Gene),
    fill = MedianTPM
  )
) +
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

#### Manual splicing table ####

# Create an output table for manually filling detected ZFX/Y splicing. Add the
# median TPM for RBMX and RBMY expressions.

splicing.table <- feature.values |>
  dplyr::select(CommonName, sex, Timepoint, Organism_part, Group, GeneId, MedianTPM) |>
  dplyr::distinct() |>
  dplyr::group_by(CommonName, sex, Timepoint, Organism_part, Group) |>
  dplyr::summarise(MaxGroupTPM = max(MedianTPM)) |>
  tidyr::pivot_wider(names_from = Group, values_from = MaxGroupTPM) |>
  dplyr::select(CommonName, sex, Timepoint, Organism_part, Max_RBMX_TPM = RBMX, Max_RBMY_TPM = RBMY, Max_ZFX_TPM = ZFX, Max_ZFY_TPM = ZFY) |>
  dplyr::mutate(ZFX_Splice_Junctions = "", ZFY_Splice_Junctions = "")

create.xlsx(splicing.table, "./report/splicing_table.xlsx")

# TODO: fill in the rest of the table manually based on the sashimi data
