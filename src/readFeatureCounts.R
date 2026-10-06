# Read featureCounts outputs, combine per species and calculate TPM for genes
source("src/functions.R")

dir.create("report/tpm")

sample.groups <- SELECTED.SAMPLES %>%
  dplyr::group_by(Organism, Tissue, Timepoint, CommonName, Sex) |>
  dplyr::summarise(
    Count = n(),
    TotalBases = sum(Bases),
    .groups = "drop_last"
  )

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

if(length(feature.files)==0) { stop("No feature counts files available to process") }

feature.values <- do.call(rbind, parallel::mclapply(feature.files,
  read.feature.count.file,
  mc.cores = DEFAULT.MC.CORES
))

readr::write_csv(feature.values, "report/tpm.csv", quote = "needed")

#### Read the TPM values ####

feature.values <- readr::read_csv("report/tpm.csv", show_col_types = FALSE) |>
  merge(SELECTED.SAMPLES, by = "Run") |>
  merge(GENE.LOCATIONS, by = c("GeneId", "CommonName", "GTF_FILE")) |>
  merge(GENOME.DATA, by = c("CommonName", "Genome", "GTF_FILE")) |>
  dplyr::mutate(
    Clade = fct_relevel(Clade, "Outgroup", "Birds", "Monotremes", "Marsupials", "Artiodactyls", "Primates", "Rodents"),
    Timepoint = fct_relevel(as.factor(Timepoint), "birth", "mid-meiosis", "adult", "Day_00-06", "Day_07-13", "Day_14-20", "Day_21-27")
  ) |>
  dplyr::group_by(Tissue, Timepoint, CommonName, Group, Sex, GeneId) |>
  dplyr::mutate(
    MedianTPM = median(TPM, na.rm = TRUE),
    MeanTPM = mean(TPM, na.rm = TRUE), nSamples = n(),
    TotalBases = sum(Bases)
  ) |>
  dplyr::ungroup()


##### Standard timepoints #####

###### Tissue by tissue ######
for (tissue in unique(feature.values$Tissue)) {
  tissue.data <- feature.values[feature.values$Tissue == tissue, ]

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
    n.rows <- length(unique(interaction(zfxy.data$Sex, zfxy.data$CommonName, zfxy.data$Gene)))

    plt <- ggplot(
      zfxy.data,
      aes(x = TPM, y = interaction(Gene, Sex, CommonName), fill = Group, col = Group)
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
    save.double.width(paste0("report/tpm/tissue/", tissue, ".zfxy.tpm.png"), plt, height = 50 + n.rows * 5)
  }

  # Mouse specific timepoints
  if (nrow(zfxy.mouse.data) > 0) {
    n.rows <- length(unique(interaction(zfxy.mouse.data$Sex, zfxy.mouse.data$Tissue, zfxy.mouse.data$Gene)))
    plt <- ggplot(
      zfxy.mouse.data,
      aes(x = TPM, y = interaction(Gene, Sex, CommonName), fill = Group, col = Group)
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
    save.double.width(paste0("report/tpm/tissue/", tissue, ".zfxy.mouse.tpm.png"), plt, height = 50 + n.rows * 5)
  }

  # RBMX/Y plot by clade
  if (nrow(rbmy.data) > 0) {
    n.rows <- length(unique(interaction(rbmy.data$Sex, rbmy.data$Tissue, rbmy.data$Gene)))
    plt <- ggplot(
      rbmy.data,
      aes(x = TPM, y = interaction(Gene, Sex, CommonName), fill = Group, col = Group)
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
    save.double.width(paste0("report/tpm/tissue/", tissue, ".rbmxy.tpm.png"), plt, height = 50 + n.rows * 5)
  }

  # Mouse specific timepoints
  if (nrow(rbmy.mouse.data) > 0) {
    n.rows <- length(unique(interaction(rbmy.mouse.data$Sex, rbmy.mouse.data$Tissue, rbmy.mouse.data$Gene)))
    plt <- ggplot(
      rbmy.mouse.data,
      aes(x = TPM, y = interaction(Gene, Sex, CommonName), fill = Group, col = Group)
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
    save.double.width(paste0("report/tpm/tissue/", tissue, ".rbmy.mouse.tpm.png"), plt, height = 50 + n.rows * 5)
  }
}


###### Species by species ######

for (species in unique(feature.values$CommonName)) {
  species.data <- feature.values[feature.values$CommonName == species, ]

  zfxy.data <- species.data |>
    dplyr::filter(
      Timepoint %in% c("adult", "mid-meiosis", "birth"),
      Group %in% c("ZFX", "ZFY", "RBMY")
    )
  zfxy.mouse.data <- species.data |>
    dplyr::filter(
      Timepoint %in% c("Day_00-06", "Day_07-13", "Day_14-20", "Day_21-27"),
      Group %in% c("ZFX", "ZFY", "RBMY")
    )

  # ZFX/Y plot by clade
  if (nrow(zfxy.data) > 0) {
    n.rows <- length(unique(interaction(zfxy.data$Sex, zfxy.data$Tissue, zfxy.data$Gene)))

    plt <- ggplot(
      zfxy.data,
      aes(x = TPM, y = interaction(Gene, Sex), fill = Group, col = Group)
    ) +
      geom_boxplot() +
      scale_fill_manual(values = c("ZFX" = "lightgreen", "ZFY" = "lightblue", "RBMX" = "#F4EA56", "RBMY" = "#c8a2c8")) +
      scale_color_manual(values = c("ZFX" = "darkgreen", "ZFY" = "darkblue", "RBMX" = "#F6BE00", "RBMY" = "purple")) +
      facet_grid(interaction(Tissue, sep = "\n") ~ Timepoint, scales = "free_y", space = "free_y") +
      theme_bw() +
      theme(
        axis.title.y = element_blank(),
        legend.position = "top"
      )
    save.double.width(paste0("report/tpm/species/", species, ".tpm.png"), plt, height = 50 + n.rows * 5)
  }

  # Mouse specific timepoints
  if (nrow(zfxy.mouse.data) > 0) {
    n.rows <- length(unique(interaction(zfxy.mouse.data$Sex, zfxy.mouse.data$Tissue, zfxy.mouse.data$Gene)))
    plt <- ggplot(
      zfxy.mouse.data,
      aes(x = TPM, y = interaction(Gene, Sex), fill = Group, col = Group)
    ) +
      geom_boxplot() +
      scale_fill_manual(values = c("ZFX" = "lightgreen", "ZFY" = "lightblue", "RBMX" = "#F4EA56", "RBMY" = "#c8a2c8")) +
      scale_color_manual(values = c("ZFX" = "darkgreen", "ZFY" = "darkblue", "RBMX" = "#F6BE00", "RBMY" = "purple")) +
      facet_grid(interaction(Tissue, sep = "\n") ~ Timepoint, scales = "free_y", space = "free_y") +
      theme_bw() +
      theme(
        axis.title.y = element_blank(),
        legend.position = "top"
      )
    save.double.width(paste0("report/tpm/species/", species, ".weeks.tpm.png"), plt, height = 50 + n.rows * 5)
  }
}

#### Create splicing table ####

# Create an output table for manually filling detected ZFX/Y splicing. Add the
# median TPM for ZFX and ZFY expression.

junction.data <- readr::read_tsv("report/coding_exon_splice_junctions.tsv", show_col_types = FALSE) |>
  dplyr::filter(E1E3>0) |>
  merge(SELECTED.SAMPLES, by=c("CommonName", "Sex", "Tissue", "Timepoint", "GeneId")) |>
  merge(GENE.LOCATIONS, by=c("CommonName", "GeneId", "GTF_FILE")) |>
  dplyr::select(CommonName, Sex, Tissue, Timepoint, GeneId, Gene, E1E2, E2E3, E1E3, pctE2Spliced, Group) |>
  dplyr::group_by(CommonName, Sex, Tissue, Timepoint, Group) |>
  dplyr::distinct() |>
  dplyr::mutate(Splicing = paste0(Gene, " (", E1E3, " reads, ", sprintf("%.2f%%)",  pctE2Spliced)),
                SplicingDetected = case_when( E1E3>=5 & pctE2Spliced>=5 ~ "Clear",
                                              E1E3>=2 & pctE2Spliced>=10 ~ "Marginal",
                                              E1E3>=1 & pctE2Spliced>=15 ~ "Marginal",
                                              E1E3>=5 & pctE2Spliced>=4 ~ "Marginal",
                                              E1E2 < 10 | E2E3 < 10 ~ "Insufficient coverage",
                                              .default = "")) |>
  dplyr::summarise(Splice_junctions = paste(Splicing, collapse = ", "),
                   Splice_detected = case_when( any(str_detect(SplicingDetected, "Clear")) ~"Clear",
                                                any(str_detect(SplicingDetected, "Marginal")) ~"Marginal",
                                                any(str_detect(SplicingDetected, "Insufficient coverage")) ~"Insufficient coverage",
                                                .default = ""
                                                ),
                   .groups = "drop_last") |>
  tidyr::pivot_wider(names_from = Group, values_from = c(Splice_junctions, Splice_detected)) |>
  dplyr::arrange(CommonName, Tissue, Sex)

splicing.table <- feature.values |>
  dplyr::select(CommonName, Sex, Timepoint, Tissue, Group, GeneId, MedianTPM, TotalBases) |>
  dplyr::distinct() |>
  dplyr::group_by(CommonName, Sex, Timepoint, Tissue, Group) |>
  dplyr::summarise(
    MaxGroupTPM = max(MedianTPM),
    TotalBases = unique(TotalBases),
    .groups = "drop_last"
  ) |>
  tidyr::pivot_wider(names_from = Group, values_from = MaxGroupTPM) |>
  merge(junction.data,
    by = c("CommonName", "Sex", "Timepoint", "Tissue"),
    all.x = TRUE
  ) |>
  merge(GENOME.DATA, by=c("CommonName")) |>
  dplyr::select(Clade, CommonName,
    Sex, Timepoint, Tissue,
    Max_ZFX_TPM = ZFX, Max_ZFY_TPM = ZFY,
    Splice_junctions_ZFX, Splice_junctions_ZFY, Splice_detected_ZFX, Splice_detected_ZFY
  ) |>
  dplyr::mutate(
    Clade = fct_relevel(Clade, "Outgroup", "Birds", "Monotremes", "Marsupials", "Artiodactyls", "Primates", "Rodents"),
    Timepoint = fct_relevel(as.factor(Timepoint), "birth", "mid-meiosis", "adult", "Day_00-06", "Day_07-13", "Day_14-20", "Day_21-27")
  ) |>
  dplyr::arrange(Clade, CommonName, Tissue)


info <- data.frame(
  Column = c("TotalBases", "Max_<gene>_TPM", "Splice_junctions_<gene>","","","" ,"",""),
  Contents = c(
    "The total number of bases in the samples selected for mapping (not the number of mapped reads).",
    "Gene TPM was calculated for each sample. Median TPM was taken across all samples. Paralogues were grouped, and max median TPM was selected, i.e whichever paralogue is most highly expressed",
    "Shows genes with coding exon 2 spliced out. The number of E1E3 junction-spanning reads is shown in parentheses, plus the percentage of E1E3 reads versus the mean of E1E2 and E2E3 spanning reads.",
    "Splicing categories:",
    "Clear: >=5 E1E3 spanning reads and E1E3 reads are >=5% of mean E1E2 + E2E3 spanning reads",
    "Marginal: >=2 E1E3 spanning reads and E1E3 reads are >=10% of mean E1E2 + E2E3 spanning reads",
    "OR >=1 E1E3 spanning reads and E1E3 reads are >=15% of mean E1E2 + E2E3 spanning reads",
    "OR >=5 E1E3 spanning reads and E1E3 reads are >=4% of mean E1E2 + E2E3 spanning reads"
  )
)

sample.summary <- SELECTED.SAMPLES |>
  merge(GENOME.DATA, by = c("CommonName", "Genome", "GTF_FILE")) |>
  dplyr::select(Clade, CommonName, Species, Genome, FASTA_URL, GTF_URL,
    Sex, Timepoint, Tissue,
    Run, DevStage, LibrarySelection, Bases
  ) |>
  dplyr::arrange(Clade, CommonName, Tissue)


expression.summary <- feature.values |>
  dplyr::select(Clade, CommonName,
    Sex, Timepoint, Tissue,
    Run, Gene, CanonicalTranscriptId, Location, TPM, MedianTPM, MeanTPM
  ) |>
  dplyr::arrange(Clade, CommonName, Tissue)

wb <- openxlsx2::wb_workbook() |>
  add.and.freeze(info, "Description") |>
  add.and.freeze(splicing.table, "Splicing summary") |>
  add.and.freeze(sample.summary, "Samples") |>
  add.and.freeze(expression.summary, "Gene Expression") |>
  openxlsx2::wb_add_dxfs_style(
    name = "zfx_splice", font_color = wb_color("darkgreen"),  #"#ccffcc" "#ccccff"
    bg_fill = wb_color("lightgreen")
  ) |>
  openxlsx2::wb_add_dxfs_style(
    name = "zfx_maginal_splice", font_color = wb_color("darkgreen"),
    bg_fill = wb_color("#ccffcc")
  ) |>
  openxlsx2::wb_add_dxfs_style(
    name = "zfy_splice", font_color = wb_color("darkblue"),
    bg_fill = wb_color("#99bbff")
  ) |>
  openxlsx2::wb_add_dxfs_style(
    name = "zfy_maginal_splice", font_color = wb_color("darkblue"),
    bg_fill = wb_color("#bbddff")
  ) |>
  openxlsx2::wb_add_conditional_formatting("Splicing summary", dims = "F1:F200", type = "dataBar", style = c("darkgreen")) |>
  openxlsx2::wb_add_conditional_formatting("Splicing summary", dims = "G1:G200", type = "dataBar", style = c("darkblue")) |>
  openxlsx2::wb_add_conditional_formatting("Splicing summary", dims = "H2:H200", type = "expression",  rule='$J2="Clear"', style = "zfx_splice") |>
  openxlsx2::wb_add_conditional_formatting("Splicing summary", dims = "I2:I200", type = "expression",  rule='$K2="Clear"', style = "zfy_splice") |>
  openxlsx2::wb_add_conditional_formatting("Splicing summary", dims = "H2:H200", type = "expression",  rule='$J2="Marginal"', style = "zfx_maginal_splice") |>
  openxlsx2::wb_add_conditional_formatting("Splicing summary", dims = "I2:I200", type = "expression",  rule='$K2="Marginal"', style = "zfy_maginal_splice") |>
  openxlsx2::wb_save(file = "./report/ZFX_ZFY_splicing_summary_tables.xlsx")


#### Tables to get values for presentations ####

zfy.expn.table <- splicing.table |>
  dplyr::select(Clade:Tissue, Max_ZFY_TPM) |>
  dplyr::mutate(Max_ZFY_TPM = round(Max_ZFY_TPM, digits = 2)) |>
  tidyr::pivot_wider(names_from = Timepoint, values_from = Max_ZFY_TPM) |>
  dplyr::filter(Sex=="male")|>
  dplyr::select(Clade, CommonName, Sex, Tissue, birth, `mid-meiosis`, adult)

zfx.expn.table <- splicing.table |>
  dplyr::select(Clade:Tissue, Max_ZFX_TPM) |>
  dplyr::mutate(Max_ZFX_TPM = round(Max_ZFX_TPM, digits = 2)) |>
  tidyr::pivot_wider(names_from = Timepoint, values_from = Max_ZFX_TPM) |>
  dplyr::filter(Sex=="male")|>
  dplyr::select(Clade, CommonName, Sex, Tissue, birth, `mid-meiosis`, adult)


splicing.zfy.presentation.table <- splicing.table |>
  dplyr::select(Clade:Tissue, Splice_junctions_ZFY) |>
  tidyr::pivot_wider(names_from = Timepoint, values_from = Splice_junctions_ZFY) |>
  dplyr::filter(Sex=="male")|>
  dplyr::select(Clade, CommonName, Sex, Tissue, birth, `mid-meiosis`, adult)

splicing.zfy.presentation.table.mouse <- splicing.table |>
  dplyr::select(Clade:Tissue, Splice_junctions_ZFY) |>
  tidyr::pivot_wider(names_from = Timepoint, values_from = Splice_junctions_ZFY) |>
  dplyr::filter(Sex=="male", CommonName=="mouse") |>
  dplyr::select(Clade, CommonName, Sex, Tissue, `Day_00-06`, `Day_07-13`, `Day_14-20`, `Day_21-27`)

splicing.zfx.presentation.table <- splicing.table |>
  dplyr::select(Clade:Tissue, Splice_junctions_ZFX) |>
  tidyr::pivot_wider(names_from = Timepoint, values_from = Splice_junctions_ZFX) |>
  dplyr::select(Clade, CommonName, Sex, Tissue, birth, `mid-meiosis`, adult)
