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

feature.values <- feature.values |>
  merge(SELECTED.SAMPLES, by = "Run") |>
  merge(GENE.LOCATIONS, by = c("GeneId", "CommonName", "GTF_FILE")) |>
  dplyr::mutate(Group = case_when(str_detect(Gene, "Z[Ff][Xx]") ~ "ZFX",
    str_detect(Gene, "Z[Ff][Xy]") ~ "ZFY",
    str_detect(Gene, "R[Bb][Mm][Yy]") ~ "RBMY",
    .default = "NA"
  )) |>
  dplyr::group_by(Organism_part, Timepoint, CommonName, Group) |>
  dplyr::mutate(MedianTPM = median(TPM))


plt <- ggplot(feature.values, aes(x = Organism_part, y = TPM, fill = Group)) +
  geom_boxplot() +
  geom_beeswarm() +
  facet_wrap(~CommonName) +
  theme_bw()

save.double.width("report/tpm", plt)
