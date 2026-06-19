# Read featureCounts outputs, combine per species and calculate TPM for genes
source("src/functions.R")

sample.groups <- SELECTED.SAMPLES %>%
  dplyr::group_by(Organism, Organism_part, Timepoint, CommonName) |>
  dplyr::summarise(Count = n())


feature.data <- SELECTED.SAMPLES |>
  dplyr::mutate(featureCountsFile = paste0("data/", CommonName, "/", Run, ".counts.txt"))

read.feature.count.file <- function(file) {
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


feature.files <- list.files(path = "data", pattern = ".*.counts.txt", full.names = TRUE, recursive = TRUE)

feature.values <- do.call(rbind, lapply(feature.files, read.feature.count.file)) |>
  merge(SELECTED.SAMPLES, by = "Run") |>
  merge(GENE.LOCATIONS, by = c("GeneId", "CommonName")) |>
  dplyr::mutate(Group = case_when(str_detect(Gene, "Z[Ff][Xx]") ~ "ZFX",
    str_detect(Gene, "Z[Ff][Xy]") ~ "ZFY",
    str_detect(Gene, "R[Bb][Mm][Yy]") ~ "RBMY",
    .default = "NA"
  )) |>
  dplyr::group_by(Organism_part, CommonName, Group) |>
  dplyr::mutate(MedianTPM = median(TPM))


write_csv(feature.values, "report/tpm.csv", quote = FALSE)
create.xlsx(feature.values, "report/tpm.xlsx")

# ggplot(feature.values |> dplyr::filter(Group == "RBMY"), aes(x = Organism_part, y = TPM, fill = )) +
# geom_boxplot() +
#   geom_beeswarm() +
#   facet_wrap(~CommonName) +
#   theme_bw()
