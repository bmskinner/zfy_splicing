# Filter metadata from SRA searches to get samples of interest

cat("Sample selection: Reading and filtering sample data\n")

# These are the filtered files used for sample selection
filt.files <- list.files("metadata", pattern = "*.filt.csv", full.names = TRUE)
file.remove(filt.files)

#### Samples from human not from PRJEB26695 ####
# None of the reads from SRR6253462 - SRR6253470 mapped successfully.
# PRJNA597586 is single cell data - skip these too.

# Skip these and select more.
read.csv("metadata/human_generic.csv") |>
  dplyr::filter(
    BioProject == "PRJNA301742",
    Assay.Type == "RNA-Seq",
    source_name == "Total testis",
    Bases > 1e8
  ) |>
  dplyr::mutate(
    DevStage = "adult",
    Timepoint = "adult",
    CommonName = "human",
    Organism_part = "testis",
    sex = "male"
  ) |>
  dplyr::slice_head(n = 10) |> # we don't need all of them
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  write.csv(file = "metadata/human_generic.filt.csv", row.names = FALSE, quote = TRUE)

#### Duck testis ####

read.csv("metadata/duck.csv") |>
  dplyr::filter(
    Assay.Type == "RNA-Seq", Organism == "Anas platyrhynchos",
    str_detect(tissue, "[T|t]estis"), Bases > 5e8
  ) |>
  dplyr::slice_head(n = 5) |> # we don't need all of them
  dplyr::mutate(
    DevStage = "adult",
    Organism_part = "testis",
    Timepoint = "adult",
    CommonName = "duck",
  ) |>
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  write.csv(file = "metadata/duck.filt.csv", row.names = FALSE, quote = TRUE)


#### Turkey testis ####
# From PRJNA597008, 38 weeks

read.csv("metadata/turkey.csv") |>
  dplyr::mutate(
    DevStage = "adult",
    Organism_part = "testis",
    Timepoint = "adult",
    CommonName = "turkey",
  ) |>
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  write.csv(file = "metadata/turkey.filt.csv", row.names = FALSE, quote = TRUE)



#### Samples from PRJEB26695 - chicken E-MTAB-6769 ####

# Chicken E-MTAB-6769
read.csv("metadata/chicken.csv") |>
  dplyr::filter(sex == "male") |>
  dplyr::rename(
    DevStage = Experimental_Factor._developmental_stage..exp.,
    OrganismPart = Experimental_Factor._organism_part..exp.
  ) |>
  dplyr::mutate(
    Timepoint = case_when(Developmental_stage == "embryo" ~ "embryo",
      DevStage == "postnatal day 0" ~ "birth",
      DevStage == "postnatal day 7" ~ "birth",
      DevStage == "postnatal day 70" ~ "mid-meiosis",
      DevStage == "postnatal day 155" ~ "adult",
      .default = "other"
    ),
    CommonName = "chicken"
  ) |>
  dplyr::filter(Timepoint != "other" & Timepoint != "embryo") |>
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  write.csv(file = "metadata/PRJEB26695.filt.csv", row.names = FALSE, quote = TRUE)

#### Samples from PRJEB26695 - opossum E-MTAB-6833 ####
# Opossum E-MTAB-6833
read.csv("metadata/opossum.csv") |>
  dplyr::filter(sex == "male") |>
  dplyr::rename(
    DevStage = Experimental_Factor._developmental_stage..exp.,
    OrganismPart = Experimental_Factor._organism_part..exp.
  ) |>
  dplyr::mutate(
    Timepoint = case_when(Developmental_stage == "embryo" ~ "embryo",
      DevStage == "postnatal day 0" ~ "birth",
      # DevStage == "postnatal day 28" ~ "birth",
      DevStage == "postnatal day 60" ~ "mid-meiosis",
      DevStage == "postnatal day 180" ~ "adult",
      .default = "other"
    ),
    CommonName = "opossum"
  ) |>
  dplyr::filter(Timepoint != "other" & Timepoint != "embryo") |>
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  write.table(file = "metadata/PRJEB26695.filt.csv", row.names = FALSE, quote = TRUE, append = TRUE, sep = ",", col.names = FALSE)

#### Samples from PRJEB26695 - mouse E-MTAB-6798 ####
# Mouse E-MTAB-6798
read.csv("metadata/mouse.csv") |>
  dplyr::filter(sex == "male") |>
  dplyr::rename(
    DevStage = Experimental_Factor._developmental_stage..exp.,
    OrganismPart = Experimental_Factor._organism_part..exp.
  ) |>
  dplyr::mutate(
    Timepoint = case_when(Developmental_stage == "embryo" ~ "embryo",
      DevStage == "postnatal day 0" ~ "birth",
      DevStage == "postnatal day 14" ~ "mid-meiosis",
      DevStage == "postnatal day 63" ~ "adult",
      .default = "other"
    ),
    CommonName = "mouse"
  ) |>
  dplyr::filter(Timepoint != "other" & Timepoint != "embryo") |>
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  write.table(file = "metadata/PRJEB26695.filt.csv", row.names = FALSE, quote = TRUE, append = TRUE, sep = ",", col.names = FALSE)

#### Samples from PRJEB26695 - human E-MTAB-6814 ####
# Human E-MTAB-6814
read.csv("metadata/human.csv") |>
  dplyr::rename(
    DevStage = Experimental_Factor._developmental_stage..exp.,
    OrganismPart = Experimental_Factor._organism_part..exp.
  ) |>
  dplyr::filter(sex == "male") |>
  dplyr::mutate(
    Timepoint = case_when(Developmental_stage == "embryo" ~ "embryo",
      DevStage == "neonate" ~ "birth",
      DevStage == "adolescent" ~ "mid-meiosis",
      DevStage == "middle adult" ~ "adult",
      DevStage == "elderly" ~ "adult",
      .default = "other"
    ),
    CommonName = "human"
  ) |>
  dplyr::filter(Timepoint != "other" & Timepoint != "embryo") |>
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  write.table(file = "metadata/PRJEB26695.filt.csv", row.names = FALSE, quote = TRUE, append = TRUE, sep = ",", col.names = FALSE)

#### Samples from PRJEB26889 - rat   E-MTAB-6811 ####
# Rat E-MTAB-6811
# read.csv("metadata/rat.csv") |>
#   dplyr::filter(sex == "male") |>
#   dplyr::rename(
#     DevStage = Experimental_Factor._developmental_stage..exp.,
#     OrganismPart = Experimental_Factor._organism_part..exp.
#   ) |>
#   dplyr::mutate(
#     Timepoint = case_when(Developmental_stage == "embryo" ~ "embryo",
#       DevStage == "postnatal day 0" ~ "birth",
#       DevStage == "postnatal day 14" ~ "mid-meiosis",
#       DevStage == "postnatal day 112" ~ "adult",
#       .default = "other"
#     ),
#     CommonName = "rat"
#   ) |>
#   dplyr::filter(Timepoint != "other" & Timepoint != "embryo") |>
#   merge(GENOME.DATA, by = "CommonName") |>
#   dplyr::select(
#     Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
#     Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
#   ) |>
#   write.table(file = "metadata/PRJEB26695.filt.csv", row.names = FALSE, quote = TRUE, append = TRUE, sep = ",", col.names = FALSE)

#### Samples from PRJEB26695 - macaque ####
# Rhesus macacque E-MTAB-6813
# Lifespan can reach up to 40 in captivity, median age 25
read.csv("metadata/macaque.csv") |>
  dplyr::filter(sex == "male") |>
  dplyr::rename(
    DevStage = Experimental_Factor._developmental_stage..exp.,
    OrganismPart = Experimental_Factor._organism_part..exp.
  ) |>
  dplyr::mutate(
    Timepoint = case_when(Developmental_stage == "embryo" ~ "embryo",
      DevStage == "postnatal day 0" ~ "birth",
      DevStage == "3 years postnatal" ~ "mid-meiosis", # adolescence for a macaque
      DevStage == "8 years postnatal" ~ "mid-meiosis", # adolescence for a macaque
      DevStage == "14 to 15 years postnatal" ~ "adult",
      DevStage == "20 to 26 years postnatal" ~ "adult", # reaching menopause in females
      .default = "other"
    ),
    CommonName = "macaque"
  ) |>
  dplyr::filter(Timepoint != "other" & Timepoint != "embryo") |>
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  write.table(file = "metadata/PRJEB26695.filt.csv", row.names = FALSE, quote = TRUE, append = TRUE, sep = ",", col.names = FALSE)


#### Samples from PRJEB33381 - multispecies ####

read.csv("metadata/PRJEB33381.csv") |>
  dplyr::rename(
    OrganismPart = Experimental_Factor._organism_part..exp.,
    Species = Experimental_Factor._organism..exp.,
    DevStage = Developmental_stage
  ) |>
  dplyr::filter(Species %in% c("Mus musculus", "Monodelphis domestica", "Macaca mulatta", "Rattus norvegicus", "Sus scrofa")) |>
  dplyr::mutate(
    Timepoint = DevStage,
    CommonName = case_when(
      Organism == "Mus musculus" ~ "mouse",
      Organism == "Monodelphis domestica" ~ "opossum",
      Organism == "Macaca mulatta" ~ "macaque",
      Organism == "Rattus norvegicus" ~ "rat",
      Organism == "Sus scrofa" ~ "pig"
    )
  ) |>
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  write.table(file = "metadata/PRJEB33381.filt.csv", row.names = FALSE, quote = TRUE, append = FALSE, sep = ",", col.names = TRUE)

#### Samples from PRJNA238328 (rat tissues) ####

# SRA metadata does not contain the ages. We get these from the published study
# The GEO to sample name is from
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53960 The sample name to
# age is in supplementary file 1 in Yu et al 2014:
# https://pmc.ncbi.nlm.nih.gov/articles/PMC4381750/

PRJNA238328.geo.to.sample <- read_tsv("metadata/PRJNA238328.GEO_to_sample.tsv", show_col_types = FALSE) |>
  dplyr::mutate(Sample_Name = stringr::str_remove(Sample_Name, "SEQC_"))
PRJNA238328.sample.sheet <- read_csv("metadata/PRJNA238328.sample.sheet.csv", show_col_types = FALSE)

PRJNA238328 <- read.csv("metadata/PRJNA238328.csv") |>
  merge(PRJNA238328.geo.to.sample, by.x = "GEO_Accession..exp.", by.y = "GEO_Accession", all.x = TRUE) |>
  merge(PRJNA238328.sample.sheet, by.x = "Sample_Name", by.y = "Sample_ID", all.x = TRUE) |>
  dplyr::filter(tissue %in% c("Brain", "Heart", "Kidney", "Lung", "Liver", "Muscle", "Spleen", "Testes")) |>
  dplyr::mutate(
    Organism_part = stringr::str_to_lower(tissue),
    Organism_part = stringr::str_replace(Organism_part, "testes", "testis"),
    DevStage = Age_Week,
    Timepoint = case_when(Age_Week == 2 ~ "mid-meiosis",
      Age_Week == 6 ~ "adult",
      Age_Week == 21 ~ "adult",
      Age_Week == 104 ~ "adult",
      .default = "unknown"
    ),
    CommonName = "rat",
  ) |>
  dplyr::group_by(Organism_part, sex, Timepoint) |>
  dplyr::arrange(desc(Bases)) |>
  dplyr::slice_head(n = 10) |> # we don't need all of them
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  write.table(file = "metadata/PRJNA238328.filt.csv", row.names = FALSE, quote = TRUE, append = FALSE, sep = ",", col.names = TRUE)

#### Samples from PRJEB26889 - E-MTAB-6811 rat ####
# Cardoso Moreira 2019
# Puberty in rats is ~6 weeks. Meiosis genes peak expression  ~P14

read.csv("metadata/PRJEB26889.csv") |>
  dplyr::filter(Developmental_Stage != "embryo") |>
  dplyr::filter(Organism_part %in% c("brain", "heart", "kidney", "lung", "liver", "muscle", "spleen", "testis", "ovary")) |>
  dplyr::mutate(
    DevStage = Experimental_Factor._developmental_stage..exp.,
    Timepoint = case_when(DevStage == "postnatal day 0" ~ "birth",
      DevStage == "postnatal day 3" ~ "birth",
      DevStage == "postnatal day 7" ~ "mid-meiosis",
      DevStage == "postnatal day 14" ~ "mid-meiosis",
      DevStage == "postnatal day 42" ~ "adult",
      DevStage == "postnatal day 112" ~ "adult",
      .default = "other"
    ),
    CommonName = "rat",
  ) |>
  dplyr::filter(Timepoint != "other") |>
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  write.table(file = "metadata/PRJEB26889.filt.csv", row.names = FALSE, quote = TRUE, append = FALSE, sep = ",", col.names = TRUE)



#### Samples from PRJNA889410 - Arvicanthis tissues ####
# Not enough restis coverage - skip these
# read.csv("metadata/PRJNA889410.csv") |>
#   dplyr::filter(tissue %in% c("Brain", "Heart", "Kidney", "Lung", "Liver", "Muscle", "Spleen", "Testis", "Ovary")) |>
#   dplyr::mutate(
#     Organism_part = stringr::str_to_lower(tissue),
#     sex = case_when(Organism_part == "testis" ~ "male",
#       Organism_part == "ovary" ~ "female",
#       .default = "unknown"
#     ),
#     DevStage = "adult",
#     Timepoint = "adult",
#     CommonName = "nilerat",
#   ) |>
#   merge(GENOME.DATA, by = "CommonName") |>
#   dplyr::select(
#     Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
#     Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
#   ) |>
#   write.table(file = "metadata/PRJNA889410.filt.csv", row.names = FALSE, quote = TRUE, append = FALSE, sep = ",", col.names = TRUE)
#


#### Samples from generic search mouse testis ####

# Goal here is to find non-adult WT mice with known age in days
read.csv("metadata/mouse.testis.csv") |>
  dplyr::rename(
    OrganismPart = Experimental_Factor._organism_part..exp.,
    Species = Experimental_Factor._organism..exp.,
    DevStage = Developmental_stage
  ) |>
  dplyr::filter(BioProject != "PRJNA630221") |> # project has dsRNA only
  dplyr::filter(Assay.Type == "RNA-Seq" & cell_type == "") |>
  dplyr::filter(!is.na(AGE) & AGE != "" & AGE != "not collected" & !str_starts(AGE, "E") & AGE != "adult") |>
  dplyr::filter(!str_detect(OrganismPart, "adipose") & !str_detect(OrganismPart, "brain")) |>
  dplyr::filter(!str_detect(AGE, "month") & !str_detect(AGE, "week") & !str_detect(AGE, "year")) |>
  dplyr::filter(!str_detect(tissue, "spermatid") & !str_detect(tissue, "spermatocyte") & !str_detect(tissue, "soermatid") & !str_detect(source_name, "stem cells")) |>
  dplyr::mutate(AgeDays = as.numeric(str_extract(AGE, "\\d+"))) |>
  dplyr::filter(AgeDays < 28) |>
  dplyr::filter(str_detect(genotype, "[W|w]ild[ |-][T|t]ype") | str_detect(source_name, "[W|w]ild[ |-][T|t]ype") | str_detect(genotype, "[W|w][T|t]")) |>
  dplyr::mutate(
    Timepoint = case_when(AgeDays < 7 ~ "Day_00-06",
      AgeDays < 14 ~ "Day_07-13",
      AgeDays < 21 ~ "Day_14-20",
      AgeDays < 28 ~ "Day_21-27",
      .default = "Other"
    ),
    CommonName = "mouse"
  ) |>
  # do we need all of them? Just take the 20 smallest runs in each age group
  dplyr::group_by(Timepoint) |>
  dplyr::arrange(Bases) |>
  dplyr::slice_head(n = 20) |>
  dplyr::mutate(
    Organism_part = "testis",
    sex = "male",
    DevStage = paste0("d", AgeDays)
  ) |> # ensure all consistent
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  write.table(file = "metadata/mouse.testis.filt.csv", row.names = FALSE, quote = TRUE, append = FALSE, sep = ",", col.names = TRUE)

#### Samples from generic search for platypus RNA-seq and specific testis search ####

read.csv("metadata/platypus.csv") |>
  dplyr::rename(
    OrganismPart = Experimental_Factor._organism_part..exp.,
    Species = Experimental_Factor._organism..exp.
  ) |>
  dplyr::filter(sex == "male") |>
  dplyr::filter(
    Assay.Type == "RNA-Seq",
    LibrarySelection != "size fractionation"
  ) |>
  dplyr::mutate(
    Organism_part = str_to_lower(case_when(tissue != "" ~ tissue,
      OrganismPart != "" ~ OrganismPart,
      source_name != "" ~ source_name,
      .default = NA
    )),
    DevStage = str_to_lower(case_when(AGE != "" ~ AGE,
      Stage != "" ~ Stage,
      .default = NA
    )),
    Timepoint = DevStage,
    CommonName = "platypus"
  ) |>
  dplyr::filter(DevStage != "") |>
  dplyr::filter(Organism_part != "fibroblast") |>
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  rbind(
    read.csv("metadata/platypus_testis.csv") |>
      dplyr::filter(sex == "male") |>
      dplyr::filter(
        Assay.Type == "RNA-Seq", Organism_part %in% c("testis", ""),
        Experimental_Factor._protocol..exp. != "Ribo-seq",
        LibrarySelection != "size fractionation"
      ) |>
      dplyr::mutate(
        Organism_part = "testis",
        DevStage = "adult",
        Timepoint = "adult",
        CommonName = "platypus"
      ) |>
      merge(GENOME.DATA, by = "CommonName") |>
      dplyr::select(
        Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
        Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
      )
  ) |>
  dplyr::distinct() |>
  write.table(file = "metadata/platypus.filt.csv", row.names = FALSE, quote = TRUE, append = FALSE, sep = ",", col.names = TRUE)

#### Samples from generic search for zebrafinch testis RNAseq ####

read.csv("metadata/zebrafinch.csv") |>
  dplyr::filter(
    Assay.Type == "RNA-Seq",
    str_detect(tissue, "[T|t]estis") | str_detect(tissue_type, "[T|t]estis")
  ) |>
  dplyr::mutate(
    Organism_part = "testis",
    DevStage = "adult",
    Timepoint = "adult",
    CommonName = "zebrafinch"
  ) |>
  dplyr::arrange(desc(Bases)) |>
  dplyr::slice_head(n = 5) |>
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  write.table(file = "metadata/zebrafinch.filt.csv", row.names = FALSE, quote = TRUE, append = FALSE, sep = ",", col.names = TRUE)

#### Samples from generic search Xenopus tropicalis testis ####

read.csv("metadata/xenopus.csv") |>
  dplyr::filter(Assay.Type == "RNA-Seq", tissue == "testis") |>
  dplyr::mutate(
    Organism_part = str_to_lower(tissue),
    DevStage = "adult",
    Timepoint = "adult",
    CommonName = "xenopus"
  ) |>
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  write.table(file = "metadata/xenopus.filt.csv", row.names = FALSE, quote = TRUE, append = FALSE, sep = ",", col.names = TRUE)

#### Samples from generic search Zebrafish testis ####

# Mature at ~3 months
read.csv("metadata/zebrafish.csv") |>
  dplyr::filter(
    Assay.Type == "RNA-Seq", tissue == "testis", genotype == "wild type"
  ) |>
  dplyr::filter(BioProject != "PRJNA540466") |> # skip known single cell
  dplyr::mutate(
    Organism_part = str_to_lower(tissue),
    sex = "male",
    DevStage = "adult",
    Timepoint = "adult",
    CommonName = "zebrafish"
  ) |>
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  write.table(file = "metadata/zebrafish.filt.csv", row.names = FALSE, quote = TRUE, append = FALSE, sep = ",", col.names = TRUE)

#### Samples from generic search Anole testis ####

read.csv("metadata/anole.csv") |>
  dplyr::filter(
    Assay.Type == "RNA-Seq", str_detect(tissue, "[T|t]estis"),
    str_detect(Stage, "[A|a]dult") | str_detect(dev_stage, "[A|a]dult")
  ) |>
  dplyr::mutate(
    Organism_part = "testis",
    DevStage = "adult",
    Timepoint = "adult",
    CommonName = "anole"
  ) |>
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  write.table(file = "metadata/anole.filt.csv", row.names = FALSE, quote = TRUE, append = FALSE, sep = ",", col.names = TRUE)


#### Samples from other koala tissues in PRJNA230900 ####

# Note that the male samples are a general mixed pool. Only female has the
# tissues explicit. Skip this one.

# koala <- read.csv("metadata/PRJNA230900.csv")

#### Samples from Koala testis ####
# PRJNA1158232 and PRJNA1187226 is a study of miRNA, but performed standard RNA-seq some samples
# Their own mapping efficiencies were similar (Table S7, Yu et al 10.21203/rs.3.rs-5671983/v1)
# The koala ids are given in supplementary table S1 from Y et al Cell. 2025 Mar 7;188(8):2081–2093.e16. doi: 10.1016/j.cell.2025.02.006
read.csv("metadata/PRJNA1158232.csv") |>
  dplyr::mutate(
    Organism_part = stringr::str_to_lower(tissue),
    DevStage = "adult",
    Timepoint = "adult",
    CommonName = "koala",
    sex = case_when(koala_id %in% c("K94283", "K98214") ~ "female",
      koala_id %in% c(
        "K94276", "K71362", "K98314", "K98224", "K63464",
        "K63855", "Cove", "Mario", "K98494", "Andy",
        "Burke", "Poppy"
      ) ~ "male",
      .default = "unknown"
    )
  ) |>
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  rbind(
    # Other testis samples in SRA
    read.csv("metadata/koala.csv") |>
      dplyr::rename(dev_stage = Developmental_Stage) |>
      dplyr::filter(
        Organism == "Phascolarctos cinereus",
        Assay.Type == "RNA-Seq", str_detect(tissue, "[T|t]estis"),
        str_detect(Stage, "[A|a]dult") | str_detect(dev_stage, "[A|a]dult")
      ) |>
      dplyr::mutate(
        Organism_part = "testis",
        DevStage = "adult",
        Timepoint = "adult",
        sex = "male",
        CommonName = "koala"
      ) |>
      merge(GENOME.DATA, by = "CommonName") |>
      dplyr::select(
        Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
        Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
      )
  ) |>
  dplyr::distinct() |> # have PRJNA1158232 samples in the wider testis search
  write.table(file = "metadata/koala.filt.csv", row.names = FALSE, quote = TRUE, append = FALSE, sep = ",", col.names = TRUE)


#### Samples from Echidna testis ####

read.csv("metadata/echidna.csv") |>
  dplyr::filter(
    Organism == "Tachyglossus aculeatus",
    Assay.Type == "RNA-Seq", str_detect(tissue, "[T|t]estis"),
    str_detect(dev_stage, "[A|a]dult")
  ) |>
  dplyr::mutate(
    Organism_part = "testis",
    DevStage = "adult",
    Timepoint = "adult",
    sex = "male",
    CommonName = "echidna"
  ) |>
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  write.table(file = "metadata/echidna.filt.csv", row.names = FALSE, quote = TRUE, append = FALSE, sep = ",", col.names = TRUE)



#### Samples from cattle testis ####

# PRJNA471564 - 2 day old and 18 month old

# PRJNA776655 - three groups - TY0 = prepuberty; TY1 = puberty; TY2 =
# postpuberty.
# bulls are prepuberty (at birth, n = 23), the second group
# represents the bulls are puberty (about 1 year old showing heat for the first
# time, n = 23), and the last group represents the bulls are postpuberty (about
# 2 years of age, n = 23).

read.csv("metadata/cattle.csv") |>
  dplyr::filter(
    Organism == "Bos taurus",
    Assay.Type == "RNA-Seq",
    BioProject %in% c("PRJNA471564", "PRJNA776655"),
  ) |>
  dplyr::mutate(
    Organism_part = "testis",
    Timepoint = case_when(str_detect(Sample.Name, "neonatal") ~ "birth",
      str_detect(Sample.Name, "mature") ~ "adult",
      str_detect(Sample.Name, "TY0") ~ "birth",
      str_detect(Sample.Name, "TY1") ~ "mid-meiosis",
      str_detect(Sample.Name, "TY2") ~ "adult",
      .default = "missing"
    ),
    DevStage = Timepoint,
    CommonName = "cattle"
  ) |>
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  write.table(file = "metadata/cattle.filt.csv", row.names = FALSE, quote = TRUE, append = FALSE, sep = ",", col.names = TRUE)

#### Samples from wallaby ####

# From PRJNA1218892

# Note - tissue is mislabelled. Actual tissue is in the library name
read.csv("metadata/PRJNA1218892.csv") |>
  dplyr::mutate(
    Organism_part = stringr::str_extract(Library.Name, ".*_(\\w+)$", group = 1),
    Organism_part = stringr::str_replace(Organism_part, "testes", "testis"),
    Timepoint = "adult", # assumed - no publication for this!
    DevStage = Timepoint,
    CommonName = "wallaby",
    sex = "male"
  ) |>
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  rbind(
    read.csv("metadata/wallaby.csv") |>
      dplyr::filter(
        Organism == "Notamacropus eugenii",
        BioProject == "PRJDB1934"
      ) |>
      dplyr::mutate(
        Organism_part = "testis",
        Timepoint = "adult", # assumed - no publication for this!
        DevStage = Timepoint,
        CommonName = "wallaby",
        sex = "male"
      ) |>
      merge(GENOME.DATA, by = "CommonName") |>
      dplyr::select(
        Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
        Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
      )
  ) |>
  write.table(file = "metadata/wallaby.filt.csv", row.names = FALSE, quote = TRUE, append = FALSE, sep = ",", col.names = TRUE)

#### Samples from tasmanian devil ####

# Brain samples for tasmanian devil. Mix of male and female

# PRJEB28680 samples are from healthy tissue

read.csv("metadata/PRJEB28680.csv") |>
  dplyr::filter(
    !(tissue_type %in% c("DFT1", "DFT2"))
  ) |>
  dplyr::mutate(
    Organism_part = tissue_type,
    Timepoint = "adult", # assumed, age not given in paper
    DevStage = Timepoint,
    CommonName = "tasmaniandevil",
    sex = ifelse(Organism_part == "testis", "male", "unknown")
  ) |>
  merge(GENOME.DATA, by = "CommonName") |>
  dplyr::select(
    Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
    Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
  ) |>
  rbind(
    # Other brain samples - most likely other tissue to have splicing if present
    read.csv("metadata/tasmanian_devil_brain.csv") |>
      dplyr::filter(
        LibrarySelection != "size fractionation",
        BioProject != "PRJEB28680"
      ) |>
      dplyr::mutate(
        Organism_part = "brain",
        Timepoint = "adult", # assumed - no publication for this!
        DevStage = Timepoint,
        CommonName = "tasmaniandevil"
      ) |>
      merge(GENOME.DATA, by = "CommonName") |>
      dplyr::select(
        Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
        Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
      )
  ) |>
  rbind(
    # SRA run selector down, so manually write equivalent table from ENA
    data.frame(
      Run = c("ERR3568424", "ERR3568434"),
      BioProject = "PRJEB34650",
      LibraryLayout = "PAIRED",
      CommonName = "tasmaniandevil",
      DevStage = "adult",
      sex = "male",
      Timepoint = "adult",
      Organism_part = "testis",
      Organism = "Sarcophilus harrisii",
      LibrarySelection = "cDNA",
      LibrarySource = "TRANSCRIPTOMIC",
      Bases = c(16643508000, 18824106750)
    ) |>
      merge(GENOME.DATA, by = "CommonName") |>
      dplyr::select(
        Run, BioProject, LibraryLayout, CommonName, Genome, GTF_FILE, DevStage, sex,
        Timepoint, Organism_part, Organism, LibrarySelection, LibrarySource, Bases
      )
  ) |>
  write.table(
    file = "metadata/tasmaniandevil.filt.csv", row.names = FALSE, quote = TRUE,
    append = FALSE, sep = ",", col.names = TRUE
  )
