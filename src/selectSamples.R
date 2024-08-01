#!/bin/Rscript
library(tidyverse)
library(xlsx)
source("src/functions.R")
# Convert developmental stages to timepoints of interest for ZFX/Y

# Chicken E-MTAB-6769
read.csv("metadata/chicken.csv") %>%
  dplyr::rename(DevStage = Experimental_Factor._developmental_stage..exp.,
                OrganismPart = Experimental_Factor._organism_part..exp.) %>%
  dplyr::mutate(Timepoint = case_when(Developmental_stage == "embryo" ~ "embryo",
                                      DevStage == "postnatal day 0" ~ "birth",
                                      DevStage == "postnatal day 7" ~ "birth",
                                      DevStage == "postnatal day 70" ~ "mid-meiosis",
                                      DevStage == "postnatal day 155" ~ "adult",
                                      .default = "other"
  )) %>%
  dplyr::filter(Timepoint != "other" & Timepoint != "embryo")  %>%
  dplyr::select(Run, DevStage, sex, Timepoint, Organism_part, Organism) %>%
  write.csv(., file = "metadata/chicken.filt.csv", row.names = FALSE, quote = TRUE)


# Opossum E-MTAB-6833
read.csv("metadata/opossum.csv") %>%
  dplyr::rename(DevStage = Experimental_Factor._developmental_stage..exp.,
                OrganismPart = Experimental_Factor._organism_part..exp.) %>%
  dplyr::mutate(Timepoint = case_when(Developmental_stage == "embryo" ~ "embryo",
                                      DevStage == "postnatal day 0" ~ "birth",
                                      # DevStage == "postnatal day 28" ~ "birth",
                                      DevStage == "postnatal day 60" ~ "mid-meiosis",
                                      DevStage == "postnatal day 180" ~ "adult",
                                      .default = "other"
  )) %>%
  dplyr::filter(Timepoint != "other" & Timepoint != "embryo") %>%
  dplyr::select(Run, DevStage, sex, Timepoint, Organism_part, Organism) %>%
  write.csv(., file = "metadata/opossum.filt.csv", row.names = FALSE, quote = TRUE)

# Mouse E-MTAB-6798
read.csv("metadata/mouse.csv") %>%
  dplyr::rename(DevStage = Experimental_Factor._developmental_stage..exp.,
                OrganismPart = Experimental_Factor._organism_part..exp.) %>%
  dplyr::mutate(Timepoint = case_when(Developmental_stage == "embryo" ~ "embryo",
                                      DevStage == "postnatal day 0" ~ "birth",
                                      DevStage == "postnatal day 14" ~ "mid-meiosis",
                                      DevStage == "postnatal day 63" ~ "adult",
                                      .default = "other"
  )) %>%
  dplyr::filter(Timepoint != "other" & Timepoint != "embryo") %>%
  dplyr::select(Run, DevStage, sex, Timepoint, Organism_part, Organism) %>%
  write.csv(., file = "metadata/mouse.filt.csv", row.names = FALSE, quote = TRUE)

# Human E-MTAB-6814
read.csv("metadata/human.csv") %>%
  dplyr::rename(DevStage = Experimental_Factor._developmental_stage..exp.,
                OrganismPart = Experimental_Factor._organism_part..exp.) %>%
  dplyr::mutate(Timepoint = case_when(Developmental_stage == "embryo" ~ "embryo",
                                      DevStage == "neonate" ~ "birth",
                                      DevStage == "adolescent" ~ "mid-meiosis",
                                      DevStage == "middle adult" ~ "adult",
                                      DevStage == "elderly" ~ "adult",
                                      .default = "other"
  )) %>%
  dplyr::filter(Timepoint != "other" & Timepoint != "embryo") %>%
  dplyr::select(Run, DevStage, sex, Timepoint, Organism_part, Organism) %>%
  write.csv(., file = "metadata/human.filt.csv", row.names = FALSE, quote = TRUE)

# Rabbit E-MTAB-6782
read.csv("metadata/rabbit.csv") %>%
  dplyr::rename(DevStage = Experimental_Factor._developmental_stage..exp.,
                OrganismPart = Experimental_Factor._organism_part..exp.) %>%
    dplyr::mutate(Timepoint = case_when(Developmental_stage == "embryo" ~ "embryo",
                                        DevStage == "postnatal day 0" ~ "birth",
                                        DevStage == "postnatal day 14" ~ "mid-meiosis",
                                        DevStage == "postnatal day 84" ~ "adult",
                                        DevStage == "postnatal day 186 to 548" ~ "adult",
                                        .default = "other"
  )) %>%
  dplyr::filter(Timepoint != "other" & Timepoint != "embryo") %>%
  dplyr::select(Run, DevStage, sex, Timepoint, Organism_part, Organism) %>%
  write.csv(., file = "metadata/rabbit.filt.csv", row.names = FALSE, quote = TRUE)


# Rat E-MTAB-6811
read.csv("metadata/rat.csv") %>%
  dplyr::rename(DevStage = Experimental_Factor._developmental_stage..exp.,
                OrganismPart = Experimental_Factor._organism_part..exp.) %>%
  dplyr::mutate(Timepoint = case_when(Developmental_stage == "embryo" ~ "embryo",
                                      DevStage == "postnatal day 0" ~ "birth",
                                      DevStage == "postnatal day 14" ~ "mid-meiosis",
                                      DevStage == "postnatal day 112" ~ "adult",
                                      .default = "other"
  )) %>%
  dplyr::filter(Timepoint != "other" & Timepoint != "embryo") %>%
  dplyr::select(Run, DevStage, sex, Timepoint, Organism_part, Organism) %>%
  write.csv(., file = "metadata/rat.filt.csv", row.names = FALSE, quote = TRUE)

# Rhesus macacque E-MTAB-6813
# Lifespan can reach up to 40 in captivity, median age 25
read.csv("metadata/macaque.csv") %>%
  dplyr::rename(DevStage = Experimental_Factor._developmental_stage..exp.,
                OrganismPart = Experimental_Factor._organism_part..exp.) %>%
  dplyr::mutate(Timepoint = case_when(Developmental_stage == "embryo" ~ "embryo",
                                      DevStage == "postnatal day 0" ~ "birth",
                                      DevStage == "3 years postnatal" ~ "mid-meiosis", # adolescence for a macaque
                                      DevStage == "8 years postnatal" ~ "mid-meiosis", # adolescence for a macaque
                                      DevStage == "14 to 15 years postnatal" ~ "adult",
                                      DevStage == "20 to 26 years postnatal" ~ "adult", # reaching menopause in females
                                      .default = "other"
  )) %>%
  dplyr::filter(Timepoint != "other" & Timepoint != "embryo") %>%
  dplyr::select(Run, DevStage, sex, Timepoint, Organism_part, Organism) %>%
  write.csv(., file = "metadata/macaque.filt.csv", row.names = FALSE, quote = TRUE)

# Map to the reference genome
# system2("bash", "mapSamples.sh")

# Read the filtered samples, match folder names
filtered.samples <- do.call(rbind, lapply(list.files(path="metadata", pattern = "*.filt.csv", full.names = TRUE), 
                                          \(f) read.csv(f) %>% dplyr::mutate(CommonName = str_replace(str_replace(f, "metadata/", ""), ".filt.csv", ""))))

sample.classifications <- filtered.samples %>% 
  dplyr::group_by(Organism, CommonName, DevStage, Timepoint) %>%
  dplyr::rename(OriginalTimepoint = DevStage,
                MappedTimepoint = Timepoint) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::ungroup()





# Export summary tables
create.xlsx(filtered.samples, "report/analysed.samples.xlsx")
create.xlsx(sample.classifications, "report/sample.classifictions.xlsx")
