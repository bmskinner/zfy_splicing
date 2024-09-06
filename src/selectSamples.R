#!/bin/Rscript
library(tidyverse)
library(xlsx)
library(fs)
source("src/functions.R")
# Convert developmental stages to timepoints of interest for ZFX/Y

#### Samples from PRJEB26695 ####

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
  ),
  CommonName = "chicken",
  Genome = "GRCg7b") %>%
  dplyr::filter(Timepoint != "other" & Timepoint != "embryo")  %>%
  dplyr::select(Run, BioProject, LibraryLayout, CommonName, Genome, DevStage, sex, Timepoint, Organism_part, Organism) %>%
  write.csv(., file = "metadata/PRJEB26695.filt.csv", row.names = FALSE, quote = TRUE)


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
  ),
  CommonName = "opossum",
  Genome = "ASM229v1") %>%
  dplyr::filter(Timepoint != "other" & Timepoint != "embryo") %>%
  dplyr::select(Run, BioProject, LibraryLayout, CommonName, Genome, DevStage, sex, Timepoint, Organism_part, Organism) %>%
  write.table(., file = "metadata/PRJEB26695.filt.csv", row.names = FALSE, quote = TRUE, append=TRUE, sep=",", col.names = FALSE)

# Mouse E-MTAB-6798
read.csv("metadata/mouse.csv") %>%
  dplyr::rename(DevStage = Experimental_Factor._developmental_stage..exp.,
                OrganismPart = Experimental_Factor._organism_part..exp.) %>%
  dplyr::mutate(Timepoint = case_when(Developmental_stage == "embryo" ~ "embryo",
                                      DevStage == "postnatal day 0" ~ "birth",
                                      DevStage == "postnatal day 14" ~ "mid-meiosis",
                                      DevStage == "postnatal day 63" ~ "adult",
                                      .default = "other"
  ),
  CommonName = "mouse",
  Genome = "GRCm39") %>%
  dplyr::filter(Timepoint != "other" & Timepoint != "embryo") %>%
  dplyr::select(Run, BioProject, LibraryLayout, CommonName, Genome, DevStage, sex, Timepoint, Organism_part, Organism) %>%
  write.table(., file = "metadata/PRJEB26695.filt.csv", row.names = FALSE, quote = TRUE, append=TRUE, sep=",", col.names = FALSE)

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
  ),
  CommonName = "human",
  Genome = "GRCh38") %>%
  dplyr::filter(Timepoint != "other" & Timepoint != "embryo") %>%
  dplyr::select(Run, BioProject, LibraryLayout, CommonName, Genome, DevStage, sex, Timepoint, Organism_part, Organism) %>%
  write.table(., file = "metadata/PRJEB26695.filt.csv", row.names = FALSE, quote = TRUE, append=TRUE, sep=",", col.names = FALSE)

# Rabbit E-MTAB-6782 - no need to include, there is no Y assembly yet
# read.csv("metadata/rabbit.csv") %>%
#   dplyr::rename(DevStage = Experimental_Factor._developmental_stage..exp.,
#                 OrganismPart = Experimental_Factor._organism_part..exp.) %>%
#     dplyr::mutate(Timepoint = case_when(Developmental_stage == "embryo" ~ "embryo",
#                                         DevStage == "postnatal day 0" ~ "birth",
#                                         DevStage == "postnatal day 14" ~ "mid-meiosis",
#                                         DevStage == "postnatal day 84" ~ "adult",
#                                         DevStage == "postnatal day 186 to 548" ~ "adult",
#                                         .default = "other"
#   ),
# CommonName = "opossum",
# Genome = "GRCg7b") %>%
#   dplyr::filter(Timepoint != "other" & Timepoint != "embryo") %>%
#   dplyr::select(Run, LibraryLayout, DevStage, sex, Timepoint, Organism_part, Organism) %>%
#   write.csv(., file = "metadata/rabbit.filt.csv", row.names = FALSE, quote = TRUE)


# Rat E-MTAB-6811
read.csv("metadata/rat.csv") %>%
  dplyr::rename(DevStage = Experimental_Factor._developmental_stage..exp.,
                OrganismPart = Experimental_Factor._organism_part..exp.) %>%
  dplyr::mutate(Timepoint = case_when(Developmental_stage == "embryo" ~ "embryo",
                                      DevStage == "postnatal day 0" ~ "birth",
                                      DevStage == "postnatal day 14" ~ "mid-meiosis",
                                      DevStage == "postnatal day 112" ~ "adult",
                                      .default = "other"
  ),
  CommonName = "rat",
  Genome = "mRatBN7.2") %>%
  dplyr::filter(Timepoint != "other" & Timepoint != "embryo") %>%
  dplyr::select(Run, BioProject, LibraryLayout, CommonName, Genome, DevStage, sex, Timepoint, Organism_part, Organism) %>%
  write.table(., file = "metadata/PRJEB26695.filt.csv", row.names = FALSE, quote = TRUE, append=TRUE, sep=",", col.names = FALSE)

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
  ),
  CommonName = "macaque",
  Genome = "Mmul_10") %>%
  dplyr::filter(Timepoint != "other" & Timepoint != "embryo") %>%
  dplyr::select(Run, BioProject, LibraryLayout, CommonName, Genome, DevStage, sex, Timepoint, Organism_part, Organism) %>%
  write.table(., file = "metadata/PRJEB26695.filt.csv", row.names = FALSE, quote = TRUE, append=TRUE, sep=",", col.names = FALSE)


#### Samples from PRJEB33381 #### 

read.csv("metadata/PRJEB33381.csv") %>%
  dplyr::rename(OrganismPart = Experimental_Factor._organism_part..exp.,
                Species = Experimental_Factor._organism..exp.,
                DevStage = Developmental_stage) %>%
  dplyr::filter(Species %in% c("Mus musculus", "Monodelphis domestica", "Macaca mulatta", "Rattus norvegicus")) %>%
  dplyr::mutate(Timepoint = DevStage,
                CommonName = case_when(Organism=="Mus musculus"~"mouse",
                                       Organism=="Monodelphis domestica"~"opossum",
                                       Organism=="Macaca mulatta"~"macaque",
                                       Organism=="Rattus norvegicus"~"rat")) %>%
  merge(., get.genome.data(), by="CommonName") %>%
  dplyr::select(Run, BioProject, LibraryLayout, CommonName, Genome, DevStage, sex, Timepoint, Organism_part, Organism) %>%
  write.table(., file = "metadata/PRJEB33381.filt.csv", row.names = FALSE, quote = TRUE, append=FALSE, sep=",", col.names = TRUE)

#### Samples from generic search mouse testis ####

# Goal here is to find non-adult WT mice with known age in days
read.csv("metadata/mouse.testis.csv") %>%
  dplyr::rename(OrganismPart = Experimental_Factor._organism_part..exp.,
                Species = Experimental_Factor._organism..exp.,
                DevStage = Developmental_stage) %>%
  dplyr::filter(Assay.Type == "RNA-Seq" & cell_type == "") %>%
  dplyr::filter(!is.na(AGE) & AGE!="" & AGE!="not collected" & !str_starts(AGE, "E")& AGE!="adult") %>%
  dplyr::filter(!str_detect(OrganismPart, "adipose") & !str_detect(OrganismPart, "brain"))%>%
  dplyr::filter(!str_detect(AGE, "month") & !str_detect(AGE, "week")& !str_detect(AGE, "year")) %>%
  dplyr::filter(!str_detect(tissue, "spermatid") & !str_detect(tissue, "spermatocyte")& !str_detect(tissue, "soermatid")  & !str_detect(source_name, "stem cells")) %>%
  dplyr::mutate(AgeDays = as.numeric(str_extract(AGE, "\\d+"))) %>%
  dplyr::filter(AgeDays < 28) %>%
  dplyr::filter(str_detect(genotype, "[W|w]ild[ |-][T|t]ype") | str_detect(source_name, "[W|w]ild[ |-][T|t]ype") | str_detect(genotype, "[W|w][T|t]")) %>%
  dplyr::mutate(Timepoint = case_when(AgeDays < 7 ~ "Day_00-06",
                                     AgeDays < 15 ~ "Day_07-13",
                                     AgeDays < 21 ~ "Day_14-20",
                                     AgeDays < 28 ~ "Day_21-27",
                                     .default = "Other"),
                CommonName = "mouse") %>%
  # do we need all of them? Just take the 10 smallest runs in each age group
  dplyr::group_by(Timepoint) %>%
  dplyr::arrange(Bases) %>%
  dplyr::slice_head(n=10) %>%
  dplyr::mutate(Organism_part = "testis",
                sex = "male",
                DevStage = paste0("d", AgeDays)) %>% # ensure all consistent
  merge(., get.genome.data(), by="CommonName") %>%
  dplyr::select(Run, BioProject, LibraryLayout, CommonName, Genome, DevStage, sex, Timepoint, Organism_part, Organism) %>%
  write.table(., file = "metadata/mouse.testis.filt.csv", row.names = FALSE, quote = TRUE, append=FALSE, sep=",", col.names = TRUE)

#### Samples from generic search for platypus RNAseq ####

platypus <- read.csv("metadata/platypus.csv") %>%
  dplyr::rename(OrganismPart = Experimental_Factor._organism_part..exp.,
                Species = Experimental_Factor._organism..exp.) %>%
  dplyr::filter(Assay.Type == "RNA-Seq") %>%
  dplyr::mutate(Organism_part = str_to_lower( case_when(tissue != "" ~ tissue,
                                   OrganismPart != "" ~ OrganismPart,
                                   source_name != "" ~ source_name,
                                   .default = NA)),
                DevStage = str_to_lower( case_when(AGE != "" ~ AGE,
                                              Stage != "" ~ Stage,
                                              .default = NA)),
                Timepoint = DevStage,
                CommonName = "platypus") %>%
  dplyr::filter(DevStage!="") %>%
  dplyr::filter(Organism_part!="fibroblast") %>%
  merge(., get.genome.data(), by="CommonName") %>%
  dplyr::select(Run, BioProject, LibraryLayout, CommonName, Genome, DevStage, sex, Timepoint, Organism_part, Organism) %>%
  write.table(., file = "metadata/platypus.filt.csv", row.names = FALSE, quote = TRUE, append=FALSE, sep=",", col.names = TRUE)


#### Samples from generic search for zebrafinch RNAseq ####

zebrafinch <- read.csv("metadata/zebrafinch.csv") %>%
  dplyr::filter(Assay.Type == "RNA-Seq", tissue=="Testis") %>%
  dplyr::mutate(Organism_part = tissue,
                DevStage = "adult",
                Timepoint = "adult",
                CommonName = "zebrafinch") %>%
  merge(., get.genome.data(), by="CommonName") %>%
  dplyr::select(Run, BioProject, LibraryLayout, CommonName, Genome, DevStage, sex, Timepoint, Organism_part, Organism) %>%
  write.table(., file = "metadata/zebrafinch.filt.csv", row.names = FALSE, quote = TRUE, append=FALSE, sep=",", col.names = TRUE)

#### Make summary tables ####

# Read the filtered samples, match folder names
filtered.samples <- read.filtered.samples()

sample.classifications <- filtered.samples %>% 
  dplyr::group_by(Organism, CommonName, DevStage, Timepoint) %>%
  dplyr::rename(OriginalTimepoint = DevStage,
                MappedTimepoint = Timepoint) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::ungroup()

# Export summary tables
fs::dir_create("report")
create.xlsx(filtered.samples, "report/analysed.samples.xlsx")
create.xlsx(sample.classifications, "report/sample.classifictions.xlsx")
