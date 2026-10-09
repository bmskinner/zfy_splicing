# Group the individual sashimi plots into images so we can look over individual
# and sample variation

source("src/functions.R")

make.aggregate.plot <- function(paths, out.file){
  if(file.exists(out.file)) return(FALSE)
  
  paths <- paths[file.size(paths)>0] # there are some placeholder images of zero size when no reads were in a sample
  
  if(length(paths)==0) return(FALSE)
  
  plots <- lapply(paths, \(x){
    img <- as.raster(readPNG(x))
    rasterGrob(img, interpolate = FALSE, default.units = "mm", width = 170, height = 50)
  })
  
  ggsave(
    out.file,
    arrangeGrob(grobs = plots, nrow = length(plots), ncol = 1),
    dpi = 300, units = "mm", width = 170, height = min(1000, 50 * length(plots))
  )
}

#### Individual runs ####

flog.info("Making aggregate plots of individual samples\n")
png.files <- data.frame(path = list.files(path = "report/raw_sashimi", pattern = "[SDE]RR.*.png$", 
                                          full.names = TRUE, recursive=TRUE)) |>
  dplyr::mutate(file = basename(path)) |>
  tidyr::separate_wider_delim(file,
                              delim = ".", names = c("Run", "GeneId", "condensed", "ext")
  ) |>
  merge(SELECTED.SAMPLES, by = c("Run", "GeneId"))


unique.combos <- png.files |>
  dplyr::group_by(CommonName, Tissue, Timepoint, GeneId, Sex) |>
  dplyr::summarise(nPNGs = n(), paths = list(path), .groups = "drop_last") |>
  dplyr::mutate(out.file = paste0("report/tissues/", CommonName, ".", Sex, ".", Tissue, ".", Timepoint, ".", GeneId , ".png"))

invisible(mcmapply(make.aggregate.plot, unique.combos$paths, unique.combos$out.file))

#### Merged images ####
flog.info("Making aggregate plots of merged species samples\n")
png.files <- data.frame(path = list.files(path = "report/merged_sashimi", pattern = ".*.png$", 
                                          full.names = TRUE, recursive=TRUE)) |>
  dplyr::mutate(file = basename(path)) |>
  tidyr::separate_wider_delim(file,
                              delim = ".", names = c("CommonName", "Tissue", "Timepoint", "Sex", "GeneId", "Gene","condensed", "ext")
  )

unique.merged.species.combos <- png.files |>
  dplyr::group_by(CommonName) |>
  dplyr::summarise(nPNGs = n(), paths = list(path), .groups = "drop_last") |>
  dplyr::mutate(out.file = paste0("report/species/", CommonName, ".png"))

invisible(mcmapply(make.aggregate.plot, unique.merged.species.combos$paths, unique.merged.species.combos$out.file))

flog.info("Making aggregate plots of merged tissues samples\n")
unique.merged.tissues.combos <- png.files |>
  dplyr::group_by(Tissue) |>
  dplyr::summarise(nPNGs = n(), paths = list(path), .groups = "drop_last") |>
  dplyr::mutate(out.file = paste0("report/tissues/", Tissue, ".png"))

invisible(mcmapply(make.aggregate.plot, unique.merged.tissues.combos$paths, unique.merged.tissues.combos$out.file))
