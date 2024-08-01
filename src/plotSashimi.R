# Make sashimi plot based on ggsashimi output.
# ggsashimi.py modified to write data objects to Rds when run
# This custom sashimi plot ensures the transcripts are always left to right 
# irrespective of strand
library(tidyverse)
library(data.table)
library(patchwork)
library(grid)
library(scales)
source("src/functions.R")

gtf.data <- read.gtf.data()
# gtf.data <- list(chicken=rtracklayer::import("genomes/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf"))

# Given a canonical transcript id, find the splice junctions
get.canonical.junctions <- function(transcript.id, gtf.data){
  exons <- gtf.data[gtf.data$type=="exon" & gtf.data$transcript_id==transcript.id]
  
  if( all(exons$strand=="-") ) {
    return( na.omit(data.frame(
      j1 = lead(GenomicRanges::end(exons))+1, # end of intron is start of next exon (end coordinate rev strand)
      j2 = GenomicRanges::start(exons)        # first base of intron is start of exon coordinate (rev strand)
    )))  
  }
  
  data.frame(
    j1 = GenomicRanges::end(exons),         # start of the intron is the end of the exon
    j2 = lead(GenomicRanges::start(exons))-1) # end of the intron is the start of the next exon
}

# given a splice junction, determine if this is canonical. i.e. is the junction
# found in the canonical transcript?
# is.junction.canonical("ENSGALT00010007119",chicken.gtf, 118299992 ,118301364)
is.junction.canonical <- function(transcript.id, gtf.data, start, end){
  junctions <- get.canonical.junctions(transcript.id,gtf.data)
  any(junctions$j1==start & junctions$j2==end)
}

# Create the panel of transcripts
make.gene.track <- function(sashimi.data){
  cat("Making gene track\n")
  data <- sashimi.data$density_list
  junctions <- sashimi.data$junction_list
  anns <- sashimi.data$ann_list
  
  # Set coordinates for the x axis
  is.minus.strand <- any(anns$exons$strand=="-")
  xmin <- min(data$x)
  xmax <- max(data$x)
  
  xtmp <- xmin
  xmin <- ifelse(is.minus.strand, xmax, xmin)
  xmax <- ifelse(is.minus.strand, xtmp, xmax)
  
  # Create the arrows to draw on introns
  make.tx.arrows <- function(){
    # Create data table for strand arrows
    txarrows = data.table()
    introns = sashimi.data$ann_list[['introns']]
    # Add right-pointing arrows for plus strand
    if ("+" %in% introns$strand && nrow(introns[strand=="+" & end-start>5, ]) > 0) {
      txarrows = rbind(
        txarrows,
        introns[strand=="+" & end-start>5, list(
          seq(start+4,end,by=453)-1,
          seq(start+4,end,by=453)
        ), by=.(tx,start,end)
        ]
      )
    }
    # Add left-pointing arrows for minus strand
    if ("-" %in% introns$strand && nrow(introns[strand=="-" & end-start>5, ]) > 0) {
      txarrows = rbind(
        txarrows,
        introns[strand=="-" & end-start>5, list(
          seq(start,max(start+1, end-4), by=453),
          seq(start,max(start+1, end-4), by=453)-1
        ), by=.(tx,start,end)
        ]
      )
    }
    txarrows
  }
  
  # Make the gene track
  ggplot()+
    # Introns and arrows
    geom_segment(data=sashimi.data$ann_list$introns, aes(x=start, xend=end, y=tx, yend=tx), linewidth=0.3)+
    geom_segment(data=make.tx.arrows(), aes(x=V1,xend=V2,y=tx,yend=tx), arrow=arrow(length=unit(0.02,"npc")))+
    
    # Exons
    geom_segment(data=sashimi.data$ann_list$exons, aes(x=start, xend=end, y=tx, yend=tx), size=5, alpha=1)+
    
    scale_y_discrete(expand=c(0,0.5))+
    coord_cartesian(xlim = c(xmin,xmax))+
    scale_x_continuous(expand=c(0,0.25), labels=comma)+
    theme_minimal()+
    theme(axis.line.y = element_blank(),
          axis.line.x = element_line(),
          axis.ticks.x = element_line(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank())
}

# Create a sashimi panel for the given ggsashimi data
# min.spanning.reads - the minimum number of junction-spanning reads needed to display on the chart
make.sashimi.panel <- function(sashimi.data, min.spanning.reads=5, label="tissue"){
  cat("Making sashimi panel\n")
  data <- sashimi.data$density_list
  junctions <- sashimi.data$junction_list
  anns <- sashimi.data$ann_list

  # Set coordinates for the x axis
  is.minus.strand <- any(anns$exons$strand=="-")
  xmin <- min(data$x)
  xmax <- max(data$x)
  
  xtmp <- xmin
  xmin <- ifelse(is.minus.strand, xmax, xmin)
  xmax <- ifelse(is.minus.strand, xtmp, xmax)

  ymax <- max(data$y)
  
  junctions <- junctions %>% 
    dplyr::filter(count>=min.spanning.reads)

  if(nrow(junctions)>0){
    # Calculate charting coordinates
    junctions$total.count <- sum(junctions$count)
    junctions$isEven <- sapply(1:nrow(junctions), \(x) x%%2==0)
  }
  
  
  # Create the coverage plot
  splot <- ggplot() + 
    geom_bar(data = data, aes(x, y), width=1, position='identity', stat='identity', alpha=0.5)+
    coord_cartesian(ylim = c(0-(ymax*0.5), ymax*1.5), 
                    xlim = c(xmin, xmax))+
    scale_x_continuous(expand=c(0,0.25))+
    labs(y = sashimi.data[[label]])
  
  
  add.junction <- function(splot, xmin, xmax, ymin, ymax, is.even, count, is.canonical){
    
    # Define the spline shapes that make the junction lines
    spline.color <- ifelse(is.canonical, "grey", "black")
    spline.size  <- ifelse(is.canonical, 1, 2)
    
    l.spline.btm <- grid::xsplineGrob(x=c(0, 0, 1, 1), y=c(1, 0, 0, 0), shape=1, gp=gpar(lwd=spline.size, col=spline.color))
    l.spline.top <- grid::xsplineGrob(x=c(0, 0, 1, 1), y=c(0, 1, 1, 1), shape=1, gp=gpar(lwd=spline.size, col=spline.color))
    r.spline.btm <- grid::xsplineGrob(x=c(1, 1, 0, 0), y=c(1, 0, 0, 0), shape=1, gp=gpar(lwd=spline.size, col=spline.color))
    r.spline.top <- grid::xsplineGrob(x=c(1, 1, 0, 0), y=c(0, 1, 1, 1), shape=1, gp=gpar(lwd=spline.size, col=spline.color))
    
    # Determine which splines to use for minus strand versus plus strand transcripts
    l.grob.btm <- l.spline.btm 
    if(is.minus.strand) l.grob.btm <- r.spline.btm
    
    l.grob.top <- l.spline.top
    if(is.minus.strand) l.grob.top <- r.spline.top
    
    r.grob.btm <- r.spline.btm 
    if(is.minus.strand) r.grob.btm <- l.spline.btm
    
    r.grob.top <- r.spline.top
    if(is.minus.strand) r.grob.top <- l.spline.top
    
    xmid <- (xmin + xmax) / 2
    ymid <- ((ymin + ymax) / 2) * 1.25

    # Left arc
    if(is.even){ # Junctions below zero
      splot <- splot+annotation_custom(grob = l.grob.btm, 
                                       xmin = xmin, 
                                       xmax = xmid, 
                                       ymin = -ymid*0.25, 
                                       ymax = 0)
    } else {
      splot <- splot+annotation_custom(grob = l.grob.top, 
                                       xmin = xmin, 
                                       xmax = xmid, 
                                       ymax = ymid, 
                                       ymin = 0)
    }
    
    # Right arc    
    if(is.even){ # Junctions below zero
      splot <- splot+annotation_custom(grob = r.grob.btm, 
                                       xmin = xmid, 
                                       xmax = xmax, 
                                       ymin = -ymid*0.25, 
                                       ymax = 0)
    } else {
      splot <- splot+annotation_custom(grob = r.grob.top, 
                                       xmin = xmid, 
                                       xmax = xmax, 
                                       ymax = ymid, 
                                       ymin = 0)
    }
    
    splot <- splot + annotate("label", x = xmid, 
                              y = ifelse(is.even, -ymid*0.25, ymid), 
                              label = as.character(count),
                              size=2)
    
    return(splot)
  }
  
  # Add the junctions, adjusting for plus vs minus strand
  if(nrow(junctions)>0){
    for(i in 1:nrow(junctions)){
      splot <- add.junction(splot, junctions[i,]$x, junctions[i,]$xend, junctions[i,]$y, 
                            junctions[i,]$yend, junctions[i,]$isEven, junctions[i,]$count,
                            junctions[i,]$is.canonical)
    }
  }
  
  # Format the final plot
  splot <- splot +
    
    theme_minimal()+
    theme(axis.line = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  list(plot=splot, junctions=junctions)
}

# Read a junction file and note if each junction is in the canonical transcript
read.rds.file <- function(rds.file){
  rds.data <- readRDS(rds.file)
  
  file.name.parts <- str_split_1(basename(rds.file), "\\.")
  rds.data$species <- file.name.parts[1]
  rds.data$tissue <- file.name.parts[2]
  rds.data$timepoint <- file.name.parts[3]
  rds.data$gene.id <- file.name.parts[4]
  rds.data$filename <- basename(rds.file)
  
  rds.data$density_list <- rds.data$density_list[[1]]
  rds.data$junction_list <- rds.data$junction_list[[1]]
  
  canonical.transcript.id <- GENE.LOCATIONS[GENE.LOCATIONS$CommonName==rds.data$species & GENE.LOCATIONS$EnsemblId==rds.data$gene.id,"CanonicalTranscript"]
  
  rds.data$junction_list$is.canonical <- mapply(is.junction.canonical, 
                                                start = rds.data$junction_list$x, 
                                                end   = rds.data$junction_list$xend,
                                                MoreArgs = list(transcript.id = canonical.transcript.id,
                                                                gtf.data      = gtf.data[[rds.data$species]]))
  
  
  rds.data
}

# Create a sashimi panel plot for all tissues of the given species and timepoint
make.species.panel <- function(species, timepoint, gene.id){
  data.files <- list.files(path = "data/merged", pattern = paste0(species, ".*", timepoint, ".*", gene.id, ".*.Rds"), full.names = TRUE)
  if(length(data.files)==0) return()
  data <- lapply(data.files, read.rds.file)

  plots <- lapply(data, \(x) make.sashimi.panel(x, min.spanning.reads = 10,  label="tissue")$plot)
  track <- make.gene.track(data[[1]]) # only one gene, only need one track
  plots[[length(plots)+1]] <- track
  
  patchwork::wrap_plots(plots, nrow = length(plots))
  ggsave(plot = last_plot(), filename = paste0("report/species/", species, ".", timepoint, ".", gene.id, ".png"), dpi = 300, units = "mm", width = 170, height = 240)
}

# Create a sashimi panel plot for all species of the given tissue and timepoint
make.tissue.panel <- function(tissue, timepoint){
  data.files <- list.files(path = "data/merged", pattern = paste0(".*", tissue, ".*", timepoint, ".*.Rds"), full.names = TRUE)
  if(length(data.files)==0) return()
  data <- lapply(data.files, read.rds.file)
  
  plots <- lapply(data, \(x) make.sashimi.panel(x, label="species")$plot)
  tracks <- lapply(data, make.gene.track)
  plots <- c(rbind(plots, tracks))
  
  patchwork::wrap_plots(plots, nrow = length(plots))
  ggsave(plot = last_plot(), filename = paste0("report/tissues/", tissue, ".", timepoint, ".png"), dpi = 300, units = "mm", width = 170, height = 240)
}

# Create a sashimi panel plot for all timepoint of the given tissue and species
make.timepoint.panel <- function(species, tissue, gene.id){
  data.files <- list.files(path = "data/merged", pattern = paste0(species, ".*", tissue, ".*", gene.id, ".*.Rds"), full.names = TRUE)
  if(length(data.files)==0) return()
  
  # Ensure files are plotted in time order
  ordered.files <- list()
  for(t in TIME.ORDER){
    for(f in data.files){
      if(str_detect(f, t)){
        ordered.files <- c(ordered.files, f)
      }
    }
  }
  
  data <- lapply(ordered.files, read.rds.file)

  plots <- lapply(data, \(x) make.sashimi.panel(x, label="timepoint")$plot)
  track <- make.gene.track(data[[1]])
  plots[[length(plots)+1]] <- track
  
  patchwork::wrap_plots(plots, nrow = length(plots))
  ggsave(plot = last_plot(), filename = paste0("report/timepoints/", species, ".", tissue, ".", gene.id, ".png"), dpi = 300, units = "mm", width = 170, height = 240)
}

# find all species combinations for plotting
all.samples <- merge(make.sample.groups(), GENE.LOCATIONS, by="CommonName") 

species.groups <- all.samples %>% 
  dplyr::group_by(CommonName, Timepoint, EnsemblId) %>%
  dplyr::summarise(SampleCount = n())

# Make species plots showing variation over tissues as a specific timepoint
mapply(make.species.panel, species.groups$CommonName, species.groups$Timepoint,species.groups$EnsemblId)

tissue.groups <- all.samples %>% 
  dplyr::group_by(Organism_part, Timepoint) %>%
  dplyr::summarise(SampleCount = n())

# Make tissue plots showing variation over species as a specific timepoint
mapply(make.tissue.panel, tissue.groups$Organism_part, tissue.groups$Timepoint)

timepoint.groups <- all.samples %>% 
  dplyr::group_by(CommonName, Organism_part, EnsemblId) %>%
  dplyr::summarise(SampleCount = n())

# Make timpoint plots showing variation over times in a specific tissue
mapply(make.timepoint.panel, timepoint.groups$CommonName, timepoint.groups$Organism_part, timepoint.groups$EnsemblId)
# make.timepoint.panel("chicken", "testis", "ENSGALG00010003052")
