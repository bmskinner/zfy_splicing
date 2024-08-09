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

#### Main functions #### 

gtf.data <- read.gtf.data()
# gtf.data <- list(chicken=rtracklayer::import("genomes/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf"))

is.reverse.strand <- function(transcript.id, gtf.data){
  exons <- gtf.data[gtf.data$type=="exon" & gtf.data$transcript_id==transcript.id]
  any(GenomicRanges::strand(exons)=="-")
}

# Given a canonical transcript id, find the splice junctions
get.canonical.junctions <- function(canonical.transcript.id, gtf.data){
  exons <- gtf.data[gtf.data$type=="exon" & gtf.data$transcript_id==canonical.transcript.id]

  if( is.reverse.strand(canonical.transcript.id, gtf.data) ) {
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
is.junction.canonical <- function(canonical.transcript.id, gtf.data, start, end){
  junctions <- get.canonical.junctions(canonical.transcript.id,gtf.data)
  any(junctions$j1==start & junctions$j2==end)
}

# Create the panel of transcripts
make.gene.track <- function(sashimi.data){
  # cat("Making gene track\n")
  data <- sashimi.data$density_list
  junctions <- sashimi.data$junction_list
  anns <- sashimi.data$ann_list
  
  # Set coordinates for the x axis, reversing if on reverse strand
  xmin <- ifelse(sashimi.data$is.reverse.strand, max(data$x), min(data$x))
  xmax <- ifelse(sashimi.data$is.reverse.strand, min(data$x), max(data$x))
  
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
    geom_segment(data=sashimi.data$ann_list$exons, aes(x=start, xend=end, y=tx, yend=tx, col=tx==sashimi.data$canonical.transcript.id), size=5, alpha=1)+
    scale_color_manual(values = c(`TRUE`="blue", `FALSE`="grey"))+
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
          panel.grid = element_blank(),
          legend.position = "none")
}

# Create a sashimi panel for the given ggsashimi data
# min.spanning.reads - the minimum number of junction-spanning reads needed to display on the chart
# y.label - the data value to be used for the y axis label
make.sashimi.panel <- function(sashimi.data, min.spanning.reads=5, label="tissue"){
  # cat("Making sashimi panel\n")
  data <- sashimi.data$density_list
  junctions <- sashimi.data$junction_list
  anns <- sashimi.data$ann_list

  # Set coordinates for the x axis
  is.minus.strand <- sashimi.data$is.reverse.strand
  xmin <- min(data$x)
  xmax <- max(data$x)
  
  xtmp <- xmin
  xmin <- ifelse(is.minus.strand, xmax, xmin)
  xmax <- ifelse(is.minus.strand, xtmp, xmax)

  ymax <- max(data$y)
  
  # Create the coverage plot
  splot <- ggplot() + 
    geom_bar(data = data, aes(x, y), width=1, position='identity', stat='identity', alpha=0.5)+
    coord_cartesian(ylim = c(0-(ymax), ymax*1.5), 
                    xlim = c(xmin, xmax))+
    scale_x_continuous(expand=c(0,0.25))+
    labs(y = label)
  
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
                                       ymin = -ymid*0.5, 
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
                                       ymin = -ymid*0.5, 
                                       ymax = 0)
    } else {
      splot <- splot+annotation_custom(grob = r.grob.top, 
                                       xmin = xmid, 
                                       xmax = xmax, 
                                       ymax = ymid, 
                                       ymin = 0)
    }
    
    splot <- splot + annotate("label", x = xmid, 
                              y = ifelse(is.even, -ymid*0.5, ymid), 
                              label = as.character(count),
                              size=2)
    
    return(splot)
  }
  
  if(nrow(junctions)>0){
    junctions <- junctions %>% 
      dplyr::filter(count>=min.spanning.reads)
  }
  
  # Add the junctions, adjusting for plus vs minus strand
  if(nrow(junctions)>0){
    
    cat("minus strand: ", is.minus.strand, "\n")
    print(junctions)

    
    # Calculate charting coordinates
    junctions$total.count <- sum(junctions$count)
    junctions$isEven <- sapply(1:nrow(junctions), \(x) x%%2==0)
    
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
  rds.data$gene.name <- GENE.LOCATIONS[GENE.LOCATIONS$EnsemblId==rds.data$gene.id, "Gene"]
  rds.data$filename <- basename(rds.file)
  
  rds.data$density_list <- rds.data$density_list[[1]]
  rds.data$junction_list <- rds.data$junction_list[[1]]
  
  rds.data$canonical.transcript.id <- GENE.LOCATIONS[GENE.LOCATIONS$EnsemblId==rds.data$gene.id, "CanonicalTranscript"]
  rds.data$gtf.data <- gtf.data[[rds.data$species]]
  rds.data$is.reverse.strand <- is.reverse.strand(rds.data$canonical.transcript.id, rds.data$gtf.data)
  
  rds.data$junction_list$is.canonical <- mapply(is.junction.canonical, 
                                                start = rds.data$junction_list$x, 
                                                end   = rds.data$junction_list$xend,
                                                MoreArgs = list(transcript.id = rds.data$canonical.transcript.id,
                                                                gtf.data      = rds.data$gtf.data))
  
  
  rds.data
}

# Create a sashimi panel plot for all tissues of the given species and timepoint
make.species.panel <- function(species, timepoint, gene.id){
  data.files <- list.files(path = "data/merged", pattern = paste0(species, ".*", timepoint, ".*", gene.id, ".*.Rds"), full.names = TRUE)
  if(length(data.files)==0) return()
  data <- lapply(data.files, read.rds.file)
  
  gene.name <- GENE.LOCATIONS[GENE.LOCATIONS$EnsemblId==gene.id,]$Gene
  
  out.png.file <- paste0("report/species/", species, ".", timepoint, ".", gene.id, ".", gene.name, ".png")
  
  # if(!file.exists(out.png.file)){
    plots <- lapply(data, \(x)  make.sashimi.panel(x, min.spanning.reads = 10,  label=x$tissue)$plot)
    track <- make.gene.track(data[[1]]) # only one gene, only need one track
    plots[[length(plots)+1]] <- track
    
    patchwork::wrap_plots(plots, nrow = length(plots))
    ggsave(plot = last_plot(), filename =out.png.file, dpi = 300, units = "mm", width = 170, height = 240)
  # }
}

# Create a sashimi panel plot for all species of the given tissue and timepoint
make.tissue.panel <- function(tissue, timepoint){
  data.files <- list.files(path = "data/merged", pattern = paste0(".*\\.", tissue, "\\..*", timepoint, ".*.Rds"), full.names = TRUE)
  if(length(data.files)==0) return()
  data <- lapply(data.files, read.rds.file)
  
  out.png.file <- paste0("report/tissues/", tissue, ".", timepoint, ".png")
  
  plots <- lapply(data, \(x) make.sashimi.panel(x, label= paste0(x$species, "\n", x$gene.name))$plot)
  tracks <- lapply(data, make.gene.track)
  plots <- c(rbind(plots, tracks))
  
  patchwork::wrap_plots(plots, nrow = length(plots))
  ggsave(plot = last_plot(), filename = out.png.file, dpi = 300, units = "mm", width = 170, height = 240)
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
  
  gene.name <- GENE.LOCATIONS[GENE.LOCATIONS$EnsemblId==gene.id,]$Gene
  
  out.png.file <-  paste0("report/timepoints/", species, ".", tissue, ".", gene.id, ".", gene.name, ".png")

  plots <- lapply(data, \(x) make.sashimi.panel(x, label=x$timepoint)$plot)
  track <- make.gene.track(data[[1]])
  plots[[length(plots)+1]] <- track
  
  patchwork::wrap_plots(plots, nrow = length(plots))
  ggsave(plot = last_plot(), filename =out.png.file, dpi = 300, units = "mm", width = 170, height = 240)
}

# find all species combinations for plotting
all.samples <- merge(make.sample.groups(), GENE.LOCATIONS, by="CommonName") 


#### Generate results figures ####

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

# Make timepoint plots showing variation over times in a specific tissue
mapply(make.timepoint.panel, timepoint.groups$CommonName, timepoint.groups$Organism_part, timepoint.groups$EnsemblId)
# make.timepoint.panel("chicken", "testis", "ENSGALG00010003052")

#### Condense info #### 

# given a gtf file and a gene to be plotted, make every exon the same width and 
# every intron the same width. Apply these new coordinates to the given junctions also.
# However, exons can overlap. Ensure to correct for this.
condense.exon.bounds <- function(gtf.data, gene.id){
  exons <- gtf.data[gtf.data$type=="exon" & gtf.data$gene_id==gene.id]
  
  reduced <- GenomicRanges::shift(GenomicRanges::reduce(exons), shift=-min(GenomicRanges::start(exons)))
  reduced$Group <- 1:length(reduced)
  print(reduced)
  
  exons <- GenomicRanges::shift(exons, shift=-min(GenomicRanges::start(exons)))
  
  exons <- exons %>%
    as.data.frame() %>%
    dplyr::select(start, end,  exon_id, strand) %>%
    # dplyr::mutate(ordered.start = ifelse(strand=="-", end, start),
    #               ordered.end  =  ifelse(strand=="-", start, end) ) %>%

    dplyr::distinct()
    # dplyr::mutate(condensed.start = dplyr::row_number()*2-1,
    #               condensed.end   = condensed.start+1) %>%
    # dplyr::select(exon.start = start, exon.end=end, condensed.start, condensed.end, exon_id, strand)
  exons
}

map.condensed.exon.bounds <- function(junctions, exons){
  
  s <- exons[,c(1, 3)]
  e <- exons[,c(2, 4)]
  
  colnames(s) <- c("original", "condensed")
  colnames(e) <- c("original", "condensed")
  ex <- rbind(s, e)
  
  junctions %>%
    dplyr::mutate(exon.end = x-1, exon.start = xend) %>%
    merge(., ex, by.x = "exon.end", by.y="original" ) %>%

    merge(., ex, by.x = "exon.start", by.y="original" ) %>%
    dplyr::select(-x, -xend) %>%
    dplyr::rename(xend = condensed.x, x=condensed.y)
}

condense.exon.bounds(gtf.data$chicken, "ENSGALG00010003052")
map.condensed.exon.bounds(test[[1]]$junction_list, condense.exon.bounds(gtf.data$chicken, "ENSGALG00010003052"))

make.condensed.panel <- function(sashimi.data, min.spanning.reads=5, label="tissue"){
  cat("Making sashimi panel\n")
  exons <- condense.exon.bounds(gtf.data[[sashimi.data$species]], sashimi.data$gene.id)
  junctions <- map.condensed.exon.bounds(sashimi.data$junction_list, exons)
  anns <- sashimi.data$ann_list
  
  # Set coordinates for the x axis
  is.minus.strand <- any(anns$exons$strand=="-")
  xmin <- 0
  xmax <- max(junctions$x, junctions$xend)
  
  xtmp <- xmin
  xmin <- ifelse(is.minus.strand, xmax, xmin)
  xmax <- ifelse(is.minus.strand, xtmp, xmax)
  
  ymax <- max(junctions$count)
  
  junctions <- junctions %>% 
    dplyr::filter(count>=min.spanning.reads)
  
  if(nrow(junctions)>0){
    # Calculate charting coordinates
    junctions$total.count <- sum(junctions$count)
    junctions$isEven <- sapply(1:nrow(junctions), \(x) x%%2==0)
  }
  
  
  # Create the coverage plot
  splot <- ggplot() + 
    coord_cartesian(ylim = c(0-(ymax*0.5), ymax*1.5), 
                    xlim = c(xmin, xmax))+
    scale_x_continuous(expand=c(0,0.25))+
    labs(y = sashimi.data[[label]])+
    # Draw the exon squares
    geom_rect(data=exons, aes(xmin=condensed.start, xmax=condensed.end, ymin=-5, ymax=5), fill="black")
  
  
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



make.condensed.panel(test[[1]])
