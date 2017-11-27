list.of.packages <- c('tidyverse', 'RColorBrewer', 'ggpubr', 'ggalt', 'emdbook')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  cat('Installing missing packages...\n')
  install.packages(new.packages)
}
cat('Silently loading packages...')
suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggpubr))
suppressMessages(library(ggalt))
library("emdbook")

lseq <- function(from=0.1, to=10000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

cleanTheme <- function(base_size = 12){
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.text = element_text(size=12),
    axis.title = element_text(size=30)
  )
}


parseDels <- function(infile='data/deletions.txt'){
  deletions <- read.delim(infile, header = T)
  deletions$insCount<-nchar(as.character(deletions$Inserted_bases))
  deletions$log10length <- as.numeric(log10(deletions$Length*1000))
  
  return(deletions)
}

mechanisms <- function(){
  
  dels <- parseDels()
  p <- ggplot(dels)
  p <- p + geom_histogram(aes(Mechanism, (..count..)), stat='count')
  p
}

sizeGroups <- function(){
  dels <- parseDels()
  
  dels <- transform(dels, Sample = reorder(Sample, log10length))
  
  
  
  p <- ggplot(dels, aes(Sample, log10length))
  p <- p + geom_jitter(aes(colour = Mechanism), size = 3)
  # p <- p + geom_smooth(aes(group = Mechanism, colour=Mechanism), se=F)
  p <- p + geom_encircle(aes(group = Mechanism, colour = Mechanism, fill = Mechanism),alpha = 0.2)
  p <- p + scale_y_continuous("Size (Kb)", labels = lseq())
  p <- p + cleanTheme() +
    theme(axis.text.x=element_blank())
  # p <- p + facet_wrap(~Mechanism)
  p
}


sizeDist <- function(){
  dels <- parseDels()
  dels <- transform(dels, Mechanism = reorder(Mechanism, log10length))
  
  p <- ggplot(dels, aes(Mechanism, log10length))
  p <- p + geom_violin(aes(fill=Mechanism),alpha=0.6)
  p <- p + geom_point()
  p <- p + scale_y_continuous("Size (Kb)", labels = lseq())
  p <- p + cleanTheme() +
    theme(
      axis.title.x=element_blank()
    )
  
  p
  
}