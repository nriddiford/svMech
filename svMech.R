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

lseq <- function(from=0.1, to=1000, length.out=5) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
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


slideTheme <- function(base_size = 25){
  theme(
    plot.title = element_text(hjust = 0.5, size = 50),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.text = element_text(size=30),
    axis.title = element_text(size=50),
    strip.text = element_text(size=25)
  )
}

setCols <- function(df, col, fill='Y',set="Pastel1"){
  names<-levels(as.factor(df[[col]]))
  names<-sort(names)
  cat("Setting colour levels:", names, "\n")
  level_number<-length(names)
  mycols <- gg_color_hue(level_number)
  # mycols<-brewer.pal(level_number, set)
  names(mycols) <- names
  fillScale <- scale_fill_manual(name = col,values = mycols)
  colScale <- scale_colour_manual(name = col,values = mycols)
  
  if(fill == 'Y') return(fillScale)
  if(fill == 'N') return(colScale)
}

parseDels <- function(infile='data/deletions.txt'){
  deletions <- read.delim(infile, header = T)
  deletions$Inserted_bases <- ifelse(deletions$Inserted_bases == "None", '', as.character(deletions$Inserted_bases))
  deletions$inslength<-nchar(as.character(deletions$Inserted_bases))
  deletions$log10length <- as.numeric(log10(deletions$Length*1000))
  deletions$mhlength<-nchar(as.character(deletions$microhomolgy))
  
  deletions$insertion <- ifelse(deletions$inslength>0,1,0)
  deletions$mh<-ifelse(nchar(as.character(deletions$microhomolgy))>0,1,0)
  
  
  dir.create(file.path("plots"), showWarnings = FALSE)
  
  
  return(deletions)
}

mechStats <- function(){
  dels <- parseDels()

  ins <- dels %>%
    group_by(insertion) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  cat(ins$freq[2], "% deletions have sequence insertions\n")
  
  mhLen <- dels %>%
    group_by(mh) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  cat(mhLen$freq[2], "% deletions have at least 1 bp microhomology\n")
  
  dels$insMH <- ifelse(dels$mh & dels$insertion == 1, 1,0)
  
  coincidence <- dels %>%
    group_by(insMH) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  cat(coincidence$freq[2], "% deletions have at least 1 bp microhomology and insertion\n")
  
  
  
}

delSize <- function(){
  dels <- parseDels()
  
  cols<-setCols(dels, "mh")
  
  dels <- transform(dels, Sample = reorder(Sample, log10length))
  
  p <- ggplot(dels, aes(Sample, log10length))
  p <- p + geom_bar(aes(fill=as.factor(mh)),stat='identity')
  # p <- p + scale_y_continuous("Size (Kb)", labels = lseq())
  p <- p + scale_y_continuous("Size (Kb)", breaks=seq(2,as.integer(max(dels$log10length)),by=1),  labels = lseq())
  
  p <- p + slideTheme() + 
    theme(
      panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 45, hjust=1)
    )
  p <- p + cols
  
  
  dels_out<-paste("Del_lengths.png")
  cat("Writing file", dels_out, "\n")
  ggsave(paste("plots/", dels_out, sep=""), width = 20, height = 10)
  p
  
}

sizeDist <- function(){
  dels <- parseDels()
  cols<-setCols(dels, "Mechanism")
  
  dels <- transform(dels, Mechanism = reorder(Mechanism, log10length))
  
  p <- ggplot(dels, aes(Mechanism, log10length))
  p <- p + geom_violin(aes(fill=Mechanism),alpha=0.6)
  p <- p + geom_jitter(width=0.2)
  # scale_y_log10(breaks=c(0.1,1,100,1000,10000)) 
  p <- p + scale_y_continuous("Size (Kb)", labels = lseq(from=0.1, to=10000, length.out=6))
  p <- p + slideTheme() +
    theme(
      axis.title.x=element_blank(),
      panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted")
    )
  p <- p + cols
  
  mech_out<-paste("Mech_lengths.pdf")
  cat("Writing file", mech_out, "\n")
  ggsave(paste("plots/", mech_out, sep=""), width = 15, height = 10)
  p
  
}


microhomologyLength <- function(){
  dels <- parseDels()
  
  p1 <- ggplot(dels)
  p1 <- p1 + geom_bar(aes(mhlength, (..count..)),stat='count')
  p1 <- p1 + scale_y_continuous("Count",limits = c(0, 20), breaks=seq(0,20,by=2), expand = c(0.01, 0.01))
  p1 <- p1 + scale_x_continuous("MH length")
  
  p1 <- p1 + slideTheme() +
    theme(
      panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted")
    )
  
  dels$extendedMH<-nchar(as.character(dels$extended_homology))
  
  
  p2 <- ggplot(dels)
  p2 <- p2 + geom_bar(aes(extendedMH),stat='count')
  p2 <- p2 + scale_y_continuous("Count",limits = c(0, 20), breaks=seq(0,20,by=2), expand = c(0.01, 0.01))

  p2 <- p2 + scale_x_continuous("Extended hom length", breaks=seq(min(dels$extendedMH),max(dels$extendedMH), by=1))
  
  p2 <- p2 + slideTheme() +
    theme(
      panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
      axis.title.y=element_blank()
    )
  
  mhCount <- dels[dels$mhlength>0,] %>%
    group_by(microhomolgy, mhlength) %>%
    summarise(n=n()) %>%
    ungroup() %>%
    transform(microhomolgy = reorder(microhomolgy, -n))
  
  p3 <- ggplot(mhCount)
  p3 <- p3 + geom_bar(aes(microhomolgy, n ), stat='identity')
  p3 <- p3 + scale_y_continuous("Count", limits = c(0, 10), expand = c(0.01, 0.01), breaks=seq(0,10,by=2))
  
  p3 <- p3 + scale_x_discrete("Microhomology")
  
  p3 <- p3 + slideTheme() +
    theme(
      panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)
    )
  # p3 <- p3 + coord_flip()
  

  extMhCount <- dels[nchar(as.character(dels$extended_homology))>0,] %>%
    group_by(extended_homology, nchar(as.character(dels$extended_homology))) %>%
    summarise(n=n()) %>%
    ungroup() %>%
    transform(extended_homology = reorder(extended_homology, -n))
  
  p4 <- ggplot(extMhCount)
  p4 <- p4 + geom_bar(aes(extended_homology, n ),stat='identity')
  p4 <- p4 + scale_y_continuous("Count",limits = c(0, 10), expand = c(0.01, 0.01), breaks=seq(0,10,by=2))
  
  p4 <- p4 + scale_x_discrete("Extended homology")
  
  p4 <- p4 + slideTheme() +
    theme(
      panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
      axis.title.y=element_blank(),
      axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5,size=20)
      
    )
  
  # p4 <- p4 + coord_flip()
  

  g <- ggarrange(p1, p2, p3, p4,
            ncol = 2, 
            nrow = 2)

  MH_out<-paste("MH_lengths.png")
  cat("Writing file", MH_out, "\n")
  ggsave(paste("plots/", MH_out, sep=""), width = 20, height = 15)
  
  g
  
}



mhLengthInssize <- function(){
  dels <- parseDels()

  # cor(dels$mhlength, dels$inslength)
  # 
  
  
  p <- ggplot(dels)
  p <- p + geom_bar(aes(insertion, (..count..), fill=as.factor(mh)), stat='count', position = 'stack')
  p
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
  # p <- p + scale_y_continuous("Size (Kb)", labels = lseq())
  p <- p + cleanTheme() +
    theme(axis.text.x=element_blank())
  # p <- p + facet_wrap(~Mechanism)
  p
}
