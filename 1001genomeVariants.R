##Nemhauser Lab Cytometer Settings
##load in the packages and codes that R needs to analyze cytometer data
install.packages('colorspa')
library(flowViz)
library(flowCore)
library(ggplot2)
library(colorspace)
library(plyr)
##set the directory where you have saved the general flow analysis files

##Set the title of the experiment,
##which should be the name of the folder containing your flow data
source('/Volumes/Nemlings/Cytometry/NemRfunctions/noSummaryAnalysis.R')
source('/Volumes/Nemlings/Cytometry/NemRfunctions/General.R')

experiment<-"20150916_1001genomeVariants"
setwd(paste0("/Volumes/Nemlings/Clay/Cytometer/",experiment))

plate1<-read.flowSet(path=paste(experiment,"_1/",sep=""),alter.names=T)
sampleNames(plate1)<-paste("1_",sampleNames(plate1),sep="")
plate2<-read.flowSet(path=paste(experiment,"_2/",sep=""),alter.names=T)
sampleNames(plate2)<-paste("2_",sampleNames(plate2),sep="")
dat<-rbind2(plate1,plate2)
length(dat)
#remove col12
grepl('12',pData(dat)$name)
dat<-subset(dat, !grepl('12',pData(dat)$name))
sampleNames(dat)
pData(dat)$name<-sampleNames(dat)
pData(dat)

AFB <- c(rep(c("TIR1","TIR1",'AFB2','AFB2',"TIR1-T154S","TIR1-T154S","TIR1-E239K","TIR1-E239K","TIR1-S546L","TIR1-S546L",'AFB2-E204K'),8),rep(c('AFB2-E204K','AFB2-A254N','AFB2-A254N','AFB2-T491R','AFB2-T491R','AFB2-Q169L','AFB2-Q169L','AFB2-D176E','AFB2-D176E','AFB2-T179M','AFB2-T179M'),8))
IAA<-c(rep(c(rep(c('IAA1','IAA17'),5),'IAA1'),8),rep(c(rep(c('IAA17','IAA1'),5),'IAA17'),8))
treatment <- rep(c(rep(0,11),rep(0.05,11),rep(1.5,11),rep(50,11)),4)
repl<-rep(c(rep(1,44),rep(2,44)),2)
mcols<-cbind(AFB,IAA,treatment,repl)

pData(dat)
pData(dat) <- cbind(pData(dat),mcols)
pData(dat)$strain <- paste(pData(dat)$AFB,pData(dat)$IAA,sep='x')
colnames(pData(dat))[-1]

dat.SS <- SteadyState(dat, ploidy = 'diploid', only = 'singlets')

theme_nem <- theme(
  panel.background = element_rect(fill="white"),
  #axis.ticks = element_line(colour=NA),
  panel.grid = element_line(colour="grey"),
  axis.text.y = element_text(colour='black'),
  axis.text.x = element_text(colour="black"),
  text = element_text(size=20, family="Helvetica")
)

theme_nem <-function (base_size = 16, base_family = "Arial")
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(legend.background = element_blank(), legend.key = element_blank(),
          panel.grid = element_line(colour="grey"),
          panel.background = element_blank(), panel.border = element_blank(),
          strip.background = element_blank(), plot.background = element_blank())
}

#Simple for presentations
p <- ggplot(subset(dat.SS, (treatment != 0.05 & treatment != 50 & IAA != 'IAA17')), aes(AFB, FL2.A/100, fill = AFB))
p <- p + geom_boxplot(outlier.size = 0) + scale_fill_manual(values = c(rep('navyblue',7),rep('dodgerblue1',4)))
p <- p + facet_grid(IAA~treatment)
p <- p + theme_nem(base_size = 20)
p <- p + ylim(c(-1,100)) + ylab('Fluorescence (AU)') + ggtitle('Auxin (uM)')
p <- p + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), legend.position = 'none')
p

#full set for paper
p <- ggplot(dat.SS, aes(AFB, FL2.A/100, fill = AFB))
p <- p + geom_boxplot(outlier.size = 0) + scale_fill_manual(values = c(rep('navyblue',7),rep('dodgerblue1',4)))
p <- p + facet_grid(IAA~treatment)
p <- p + theme_nem(base_size = 20)
p <- p + ylim(c(-1,100)) + ylab('Fluorescence (AU)') + ggtitle('Auxin (uM)')
p <- p + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
p

#playing around with colorspace for creating palettes
#<https://cran.r-project.org/web/packages/colorspace/vignettes/hcl-colors.pdf>
pal <- function(col, border = "light gray", ...)
{
    n <- length(col)
    plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
    axes = FALSE, xlab = "", ylab = "", ...)
    rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}
#we want to create a sequential or perhaps diverging palette
#could
#could order variants by level of substrate then set WT to hue and change chroma and luminescence for variants


p <- ggplot(dat.SS, aes(AFB, FL2.A/100, fill = AFB))
p <- p + geom_boxplot(outlier.size = 0) + scale_fill_manual(values = c(rep('navyblue',7),rep('dodgerblue1',4)))
p <- p + facet_grid(IAA~treatment)
p <- p + theme_nem(base_size = 20)
p <- p + ylim(c(-1,100)) + ylab('Fluorescence (AU)') + ggtitle('Auxin (uM)')
p <- p + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
p
