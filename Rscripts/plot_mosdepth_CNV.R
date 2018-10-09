library(ggplot2)
library(RColorBrewer)
colors1 <- colorRampPalette(brewer.pal(8, "RdYlBu"))
manualColors = c("dodgerblue2","red1","grey20")

bedwindows = read.table("coverage/mosdepth.10000bp.gg.tab.gz",header=F)
colnames(bedwindows) = c("Chr","Start","End","Depth","Group","Strain")
#bedwindows = subset(bedwindows,bedwindows$Chr != "MT_CBS_6936") # drop MT for this

bedwindows$CHR <- sub("scaffold_","",bedwindows$Chr,perl=TRUE)
chrlist = 1:8
d=bedwindows[bedwindows$CHR %in% chrlist, ]

d <- d[order(d$CHR, d$Start), ]
d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$Start,d$CHR,length)) 

d$pos=NA

nchr = length(unique(chrlist))
lastbase=0
ticks = NULL
minor = vector(,8)

for (i in 1:nchr ) {
    if (i ==1) {
        d[d$index==i, ]$pos = d[d$index==i, ]$Start
    } else {
        ## chromosome position maybe not start at 1, eg. 9999. So gaps may be produced. 
        lastbase = lastbase + max(d[d$index==(i-1),"Start"])
	      minor[i] = lastbase
        d[d$index == i,"Start"] =
             d[d$index == i,"Start"]-min(d[d$index==i,"Start"]) +1
        d[d$index == i,"End"] = lastbase
        d[d$index == i, "pos"] = d[d$index == i,"Start"] + lastbase
    }
}
ticks <-tapply(d$pos,d$index,quantile,probs=0.5)
ticks
minorB <- tapply(d$End,d$index,max,probs=0.5)
minorB
minor
d$Group = factor(d$Group, levels = c("LL", "UL", 
                                     "UUMRR","Sp1", "Sp2", "B", "C", "CL","PH"))
xmax = ceiling(max(d$pos) * 1.03)
xmin = floor(max(d$pos) * -0.03)

pdffile="plots/Genomewide_cov_by_10kb_win_mosdepth.pdf"
pdf(pdffile,width=7,height=2.5)
Title="Depth of sequence coverage"

#What about the color scheme I have for Ul/LL/Sp in Fig 1 which is Upper=bright blue, lower=red, sputum=black/dark gray


p <- ggplot(d,
            aes(x=pos,y=Depth,color=Group)) +
	        geom_vline(mapping=NULL, xintercept=minorB,alpha=0.5,size=0.1,colour='grey15')	+
    geom_point(alpha=0.8,size=0.4,shape=16) +
    scale_color_brewer(palette="RdYlBu",type="seq") +
    labs(title=Title,xlab="Position",y="Normalized Read Depth") +
    scale_x_continuous(name="Chromosome", expand = c(0, 0),
                       breaks = ticks,                      
                       labels=(unique(d$CHR))) +
    scale_y_continuous(name="Normalized Read Depth", expand = c(0, 0),
                       limits = c(0,3)) + theme_classic() + 
    guides(fill = guide_legend(keywidth = 3, keyheight = 1)) 

p


bedwindows = read.table("coverage/mosdepth.5000bp.gg.tab.gz",header=F)
colnames(bedwindows) = c("Chr","Start","End","Depth","Group","Strain")
#bedwindows = subset(bedwindows,bedwindows$Chr != "MT_CBS_6936") # drop MT for this

bedwindows$CHR <- sub("scaffold_","",bedwindows$Chr,perl=TRUE)
# reuse chrlist
d=bedwindows[bedwindows$CHR %in% chrlist, ]

d <- d[order(d$CHR, d$Start), ]
d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$Start,d$CHR,length)) 

d$pos=NA

#reuse from before
#nchr = length(unique(d$CHR))
lastbase=0
ticks = NULL
minor = vector(,8)
for (i in 1:nchr ) {
    if (i==1) {
        d[d$index==i, ]$pos=d[d$index==i, ]$Start
    } else {
        ## chromosome position maybe not start at 1, eg. 9999. So gaps may be produced. 
        lastbase = lastbase + max(d[d$index==(i-1),"Start"])
      	minor[i] = lastbase
        d[d$index == i,"Start"] =
             d[d$index == i,"Start"]-min(d[d$index==i,"Start"]) +1
        d[d$index == i,"End"] = lastbase
        d[d$index == i, "pos"] = d[d$index == i,"Start"] + lastbase
    }
}
ticks <-tapply(d$pos,d$index,quantile,probs=0.5)
ticks
minorB <- tapply(d$End,d$index,max,probs=0.5)
minorB
minor
d$Group = factor(d$Group, levels = c("LL", "UL", 
                                     "UUMRR","Sp1", "Sp2", "B", "C", "CL","PH"))

xmax = ceiling(max(d$pos) * 1.03)
xmin = floor(max(d$pos) * -0.03)

plot_strain <- function(strain,data) {
 l = subset(data,data$Strain == strain)
 Title=sprintf("Chr coverage plot for %s",strain)
 p <- ggplot(l,
            aes(x=pos,y=Depth,color=CHR))  + 
    scale_colour_brewer(palette = "Set3") +
    geom_point(alpha=0.9,size=0.8,shape=16) +
    labs(title=Title,xlab="Position",y="Normalized Read Depth") +
    scale_x_continuous(name="Chromosome", expand = c(0, 0),
                       breaks=ticks,
                       labels=(unique(d$CHR))) +
    scale_y_continuous(name="Normalized Read Depth", expand = c(0, 0),
                       limits = c(0,3)) + theme_classic() +
    guides(fill = guide_legend(keywidth = 3, keyheight = 1))
}

plts <- lapply(unique(d$Strain),plot_strain,data=d)
pdf("plots/StrainPlot_5kb.pdf")
plts

plot_chrs <-function (chrom, data) {
  Title=sprintf("Chr%s depth of coverage",chrom)
  l <- subset(data,data$CHR==chrom)
  l$bp <- l$Start
  p<-ggplot(l,
            aes(x=bp,y=Depth,color=Strain)) +
    geom_point(alpha=0.7,size=0.8,shape=16) +
    # scale_color_brewer(palette="RdYlBu",type="seq") +
    labs(title=Title,xlab="Position",y="Normalized Read Depth") +
    scale_x_continuous(expand = c(0, 0), name="Position") +
    scale_y_continuous(name="Normalized Read Depth", expand = c(0, 0),
                       limits = c(0,3)) + theme_classic() +
    guides(fill = guide_legend(keywidth = 3, keyheight = 1))
}
pdf("plots/ChrPlot_5kb.pdf")
plts <- lapply(1:nchr,plot_chrs,data=d)
plts

