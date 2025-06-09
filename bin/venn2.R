#!/usr/bin/env Rscript
library(VennDiagram)

args <- commandArgs(trailing=T)
	site_f1 <- args[1]
	site_f2 <- args[2]
	pdf_f <- args[3]

A <- read.table(site_f1, header = F)
B <- read.table(site_f2, header = F)

a <- A[,1]
b <- B[,1]
ab <- intersect(a,b)
amb <- setdiff(a,b)
bma <- setdiff(b,a)

length_a <- length(a)
length_b <- length(b)
length_ab <- length(ab)

pdf(pdf_f, w=8, h=8)
par(mar=c(4.5, 4.5, 4.5, 4.5) + 0.5)
FillCol <- c("red", "blue")
draw.pairwise.venn(
	area1=length_a,area2=length_b,cross.area=length_ab,
	category=c(site_f1,site_f2),
	cat.cex=1,cat.pos=0,cat.fontface = "bold",
    #col=c('darkred','darkgreen'),
	col="black",lty=1,lwd=2,
    #rotation=1,
	alpha=0.50,cex=1.8,#????
	fontfamily="serif", rotation.degree = 0, 
    #rotation.centre = c(0.5, 0.5),
	sep.dist=0.05,
	fill=FillCol,cat.col=FillCol,
	scaled=F,reverse=F
)
