#!/usr/bin/env Rscript

if(!require("VennDiagram", quietly = TRUE)){
    install.packages("VennDiagram", repos="https://cloud.r-project.org/")
}

library(VennDiagram)

args <- commandArgs(trailing=T)
	site_f1 <- args[1]
	site_f2 <- args[2]
	pdf_f <- args[3]

name_f1 <- sub("\\.[^.]*$", "", basename(site_f1))
name_f2 <- sub("\\.[^.]*$", "", basename(site_f2))

A <- read.table(site_f1, header = F)
B <- read.table(site_f2, header = F)

a <- A[,1]
b <- B[,1]

length_a <- length(a)
length_b <- length(b)
length_ab <- length(intersect(a, b))

pdf(pdf_f, width=10, height=10)
par(mar=c(4.5, 4.5, 4.5, 4.5) + 0.5)
FillCol <- c("red", "blue")
draw.pairwise.venn(
	area1=length_a, area2=length_b, cross.area=length_ab,
	category=c(name_f1, name_f1),
	cat.cex=1.8,
    cat.pos=0,
    cat.fontface = "bold",
    cat.col=FillCol,
    #col=c('darkred','darkgreen'),
	col="black", 
    lty=1,
    lwd=2,
    #rotation=1,
	alpha=0.50,
    cex=2,
	fontfamily="serif",
    rotation.degree = 0, 
    #rotation.centre = c(0.5, 0.5),
	sep.dist=0.05,
	fill=FillCol,
	scaled=F,
    reverse=F
)
