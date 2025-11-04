#!/usr/bin/env Rscript

if(!require("VennDiagram", quietly = TRUE)){
    install.packages("VennDiagram", repos="https://cloud.r-project.org/")
}

library(VennDiagram)

args <- commandArgs(trailing=T)
	site_f1 <- args[1]
	site_f2 <- args[2]
	site_f3 <- args[3]
	pdf_f <- args[4]

name_f1 <- sub("\\.[^.]*$", "", basename(site_f1))
name_f2 <- sub("\\.[^.]*$", "", basename(site_f2))
name_f3 <- sub("\\.[^.]*$", "", basename(site_f3))

A <- read.table(site_f1, header = F)
B <- read.table(site_f2, header = F)
C <- read.table(site_f3, header = F)

a <- A[,1]
b <- B[,1]
c <- C[,1]

length_a <- length(a)
length_b <- length(b)
length_c <- length(c)
length_ab <- length(intersect(a, b))
length_ac <- length(intersect(a, c))
length_bc <- length(intersect(b, c))
length_abc <- length(intersect(intersect(a, b), c))

pdf(pdf_f, w=8, h=8)
par(mar=c(4.5, 4.5, 4.5, 4.5) + 0.5)
FillCol <- c("dodgerblue", "goldenrod1", "darkorange1")
draw.triple.venn(
    area1=length_a, area2=length_b, area3=length_c,
    n12=length_ab, n23=length_bc, n13=length_ac, n123=length_abc,
    category=c(name_f1, name_f2, name_f3),
    cat.cex=1.8,
    cat.fontface="bold",
    cat.col=FillCol,
    col="black",
    lty=1,lwd=2,
    #rotation=1,
    alpha=0.50,
    cex=1.8,
    fontfamily="serif", 
    rotation.degree=60, 
    #rotation.centre = c(0.5, 0.5),
    fill=FillCol,
    reverse=F
)
