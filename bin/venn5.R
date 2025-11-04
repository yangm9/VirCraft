#!/usr/bin/env Rscript

if(!require("VennDiagram", quietly = TRUE)){
    install.packages("VennDiagram", repos="https://cloud.r-project.org/")
}

library(VennDiagram)

args <- commandArgs(trailing=T)
	site_f1 <- args[1]
	site_f2 <- args[2]
	site_f3 <- args[3]
	site_f4 <- args[4]
	site_f5 <- args[5]
	pdf_f <- args[6]

name_f1 <- sub("\\.[^.]*$", "", basename(site_f1))
name_f2 <- sub("\\.[^.]*$", "", basename(site_f2))
name_f3 <- sub("\\.[^.]*$", "", basename(site_f3))
name_f4 <- sub("\\.[^.]*$", "", basename(site_f4))
name_f5 <- sub("\\.[^.]*$", "", basename(site_f5))

A <- read.table(site_f1, header = F)
B <- read.table(site_f2, header = F)
C <- read.table(site_f3, header = F)
D <- read.table(site_f4, header = F)
E <- read.table(site_f5, header = F)

a <- A[,1]
b <- B[,1]
c <- C[,1]
d <- D[,1]
e <- E[,1]

length_a <- length(a)
length_b <- length(b)
length_c <- length(c)
length_d <- length(d)
length_e <- length(e)
length_ab <- length(intersect(a, b))
length_ac <- length(intersect(a, c))
length_ad <- length(intersect(a, d))
length_ae <- length(intersect(a, e))
length_bc <- length(intersect(b, c))
length_bd <- length(intersect(b, d))
length_be <- length(intersect(b, e))
length_cd <- length(intersect(c, d))
length_ce <- length(intersect(c, e))
length_de <- length(intersect(d, e))
length_abc <- length(intersect(intersect(a, b), c))
length_abd <- length(intersect(intersect(a, b), d))
length_abe <- length(intersect(intersect(a, b), e))
length_acd <- length(intersect(intersect(a, c), d))
length_ace <- length(intersect(intersect(a, c), e))
length_ade <- length(intersect(intersect(a, d), e))
length_bcd <- length(intersect(intersect(b, c), d))
length_bce <- length(intersect(intersect(b, c), e))
length_bde <- length(intersect(intersect(b, d), e))
length_cde <- length(intersect(intersect(c, d), e))
length_abcd <- length(intersect(intersect(intersect(a, b), c), d))
length_abce <- length(intersect(intersect(intersect(a, b), c), e))
length_abde <- length(intersect(intersect(intersect(a, b), d), e))
length_acde <- length(intersect(intersect(intersect(a, c), d), e))
length_bcde <- length(intersect(intersect(intersect(b, c), d), e))
length_abcde <- length(intersect(intersect(intersect(intersect(a, b), c), d), e))

pdf(pdf_f, width=10, height=10)
FillCol <- c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")

draw.quintuple.venn(
    area1=length_a, area2=length_b, area3=length_c, area4=length_d, area5=length_e, 
    n12=length_ab, n13=length_ac, n14=length_ad, n15=length_ae, n23=length_bc, 
    n24=length_bd, n25=length_be, n34=length_cd, n35=length_ce, n45=length_de, 
    n123=length_abc, n124=length_abd, n125=length_abe, n134=length_acd, n135=length_ace, 
    n145=length_ade, n234=length_bcd, n235=length_bce, n245=length_bde, n345=length_cde, 
    n1234=length_abcd, n1235=length_abce, n1245=length_abde, n1345=length_acde, n2345=length_bcde, 
    n12345=length_abcde,
    category=c(name_f1, name_f2, name_f3, name_f4, name_f5),
    cat.cex=1.8,
    cat.fontface = "bold",
    cat.col=FillCol,
    #col=c('darkred','darkgreen','darkblue'),
    lty=1,
    lwd=2,
    alpha=0.50,
    cex=2,
    #rotation=1,
    fontfamily="serif", 
    rotation.degree=60, 
    #rotation.centre = c(0.5, 0.5),
    fill=FillCol,
    reverse=F
)
