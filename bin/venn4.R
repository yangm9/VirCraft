library(grid)
library(VennDiagram)
args <- commandArgs(trailing=T)
	site_f1 <- args[1]
	site_f2 <- args[2]
	site_f3 <- args[3]
	site_f4 <- args[4]
	pdf_f <- args[5]


A=read.table(site_f1,header = F)
B=read.table(site_f2,header = F)
C=read.table(site_f3,header = F)
D=read.table(site_f4,header = F)

a=A[,1]
b=B[,1]
c=C[,1]
d=D[,1]

length_a=length(a)
length_b=length(b)
length_c=length(c)
length_d=length(d)
length_ab=length(intersect(a,b))
length_ac=length(intersect(a,c))
length_ad=length(intersect(a,d))
length_bc=length(intersect(b,c))
length_bd=length(intersect(b,d))
length_cd=length(intersect(c,d))
length_abc=length(intersect(intersect(a,b),c))
length_abd=length(intersect(intersect(a,b),d))
length_acd=length(intersect(intersect(a,c),d))
length_bcd=length(intersect(intersect(b,c),d))
length_abcd=length(intersect(intersect(intersect(a,b),c),d))

pdf(pdf_f)
draw.quad.venn(
    area1=length_a, area2=length_b, area3=length_c, area4=length_d,
    n12=length_ab, n13=length_ac, n14=length_ad, n23=length_bc, n24=length_bd, n34=length_cd,
    n123=length_abc, n124=length_abd, n134=length_acd, n234=length_bcd, n1234=length_abcd, 
    #category = c("A", "B", "C", "D"),
    category=c(site_f1,site_f2,site_f3,site_f4),
    cat.cex=1,
    cat.fontface = "bold",
    #col=c('darkred','darkgreen','darkblue'),
    lty=1,
    lwd=2,
    #rotation=1,
    alpha=0.50,#????
    cex=1.8,#????
    fontfamily = "serif", 
    rotation.degree = 0, 
    #rotation.centre = c(0.5, 0.5),
    fill=c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
    cat.col=c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
    reverse = F
)
