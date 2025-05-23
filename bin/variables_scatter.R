#!/usr/bin/env Rscript
#yangm@idsse.ac.cn
if(!require("ggplot2")) {install.packages("ggplot2")}

library(ggplot2)

argv<-commandArgs(T)

if(length(argv)<3){stop("inputs: <variables.tsv> <Variables> <plot_dir>")}

df=read.table(
    argv[1],header=T,sep="\t",
    row.names=1,check.names=F,quote=""
)

Variables=strsplit(argv[2],"~")[[1]]
X_Val=log(df[,Variables[1]])/log(10)
Y_Val=df[,Variables[2]]
n=0
if(is.na(Variables[3])){
    Variables[3]="Colour"
    df[,Variables[3]]="NA"
    n=n+1
}
if(is.na(Variables[4])){
    Variables[4]="Shape"
    df[,Variables[4]]="NA"
    n=n+1
}
Title=paste("The_",Variables[1],"_and_",Variables[2],sep='')
plot<-ggplot(df,
        aes(x=X_Val,y=Y_Val,
            color=df[,Variables[3]],
            shape=df[,Variables[4]]
        )
    )+
    labs(x=paste("log10(",Variables[1],")",sep=""),y=Variables[2],
        color=Variables[3],shape=Variables[4]
    )+ #坐标轴和图例命名
    geom_point()+
    ggtitle(Title)+
    scale_shape_manual(values=c(1,3,5))+theme_classic()

out_pdf=paste(argv[3],"/",Variables[1],"_",Variables[2],"_scatter.pdf",sep="")
pdf(out_pdf,width=10,height=8)
if(n==2){
    plot+theme(legend.position='none')
}else{
    plot
}
dev.off()
