#!/usr/bin/env Rscript

if(!require("ggplot2")) {install.packages("ggplot2")}

library(ggplot2)

argv<-commandArgs(T)

if(length(argv)<3){stop("inputs: <variables.xls> <Variables> <plot_dir>")}

df=read.table(
    argv[1],header=T,sep="\t",
    row.names=1,check.names=F,quote=""
)

Variables=strsplit(argv[2],"~")[[1]]
X_Val=log(df[,Variables[1]])
Y_Val=log(df[,Variables[2]])
if(!exists(Variables[3])){i
    Variables[3]="Colour"
    df[,Variables[3]]="NA"
}
if(!exists(Variables[4])){
    Variables[4]="Shape"
    df[,Variables[4]]="NA"
}
Title=paste("The_",Variables[1],"_and_",Variables[2],sep='')
plot<-ggplot(df,
        aes(x=X_Val,y=Y_Val,
            color=df[,Variables[3]],
            shape=df[,Variables[4]]
        )
    )+
    labs(x=Variables[1],y=Variables[2],
        color=Variables[3],shape=Variables[4]
    )+ #坐标轴和图例命名
    geom_point()+
    ggtitle(Title)+
    scale_shape_manual(values=c(1,3,5))+theme_classic()

out_pdf=paste(argv[3],"/",Variables[1],"_",Variables[2],"_scatter.pdf",sep="")
pdf(out_pdf,width=10,height=8)
plot
  #scale_color_manual(values=c('#E00303','#03E0B6','#E0BF03',"#EE9A00","#006400"))+
  #geom_point(
    #shape=21,
    #size=4,
    #stroke =1.5
  #)
