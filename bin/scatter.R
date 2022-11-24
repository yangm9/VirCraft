if (!require("ggplot2")) {install.packages("ggplot2")}

library(ggplot2)

argv<-commandArgs(T)

#<Contig_len_sum_tpm_quality.xls>

df=read.table(
    argv[1],header=TRUE,sep="\t",
    row.names=1,check.names=F,quote=""
)
plot<-ggplot(df,
       aes(x=log(Contig_Length),
           y=log(Total_Abundance),
           color=checkv_quality,
           shape=Source
          ))+ 
  geom_point()+
  ggtitle("The size of viral contigs (x axis) and their total abundance")+
  scale_shape_manual(values=c(1,3,5))+theme_classic()
pdf(paste(argv[2],"/length_abundance_scatter.pdf",width=10,height=8)
plot
  #scale_color_manual(values=c('#E00303','#03E0B6','#E0BF03',"#EE9A00","#006400"))+
  #geom_point(
    #shape=21,
    #size=4,
    #stroke =1.5
  #)
