setwd("G:\\project\\YM1XLD211101MV-西沙永乐龙洞生态系统结构和生物无氧适应机制\\virus\\7.gene_annotation")
library(ggplot2)
library(grid)
library(RColorBrewer)
#color = rainbow(24)
d=read.table('kegg_lv2_genenum.xls', sep='	')
ggplot(d,aes(x=V1,y=V2,fill=V1))+
  geom_bar(stat="identity") + coord_flip()+
  labs(title='KEGG Level2 Classification (all)', x='', y='Number of Genes')+
  theme(panel.background = element_rect(fill='transparent'),
        panel.grid=element_line(color='grey'),
        panel.border=element_rect(fill='transparent',color='black'),
        legend.position="none", axis.text=element_text(color='black', size=12))+
  theme(legend.title=element_blank(), legend.text=element_text(angle=270),
        legend.key.width=unit(1, 'mm'), legend.key.height=unit(4, 'cm'),  
        legend.text.align=0.7, plot.title=element_text(face='bold'))+
  geom_text(aes(label=V2), hjust=-0.5, vjust=0.5, size = 4)+
  scale_y_continuous(limits=c(0,max(d$V2)*1.2),trans='sqrt')