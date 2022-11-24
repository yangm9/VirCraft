setwd('G:\\project\\YM1XLD211101MV-西沙永乐龙洞生态系统结构和生物无氧适应机制\\virus\\5.abundance\\diversity')
library(vegan)
otu <- read.table(file="all_merged3.tpm.xls",sep="\t",header=T,check.names=FALSE,row.names=1)
otu_t <- t(otu)
Shannon=diversity(otu_t,"shannon") 
Simpson=diversity(otu_t,"simpson") 
InSimpson=diversity(otu_t,"inv")
specN=specnumber(otu_t)
Pielou=Shannon/log(specN)  #计算Pielou均匀度指数
alpha_df=t(rbind(specN,Pielou,Shannon,Simpson,InSimpson))
alpha_df.reset_index()
write.table(alpha_df,'alpha_diversity.xls',sep='\t',row.names = TRUE)

group_df <- read.table("groups.txt", sep='\t', header=T)
colnames(group_df) <- c("samples","group")
names(otu_t)[1]<-"samples"
df <- merge(otu_t,group_df,by="samples")
S_t <- t(t(S))
