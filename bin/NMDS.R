if(!require("vegan")){install.packages("vegan")}
if(!require("ggplot2")){install.packages("ggplot2")}

library(vegan)#计算距离时需要的包
library(ggplot2)#绘图包

argv<-commandArgs(T)

if(length(argv)<2){
    stop("inputs: <otus_tab.xls> <group.txt> <outdir>")
}

otu <- read.table(argv[1],sep="\t",header=T,check.names=FALSE,row.names=1)
otu_t <- t(otu)
otu.distance <- vegdist(otu_t,method='bray')
df_nmds <- metaMDS(otu.distance,k=2)
summary(df_nmds)
df_nmds_stress <- df_nmds$stress
pdf(paste(argv[3],'/stressplot.pdf',sep=''),width=10,height=8)
stressplot(df_nmds)
dev.off()

#提取作图数据
df_points <- as.data.frame(df_nmds$points)
#添加samp1es变量
df_points$sample <- row.names(df_points)
names(df_points)[1:2] <- c('NMDS1','NMDS2')

#添加分组
group <- read.table(argv[2],sep='\t',header=T)
group <- group[,0:2]
colnames(group) <- c("sample","group")
df <- merge(df_points,group,by="sample")
color <- c("#1597A5","#FFC24B","#FEB3AE","#FEB3AE")#颜色变量

plt <- ggplot(data=df,aes(x=NMDS1,y=NMDS2))+#指定数据、X轴、Y轴，颜色
    theme_bw()+#主题设置
    geom_point(aes(color=group),shape=19,size=3)+#绘制点图并设定大小
    theme(panel.grid=element_blank())+
    geom_vline(xintercept=0,lty="dashed",size=1,color='grey50')+
    geom_hline(yintercept=0,lty="dashed",size=1,color='grey50')+#图中虚线
    geom_text(aes(label=sample,
                y=NMDS2+0.03,x=NMDS1+0.03,
                vjust=0,color=group),
              size=3.5,show.legend=F)+#添加数据点的标签
    stat_ellipse(data=df, #添加椭圆
                 geom="polygon",level=0.95,
                 linetype=2,size=0.5,
                 aes(fill=group),
                 alpha=0.2)+
    scale_color_manual(values=color) +#点的颜色设置
    scale_fill_manual(values=color)+#椭圆颜色
    theme(axis.title.x=element_text(size=12),#修改X轴标题文本
        axis.title.y=element_text(size=12,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=10),#修改x轴刻度标签文本
        axis.text.x=element_text(size=10),#修改y轴刻度标签文本
        panel.grid=element_blank())+#隐藏网格线
    ggtitle(paste('Stress=',round(df_nmds_stress, 3)))#添加应力函数值

pdf(paste(argv[3],'/NMDS.pdf',sep=''),width=10,height=8)
plt
dev.off()
