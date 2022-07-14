if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

rm(list = ls())
options(stringsAsFactors = F)#在调用as.data.frame的时，将stringsAsFactors设置为FALSE可以避免character类型自动转化为factor类型
# 注意查看下载文件的大小，检查数据 
f='GSE54514_eSet.Rdata'
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11121
library(GEOquery)
# 这个包需要注意两个配置，一般来说自动化的配置是足够的。
#Setting options('download.file.method.GEOquery'='auto')
#Setting options('GEOquery.inmemory.gpl'=FALSE)
if(!file.exists(f)){
  gset <- getGEO('GSE54514', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE54514_eSet.Rdata')  ## 载入数据
class(gset)  #查看数据类型
length(gset)  #
class(gset[[1]])

a=gset[[1]] #
dat_54514=exprs(a) #a现在是一个对象，取a这个对象通过看说明书知道要用exprs这个函数
dim(dat)#看一下dat这个矩阵的维度

pd_54514=pData(a)

gset_141864 <- gset
pd_65682 <- pd

dat_141864 <- dat


#这一部分下载了四个数据库的数据，提取了他们的测序数据和
#临床数据，以供下一步分析使用
save(dat_65682, dat_141864, pd_65682, pd_141864,
     dat_79962, pd_79962, dat_54514, 
     pd_54514, file='basic_data.Rdata')





#热图
library(ggplot2)
library(pheatmap)
heatmap_heart <- heart[c('SMU1','NNMT','TMEM184C','TBC1D22B','DECR1','VPS45'),]
heatmap_heart <- as.matrix(scale(heatmap_heart))
heatmap_blood <- as.matrix(scale(heatmap_blood))
heatmap_heart <- t(heatmap_heart)
heatmap_heart$name <- rownames(heatmap_heart)
p <- ggplot(heatmap_heart, aes())
p + geom_tile()


col <- colorRampPalette(c("red", "blue"))(256)
pheatmap(heatmap_heart, cluster_rows= F, cluster_cols=F)
heatmap_blood <- gse65682[,c('SMU1','NNMT','TMEM184C','TBC1D22B','DECR1','VPS45')]
heatmap_blood <- as.matrix((heatmap_blood))
pheatmap(heatmap_heart)

col <- colorRampPalette(c("red", "white", "blue"))(256)
heatmap(heatmap_heart,scale = "none", col=col)


#箱式图
box_heart <- heart[c('SMU1','NNMT','TMEM184C','TBC1D22B','DECR1','VPS45'),]
box_heart <- t(box_heart)
box_heart<- as.data.frame(scale(box_heart))
box_heart$group <- group_list
ggplot(box_blood, aes(x=factor(group), y =NNMT)) +geom_boxplot()
box_blood <- gse65682[,c('SMU1','NNMT','TMEM184C','TBC1D22B','DECR1','VPS45')]
box_blood <- as.data.frame(scale(box_blood))
box_blood$group <- group_65682

View(d1)
d2 <- d1[d1$logFC >0,]
d3 <- d1[d1$logFC <0,]
d4 <- DEGs_gse65682[DEGs_gse65682$logFC >0,]
d5 <- DEGs_gse65682[DEGs_gse65682$logFC <0,]
d6 <- DEGs_gse65682[union(intersect(rownames(d2),rownames(d4)), intersect(rownames(d3),rownames(d5))),]
union(intersect(rownames(d2),rownames(d4)), intersect(rownames(d3),rownames(d5)))
write.csv(group_54514, file = "D:/数据库/Sepsis/GSE54514/group_54514.csv")         

#D:/数据库/Sepsis/GSE65682


library(AnnoProbe)
colnames(gse79962)
a0 <- annoGene(colnames(gse79962), 'SYMBOL', species = "human")
table(a0$biotypes)
