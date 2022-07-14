rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
# 每次都要检测数据
geneGroup_141864[1:4,1:4] 
table(group_heart_1) #table函数，查看group_list中的分组个数
#通过为每个数据集绘制箱形图，比较数据集中的数据分布
boxplot(gene_heart[1,]) #按照group_list分组画箱线图

#了解批次间差异
library(ggplot2)
library(reshape2)
exprSet_L = melt(combat_edata )
colnames(exprSet_L) = c('probe', 'sample', 'value')
exprSet_L$group = rep(group_heart_1, each=nrow(gene_heart))
head(exprSet_L)
p = ggplot(exprSet_L,aes(x=sample, y=value, fill = group))+
  geom_boxplot()


bp=function(g){         #定义一个函数g，函数为{}里的内容
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
bp(dat[1,]) ## 调用上面定义好的函数，避免同样的绘图代码重复多次敲。
bp(dat[2,])
bp(dat[3,])
bp(dat[4,])
dim(dat)


#消除批次间差异

table(group_heart_1)
batch <- matrix('2',nrow = 61, ncol = 1)
batch[c(1:10)] = c('1')
batch <- batch[,1]
a <- cbind(group_heart_1,batch)
a <- data.frame(a)
a1 <- colnames(gene_heart) 
a <- cbind(a,a1)


gene_heart <- data.frame(gene_heart)

modcombat = model.matrix(~1, data=a)

combat_edata = ComBat(dat=gene_heart, batch=batch, mod=modcombat, par.prior=TRUE, prior.plot=FALSE)
heart <- gene_heart2 
batchinf <- a
save(heart, batchinf, group_heart_1, group_heart_2,file = 'Heart.Rdata')



