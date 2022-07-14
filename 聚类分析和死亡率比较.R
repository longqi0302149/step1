genes <- read.csv('mortality_genes_clustering.csv')
genes <- genes[,2:213]
rownames(genes) <- geneNames$X
geneNames <- read.csv('mortality_group.csv')

#聚类分析
library(NbClust)
df <- scale(genes)
nc <- NbClust(df, min.nc = 2, max.nc = 15, method = "kmeans")
table(nc$Best.nc[1,])
table(nc$Best.partition)
group <- nc$Best.partition
a <- cbind(group, geneNames$mortalty)
table(a)
a1 <- a[a$group==1,]
table(a1$V2)
a2 <- a[a$group==2,]
table(a2$V2)
a3 <- a[a$group==3,]
table(a3$V2)
clustering <- write.csv(group, file = "clustering.csv")
barplot(table(nc$Best.nc[1,]))


#死亡率比较
library(survival)
library(ISwR)
library(ggplot2)
library(dplyr)
library(survminer)
dataFormortality <- as.data.frame(cbind(group, geneNames$mortalty, geneNames$days))
colnames(dataFormortality) <- c("group", "mortality", "days")
attach(dataFormortality)
Surv(days, mortality==1)
fit <- survfit(Surv(days, mortality)~group, data = dataFormortality)
ggsurvplot(fit, pval = TRUE,
           conf.int = T,
           conf.int.style="ribbon",
           conf.int.alpha=0.1)#显示置信区间


ggsurvplot(fit,risk.table=TRUE,#生存统计统计表
           
          # conf.int=TRUE,#添加置信区间带
           
           palette = c("skyblue","red","yellow"),#颜色设置
           
           pval=TRUE,#log-rank检验
           
           pval.method=TRUE)#添加检验text


