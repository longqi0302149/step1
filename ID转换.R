rm(list = ls())
load('basic_data.Rdata')
#GSE79962脓毒性心肌病和正常心肌，扩张性心肌病，缺血性心肌病的数据集
        #其中包括20名脓毒症患者和11例正常的donors心肌组织
        #
#GSE141864包含8例脓毒症患者的心肌组织和2例正常对照组的心肌组织
#GSE54514包含26名sepsis survivors and 9名 sepsis nonsurvivors
        #以及18名健康对照组患者
#GSE65682包含802名脓毒症患者和对照组的外周血测序和生存资料


#GSE79962和GSE54514的基因转换的数据集以前已经完成，直接复制
#到本文件夹即可，GSE141864和GSE65682的数据集需要在本次完成


#GSE65682的ID转换
if(F){
  library(GEOquery)
  #Download GPL file, put it in the current directory, and load it:
  gpl <- getGEO('GPL13667', destdir=".")
  colnames(Table(gpl))  
  head(Table(gpl)[,c(1,21)]) ## you need to check this , which column do you need
  probe2gene=Table(gpl)[,c(1,21)]
  head(probe2gene)
  library(stringr)  
  save(probe2gene,file='probe2gene.Rdata')
}
# 
# load(file='probe2gene.Rdata')
# ids=probe2gene 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("hgu219.db")
library(hgu219.db)
ids=toTable(hgu219SYMBOL) #toTable这个函数：通过看hgu133plus2.db这个包的说明书知道提取probe_id（探针名）和symbol（基因名）的对应关系的表达矩阵的函数为toTable
head(ids) #head为查看前六行

head(ids)
colnames(ids)=c('probe_id','symbol')  
ids=ids[ids$symbol != '',]
ids=ids[ids$probe_id %in%  rownames(dat_65682),]

dat_65682[1:4,1:4]   
dat=dat_65682[ids$probe_id,] 

ids$median=apply(dat,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids=ids[order(ids$symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids=ids[!duplicated(ids$symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
dat=dat[ids$probe_id,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
rownames(dat)=ids$symbol#把ids的symbol这一列中的每一行给dat作为dat的行名
dat[1:4,1:4]  #保留每个基因ID第一次出现的信息
gene_65682 <- dat
save(gene_65682,file = 'gene_65682.Rdata')




#GSE141864的ID转换
gse <- getGEO(filename = "GSE141864_family.soft.gz",
              destdir=".")       ## 平台文件
str(gse)
length(gse)



#
id_probe <- gse@gpls$GPL17586@dataTable@table

dim(id_probe)

head(id_probe)
id_probe[1:4,1:15]
View(head(id_probe))## you need to check this , which column do you need

probe2gene <- id_probe[,c(2,8)]
View(probe2gene)

##
library(stringr) 
probe2gene$symbol=trimws(str_split(probe2gene$gene_assignment,'//',simplify = T)[,2])
plot(table(table(probe2gene$symbol)),xlim=c(1,50))
head(probe2gene)

dim(probe2gene)
View(head(probe2gene))
ids2 <- probe2gene[,c(1,3)]
View(head(ids2))
ids2[1:20,1:2]#含有缺失值
table(table(unique(ids2$symbol)))#30907 ,30906个基因，一个空字符
save(ids2,probe2gene,file='gse-probe2gene.Rdata')


##
library(GEOmirror)

gset
a=exprs(gset[[1]])
a[1:4,1:4]
gset[[1]]@annotation

#过滤表达矩阵
exprSet <- a

library(dplyr)
exprSet <- exprSet[rownames(exprSet) %in% ids2$probeset_id,]
dim(exprSet)
exprSet[1:5,1:5]

#ids过滤探针
ids <- ids2[match(rownames(exprSet),ids2$probeset_id),]
dim(ids)
ids[1:5,1:2]
#ids2[1:5,1:2]

#合并表达矩阵和ids

idcombine <- function(exprSet, ids){
  tmp <- by(exprSet,
            ids$symbol,
            function(x) rownames(x)[which.max(rowMeans(x))])
  probes <- as.character(tmp)
  print(dim(exprSet))
  exprSet <- exprSet[rownames(exprSet) %in% probes,]
  
  print(dim(exprSet))
  rownames(exprSet) <- ids[match(rownames(exprSet), ids$probeset_id),2]
  return(exprSet)
}

new_exprSet <- idcombine(exprSet,ids)
new_exprSet[1:4,1:6]
dim(new_exprSet)

rownames(new_exprSet)
gene_141864 <- new_exprSet
save(gene_141864,file = "gene_141864.Rdata")

gene_79962 <- new_exprSet
save(gene_79962,file = "gene_79962.Rdata")

save(gene_141864, gene_54514, gene_65682, gene_79962,
     file = 'gene.Rdata')

##
##基因数据与临床资料的合并
#GSE79962
a1 <- gene_79962
a1 <- t(a1)
a2 <- pd_79962
a2 <- t(a2)
gse79962 <- cbind(a1, a2)

#GSE141864
a1 <- gene_141864
a2 <- pd_141864
gse141864 <- cbind(a1, a2)

#GSE54514
a1 <- gene_54514
a2 <- pd_54514
gse54514 <- cbind(a1, a2)

#GSE65682
a1 <- gene_65682
a2 <- pd_65682
gse65682 <- cbind(a1, a2)

#添加gse141864的分组信息
group_141864 <- matrix('sepsis',nrow = 10, ncol = 1)
group_141864[c(5,8)] = c('control')
a <- gene_141864[,1:10]
View(a)
a <- t(a)
a1 <- rbind(a,group)
View(a1)
geneGroup_141864 <- a1


#添加gse65682的分组信息
pd_65682$title
table(pd_65682$title)
View(pd_65682$title)
group_65682 <- matrix('sepsis',nrow = 802, ncol = 1)

group_65682[c(7,9,62,64,70:72,74:76,84,91:101,109:111,
        113:117,120,148,154,156,157,160:164,166,169)] = c('control')
table(group_65682)
a <- gene_65682
a <- t(a)
a1 <- cbind(a,group_65682)
View(a1)
geneGroup_65682 <- a1
geneGroup_65682[,10860]
geneGroup_65682 <- t(geneGroup_65682)


library(tidyverse)
rename(geneGroup_65682,group = V10860)
geneGroup_65682 <- data.frame(geneGroup_65682)

#添加gse79962的分组信息
table(pd_79962$title)
group_79962 <- matrix('sepsis',nrow = 51, ncol = 1)
group_79962[c(1:9)] = c('DM')
group_79962[c(10:20)] = c('IHD')
group_79962[c(21:31)] = c('NDH')
a <- gene_79962
head(a)
head(group_79962)
geneGroup_79962 <- a1
geneGroup_79962 <- t(geneGroup_79962)
group_79962_2<- matrix('sepsis',nrow = 51, ncol = 1)

group_79962_2[c(1:31)] = c('control')



#添加gse54514的分组信息
a <- pd_54514[,c(43,45,46)]
gene_54514[1:4,1:4]
a[1:4,]
gene_54514 <- t(gene_54514)
geneGroup_54514 <- cbind(gene_54514,a)
View(geneGroup_54514)
group_54514 <- matrix('control',nrow = 163, ncol = 1)
group_54514[c(19:145)] = c('sepsis')

save(geneGroup_141864, geneGroup_54514, geneGroup_65682, geneGroup_79962,
     file = 'geneGroup.Rdata')

save(group_141864, group_54514, group_65682, group_79962,group_79962_2,
     file = 'group.Rdata')


#心肌标本数据库的合并GSE79962和GSE141864
a1 <- rownames(gene_141864)
a2 <- rownames(gene_79962)
a3 <- intersect(a1,a2)
a4 <- gene_141864[colnames=a3,]
a4 <- a4[,1:10]
a5 <- gene_79962[colnames=a3,]
a6 <- cbind(a4,a5)
colnames(a4)
colnames(a5)
colnames(a6)
gene_heart <- a6
group_heart_1 <- rbind(group_141864,group_79962)
group_heart_2 <- rbind(group_141864,group_79962_2)


save(geneGroup_141864, geneGroup_54514, geneGroup_65682, geneGroup_79962,gene_heart,
     file = 'geneGroup.Rdata')

save(group_141864, group_54514, group_65682, group_79962,group_79962_2,group_heart_1,group_heart_2,
     file = 'group.Rdata')









