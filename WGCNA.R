rm(list = ls())
colnames(batchinf)
hclust_table <- batchinf[,1:2]
rownames(hclust_table) = colnames(heart)
colnames(hclust_table) = c('condition', 'group')


exprSet <- heart
exprSet[1:5,1:5]


group_list = as.character( batchinf[,1])
group_list <- group_list[11:61]
table(group_list)

group_list_2= as.character( batchinf[,2])
table(group_list_2)

group_list[c(5,8)] <- c('sepsis')


#心肌标本的DEGs分析
#GSE79962的差异基因分析
library( "limma" )
{
  design <- model.matrix( ~0 + factor( group_list ) )
  colnames( design ) = levels( factor( group_list  ) )
  rownames( design ) = colnames( heart )
}
design

contrast.matrix <- makeContrasts( "NDH-sepsis", levels = design )
contrast.matrix

## design和contrast.matrix都是为了下面的差异分析制作的

{
  fit <- lmFit( heart, design )
  fit2 <- contrasts.fit( fit, contrast.matrix ) 
  fit2 <- eBayes( fit2 )
  nrDEG_NDH2sepsis = topTable( fit2, coef = 1, n = Inf ) 
  write.table( nrDEG_NDH2sepsis, file = "nrDEG.out")
}
nrDEG_NDH2sepsis[rownames(hubgene_1),]

d1 <- nrDEG_NDH2sepsis[nrDEG_NDH2sepsis$P.Value < 0.05,]
d2 <- nrDEG_NDH2DM[nrDEG_NDH2DM$P.Value < 0.05,]
d3 <- nrDEG_NDH2IHD[nrDEG_NDH2IHD$P.Value < 0.05,]
rownames(d1)
d4 <- intersect(rownames(d1),rownames(d2))
d5 <- intersect(d4, rownames(d3))

head(nrDEG_65682)
nrDEG[nrDEG$P.Value < 0.01,]
DEGs_blood <- nrDEG[nrDEG$P.Value < 0.01,]
a1 <- rownames(DEGs_blood)
d6 <- DEGs_gse65682[DEGs_gse65682$P.Value < 0.05,]
dim(d1[d1$logFC < 0,] )
intersect(rownames(d1[d1$logFC < 0,]), rownames(d6[d6$logFC < 0,]))



# 计算数据集的每个基因的平均绝对离差，取排在前面最大的5000个基因。
datExpr = t(exprSet[order(apply(exprSet,1,mad), decreasing = T)[1:5000],])
# 取sepsis和control的差异表达基因进行WGCNA
x1 <- rownames(DEGs_gse65682[DEGs_gse65682$adj.P.Val < 0.05,])
x2 <- rownames(DEGs_gse65682[DEGs_gse65682$adj.P.Val < 0.05 & DEGs_gse65682$logFC < 0,])
x3 <- rownames(DEGs_heart[DEGs_heart$adj.P.Val < 0.05,])
x4 <- rownames(DEGs_heart[DEGs_heart$adj.P.Val < 0.05 & DEGs_heart$logFC < 0,])
x5 <-DEGs_gse65682[intersect(x1,x3),]
x5$group <- c('up')
x6 <-DEGs_gse65682[intersect(x2,x4),]
x6$group <- c('down')
x7 <- rbind(x5,x6)
x8 <- cbind(rownames(DEGs_heart), moduleColors)
x8 <- as.data.frame(x8)
rownames(x8) <- rownames(DEGs_heart)
x8 <- x8[,2]  
x9 <- x8[rownames(x7),]
x10 <- x8[rownames(x5),]
write.csv(x8, file = "C:/Users/Administrator/Desktop/论文/脓毒血症，脓毒性心肌病，建模/supplement-返修/ModuleColorinHeart.csv")
write.csv(x9, file = "C:/Users/Administrator/Desktop/论文/脓毒血症，脓毒性心肌病，建模/supplement-返修/ModuleColorinBlood.csv")         
write.csv(x10, file = "C:/Users/Administrator/Desktop/论文/脓毒血症，脓毒性心肌病，建模/supplement-返修/ModuleColorinBlood_all.csv")         

          
datExpr <- exprSet[rownames(d1),]
datExpr <- heart[rownames(DEGs_heart),]
dim(datExpr)
datExpr <- t(datExpr)
#第二步，绘制进化树
install.packages('WGCNA')
library('WGCNA')
{
  # 样本信息不同特征值赋予不同的颜色
  datExpr_tree<-hclust(dist(datExpr), method = "average")
  par(mar = c(0,5,2,0))
  plot(datExpr_tree, main = "Sample clustering", sub="", xlab="",
       cex.axis = 0.9, cex.main = 1.2, cex.lab=1, cex = 0.7)
  
  
  status_colors <- numbers2colors(as.numeric(factor(hclust_table$condition)), 
                                                                 colors = c("red","white","pink","blue","green"),signed = FALSE)
  par(mar = c(1,4,3,1),cex=0.8)
  plotDendroAndColors(datExpr_tree, cbind(status_colors),
                      groupLabels = colnames(hclust_table),
                      cex.dendroLabels = 0.8,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.01,
                      main = "Sample dendrogram and trait heatmap")
  
}


#确定软阈值
powers = c(c(1:10), seq(from = 12, to=20, by=2))
## [1]  1  2  3  4  5  6  7  8  9 10 12 14 16 18 20
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
# Scale independence
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices$Power,
     -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n", 
     main = paste("Scale independence"))
text(sft$fitIndices$Power,
     -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     labels=powers,
     cex=cex1,
     col="red")
abline(h=0.90,col="RED")

# Mean connectivity
plot(sft$fitIndices$Power,
     sft$fitIndices$mean.k.,
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity",
     type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices$Power,
     sft$fitIndices$mean.k.,
     labels=powers,
     cex=cex1,
     col="red")

best_beta=sft$powerEstimate
best_beta






#第四步构建共表达矩阵
net = blockwiseModules(datExpr, power = 9,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
## 这个诡异的问题必须记录一下，第一次运行blockwiseModules很顺利，后来就开始报错了，网上给的解决方案是重启电脑，实在不行
## 把WGCNA包重装一次，详情看这里https://www.biostars.org/p/305714/

table(net$colors)

moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
# 计算根据模块特征向量基因计算模块相异度：
MEDiss = 1-cor(MEs0);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result

plot(METree, 
     main = "Clustering of module eigengenes",
     xlab = "", 
     sub = "")
# 在聚类图中画出剪切线
MEDissThres = 0.25
abline(h=MEDissThres, col = "red")

merge_modules = mergeCloseModules(datExpr, moduleColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge_modules$colors;
mergedMEs = merge_modules$newMEs;
plotDendroAndColors(net$dendrograms[[1]], cbind(moduleColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)






#第五步TOM图
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
geneTree = net$dendrograms[[1]]
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6)
nSelect = 400
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
sizeGrWindow(9,9)
plotDiss = selectTOM^7
diag(plotDiss) = NA
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")



#第七部模块与性状的关系
design=model.matrix(~0+ group_list)
colnames(design)=levels(factor(group_list))
# 根据相关性重新排序
MEs = orderMEs(MEs0)
# 计算两个相关性值
moduleTraitCor = cor(MEs, design , use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(datExpr))

sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "
(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# 可视化
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels =colnames(design[,1:4]),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

sizeGrWindow(5,7.5);
par(cex = 0.9)
# 特征基因网络图
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, 
                      xLabelsAngle= 90)
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
## 聚类图
plotEigengeneNetworks(MEs, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
## 性状与模块热图
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)




#模块相关性进一步分析
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
## 算出每个模块跟特征基因的皮尔森相关系数矩阵
## MEs是每个模块在每个样本里面的值
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr,use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

module = "yellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for condition",
                   main = paste("Module membership vs. gene significance
"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "purple"
column = match(module, modNames);
moduleGenes = moduleColors==module;
table(moduleColors)
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for condition",
                   main = paste("Module membership vs. gene significance
"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)




# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = 7); 
# Select module
module = "brown";
# Select module probes
probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
View(modProbes)
## 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOMturquoise = TOM[inModule, inModule];
dimnames(modTOMturquoise) = list(modProbes, modProbes)
## 模块对应的基因关系矩阵 
cyt = exportNetworkToCytoscape(
  modTOMturquoise,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes, 
  nodeAttr = moduleColors[inModule]
)




#绘制GO图
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = 7); 
# Select module
library(enrichplot)
library(clusterProfiler)
module = "yellow";
# Select module probes
probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
## 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
## 模块对应的基因关系矩阵 
gene <- bitr(modProbes, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = "hugene10sttranscriptcluster.db")
ego_yellow <- enrichGO(gene =gene$ENTREZID,
                     keyType = 'ENTREZID',
                     OrgDb = hugene10sttranscriptcluster.db,
                     ont = "bp",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.01,
                     readable = TRUE)

cnetplot(ego_yellow,categorySize = "pvalue",foldChange = gene$ENTREZID)


cnetplot(ego_purple, foldChange=NULL, circular = TRUE, colorEdge = TRUE)

goplot(ego_purple)


save(TOM,datExpr,exprSet,design,file = 'WGCNA.Rdata')


library(clusterProfiler)
library(org.Hs.eg.db)    
library(ggplot2)    
module = "turquoise";
# Select module probes
probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
geneNames <- modProbes
#geneNames <- rownames(DEGs_blood)

gene <-  mapIds(org.Hs.eg.db, geneNames, 'ENTREZID', 'SYMBOL')  


#GO分析和作图
library(clusterProfiler)
library(org.Hs.eg.db)

geneNames <- rownames(d10[d10$module == 1,])
gene <-  mapIds(org.Hs.eg.db, geneNames, 'ENTREZID', 'SYMBOL')  
GO<-enrichGO(gene = gene, OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",ont="ALL",qvalueCutoff = 0.05,readable = T)
go<-as.data.frame(GO)
View(go)
write.csv(go, file = "C:/Users/Administrator/Desktop/论文/脓毒血症，脓毒性心肌病，建模/supplement-返修/GoTurquoiseinBlood.csv")         

table(go[,1])
go_MF<-go[go$ONTOLOGY=="MF",][1:10,]
go_CC<-go[go$ONTOLOGY=="CC",][1:10,]
go_BP<-go[go$ONTOLOGY=="BP",][1:10,]
go_enrich_df<-data.frame(ID=c(go_BP$ID, go_CC$ID, go_MF$ID),
                                                   Description=c(go_BP$Description, go_CC$Description, go_MF$Description),
                                                   GeneNumber=c(go_BP$Count, go_CC$Count, go_MF$Count),
                                                   type=factor(c(rep("biological process", 10), 
                                                   rep("cellular component", 10),
                                                   rep("molecular function",10)),
                                                   levels=c("molecular function", 
                                                   "cellular component", "biological process")))
                             
## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## shorten the names of GO terms
shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}

labels=(sapply( levels(factor(go_enrich_df$Description))[as.numeric(go_enrich_df$number)], shorten_names))
names(labels) = rev(1:nrow(go_enrich_df))

## colors for bar // green, blue, orange
CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
library(ggplot2)
p <- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  scale_fill_manual(values = CPCOLS) + theme_test() + 
  scale_x_discrete(labels=labels) +
  xlab("GO term") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")
#coord_flip(...)横向转换坐标：把x轴和y轴互换，没有特殊参数
p



#综合柱状图
p1 <- ggplot(data=goAll)+  geom_bar(aes(x=Description,y=-log10(pvalue), fill=GOType), stat='identity') + coord_flip() + scale_x_discrete(limits=goAll$Description) 

ggsave("out_bar.pdf", p1, width = 10, height=6)


p2 <- ggplot(Edata, aes(x=GeneRatio, y=`GO description`)) +
  geom_point(aes( size= Count , colour = -log10( pvalue ))  ) + scale_y_discrete(limits=Edata$`GO description`)+
  ggtitle("GO enrichment")  +  scale_color_gradient(low = 'green', high = 'red') + xlim(range(Edata$GeneRatio)) +
  theme(axis.text.x=element_text(angle=0,size=8, vjust=0.7), axis.text.y=element_text(angle=0,size=6, vjust=0.7),plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16), panel.background = element_rect(fill="white", colour='gray'), panel.grid.major = element_line(size = 0.05, colour = "gray"), panel.grid.minor.y = element_line(size=0.05, colour="gray"), panel.grid.minor.x = element_line(size=0.05, colour="gray")
  )

ggsave("out_GO.pdf", p2, width = 8, height=7)

#另一种作图
geneNames <- rownames(d4[d4$module == c("brown"),])
gene <-  mapIds(org.Hs.eg.db, geneNames, 'ENTREZID', 'SYMBOL')  
GO_brown <- enrichGO(   gene   = gene,    
            OrgDb  = org.Hs.eg.db,    
            ont   = "BP"  ,    
            pAdjustMethod = "BH",    
            pvalueCutoff  = 0.01,    
            qvalueCutoff  = 0.05) 
GO_brown <- GObrown
View(GO_brown@result)
head(GO_turquoise@result)
dotplot(GO_turquoise, showCategory = 10)
write.csv(GO_brown, file = "D:/R FILE/pathomyocarditis/GO_cobrown.csv")
save(GO_blue, GO_brown, GO_turquoise, file = "coGO.Rdata")
#绘图
library(forcats)
library(ggplot2)
f <- GO_blue@result[1:5,]
ggplot(f, aes(x = Description, y = pvalue))+
  geom_bar(stat = "identity") + 
  labs(x = "", y = "pvalue") + 
  coord_flip()





d10 <- d1[gene_name,]



#提取基因名称
module = "blue";
# Select module probes
DEGs_blood <- nrDEG_65682[nrDEG_65682$P.Value < 0.01,]
probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
geneNames <- modProbes
gene_blue <- geneNames
gene_yellow <- geneNames
gene_blood <- rownames(DEGs_blood)
gene_heart <- rbind(gene_blue,gene_brown,gene_turquoise,gene_yellow)
gene_cross_blue <- intersect(gene_blood,gene_blue)
gene_cross_yellow <- intersect(gene_blood,gene_yellow)
gene_cross_brown <- intersect(gene_blood,gene_brown)
gene_cross_turquoise <- intersect(gene_blood,gene_turquoise)
save(gene_cross_blue,gene_cross_brown,gene_cross_turquoise,gene_cross_yellow,
     file = 'gene_cross.Rdata')

gene_name <- c(gene_blue,gene_brown,gene_turquoise,gene_yellow)


#提取共同方向表达的基因
moduleColors
d1$module <- moduleColors
e1 <- d1[d1$logFC > 0,]
e2 <- d1[d1$logFC < 0,]
e3 <- DEGs_blood[DEGs_blood$logFC > 0,]
e4 <- DEGs_blood[DEGs_blood$logFC < 0,]
intersect(rownames(e1),rownames(e3))
intersect(rownames(e2),rownames(e4))
union(intersect(rownames(e1),rownames(e3)),intersect(rownames(e2),rownames(e4)))
e5 <- union(intersect(rownames(e1),rownames(e3)),intersect(rownames(e2),rownames(e4)))
d4 <- d1[c(e5),]
table(d4$module)
rownames(d4[d4$module == c("blue"),])
f3 <- rownames(d4[d4$module == c("brown"),])
write.csv(f3, file = "D:/R FILE/pathomyocarditis/f3.csv")


f1 <- rownames(d4[d4$module == c("blue"),])
write.csv(f1, file = "D:/R FILE/pathomyocarditis/f1.csv")

gene_name <- rownames(d4)
write.csv(gene_name, file = "D:/R FILE/pathomyocarditis/gene_name.csv")
save(gene_name, file = "gene_name.Rdata")


write.csv(d10, file = "C:/Users/Administrator/Desktop/论文/脓毒血症，脓毒性心肌病，建模/supplement-返修2/ModulesinBlood.csv")         
