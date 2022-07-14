rm(list = ls())

#获取心脏和血液中的共表达基因

e1 <- DEGs_blood[DEGs_blood$logFC > 0, ]
e2 <- DEGs_heart[DEGs_heart$logFC > 0, ]
co_upgene <- intersect(rownames(e1), rownames(e2))
co_upgene <- intersect(co_upgene,rownames(gene_65682))
e3 <- DEGs_blood[DEGs_blood$logFC < 0, ]
e4 <- DEGs_heart[DEGs_heart$logFC < 0, ]
co_downgene <- intersect(rownames(e3), rownames(e4))
co_downgene <- intersect(co_downgene, rownames(gene_65682))
coeprGene <- union(co_upgene, co_downgene)
#coeprGene为在脓毒症心脏和血液标本中共同上调或下调表达的标本


#合并基因和分组情况
gene_name <- read.csv(gene_name.csv)
x1 <- read.csv("D:/R FILE/GSE65682/model construction/ModuleColorinBlood_all.csv")
coeprGene <- intersect(gene_name, rownames(gene_65682))
gene_65682[coeprGene,] <- sapply(gene_65682[coeprGene,], as.numeric)
a1 <- rbind(gene_65682[coeprGene,],geneGroup_65682[,10860])
a1 <- t(a1)
a1$mortality <-  pd_65682$`mortality_event_28days:ch1`
rownames(a1) <- c(x1$geneName,"group")
a1 <- sapply(a1[1:400],as.numeric)
a1 <- data.frame(a1)




#合并基因和死亡事件
colnames(pd_65682)
b1 <- pd_65682[c(16,17)]
b1$mortality[b1$characteristics_ch1.6 == c("mortality_event_28days: 1")] <- 1
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 1")] <- c(1)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 2")] <- c(2)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 20")] <- c(20)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 21")] <- c(21)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 22")] <- c(22)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 23")] <- c(23)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 24")] <- c(24)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 25")] <- c(25)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 26")] <- c(26)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 27")] <- c(27)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 28")] <- c(28)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 10")] <- c(10)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 11")] <- c(11)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 12")] <- c(12)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 13")] <- c(13)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 14")] <- c(14)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 15")] <- c(15)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 16")] <- c(16)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 17")] <- c(17)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 18")] <- c(18)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 19")] <- c(19)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 3")] <- c(3)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 4")] <- c(4)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 5")] <- c(5)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 6")] <- c(6)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 7")] <- c(7)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 8")] <- c(8)
b1$time2death[b1$characteristics_ch1.7 == c("time_to_event_28days: 9")] <- c(9)
b1 <- b1[c(3,4)]
b1$mortality[(b1$mortality) == 0]<- c("survive")
b1$mortality[(b1$mortality) == 1]<- c("death")



#合并基因分组和死亡信息
a1 <- t(a1)
d1 <- cbind(a1, b1)
View(d1)
#d1是共表达基因在血液标本中的表达情况，合并了分组和生存情况，生存时间


#选择出脓毒血症组的病例
#d1 <- d1[d1$group==c('sepsis'),]
#出去空值
#d1 <- na.omit(d1)
#View(d1)
#dim(d1)
#colnames(d1)
data1 <- d1[,1:1049]#data1是去除了分组信息的共表达血液标本的情况
data1 <- t(data1)
#data1 <- sapply(data1,as.numeric)
rownames(data1) <- rownames(design)
#计算共表达基因中与死亡相关的基因
{
  design <- model.matrix( ~0 + factor( d1[,1051] ) )
  colnames( design ) = levels( factor( d1[,1051] ))
  rownames( design ) = rownames(d1)
}
colnames(design) <- c('survivor','nosurvivor')
library(limma)

contrast.matrix <- makeContrasts( "survivor-nosurvivor", levels = design )
contrast.matrix

## design和contrast.matrix都是为了下面的差异分析制作的

{
  fit <- lmFit( data1, design )
  fit2 <- contrasts.fit( fit, contrast.matrix ) 
  fit2 <- eBayes( fit2 )
  nrDEG = topTable( fit2, coef = 1, n = Inf ) 
  write.table( nrDEG, file = "nrDEG_gse65682_sepsis.out")
}
head(nrDEG)
View(nrDEG[nrDEG$P.Value < 0.01,])
DEGs <-nrDEG[nrDEG$P.Value < 0.01,]
write.csv(DEGs, file = "C:/Users/Administrator/Desktop/论文/脓毒血症，脓毒性心肌病，建模/supplement-返修2/DEGsForClustering.csv")         

DEGs_down <- nrDEG[nrDEG$P.Value < 0.01&nrDEG$logFC <= -1,]
DEGs_up <- nrDEG[nrDEG$P.Value < 0.01&nrDEG$logFC >= 1,]
DEGs_final <- rbind(DEGs_down,DEGs_up )

save(DEGs,file = 'DEGs_final.Rdata')

dim(DEGs)

geneformodel <- rownames(nrDEG[nrDEG$adj.P.Val < 0.05,])
geneformodel <- geneformodel[1:105]
save(geneformodel, file = "geneformodel.Rdata")
#建立模型效率检验函数performance
performance <- function(table, n = 2){
  if(!all(dim(table) == c(2,2)))
    stop("Must be a 2*2 table")
  tn = table[1,1]
  fp = table[1,2]
  fn = table[2,1]
  tp = table[2,2]
  sensitivity = tp/(tp+fn)
  specificity = tn/(tn+fp)
  ppp = tp/(tp+fp)
  npp = tn/(tn+fn)
  hitrate = (tp+tn)/(tp+tn+fp+fn)
  result <- paste("Sensitivity = ", round(sensitivity, n),
                  "\nSpecifity = ", round(specificity, n),
                  "\nPositive Predictive Value = ", round(ppp, n),
                  "\nNegative Predictive Value = ", round(npp, n),
                  "\nAccuracy = ", round(hitrate, n),"\n", sep = "")
  cat(result)
}





#用于多个模型效率的对比
#区分训练集和验证集

data1 <- data1[rownames(DEGs_2),]
d2 <-  as.data.frame(d1[,402])
colnames(d2) <- c('mortality')
data2 <- cbind(data1, d2)
colnames(data2) <- c(colnames(data1),"mortality")
data2 <- na.omit(data2)
class(data2)
data2$mortality <- factor(data2$mortality , levels=c("survive", "death"), ordered = TRUE)

name <- rownames(DEGs[DEGs$abslogFC >= 0.35,])
data3 <- data2[,c(name, 'mortality')]











set.seed(1234)
data4 <- scale(data4)
ind <- sample(2, nrow(data4), replace = TRUE, prob = c(0.3,0.7))
train <- data4[ind == 1,]
test <- data4[ind == 2,]

train <- data.frame(train)
train <- sapply(train, as.numeric)
train <- cbind(train, train_mortality)
test <- sapply(test, as.numeric)
test <- data.frame(test)

library(randomForest)
set.seed(1234)
fit.forest <- randomForest(mortality~., data = train,importance = TRUE
                           )
                          
fit.forest
summary(fit.forest)

forest.pred <- predict(fit.forest,test1)
forest.perf <- table(test1$mortality,forest.pred,
                     dnn = c('Actual', 'Predicted'))
forest.perf
importance(forest.pred)



fit.logit <- glm(mortality~CD24+LDLR, data = train, family = binomial())
summary(fit.logit)
prob <- predict(fit.logit, train, type = "response")
logit.pred <- factor(prob > .5, levels = c(FALSE, TRUE),
                     labels = c("survive", "death"))
logit.perf <- table(train$mortality, logit.pred,
                    dnn = c("Actual", "Predicted"))
logit.perf


#支持向量机
library(e1071)
set.seed(1234)
fit.svm <- svm(mortality~., data = train)
fit.svm
svm.pred <- predict(fit.svm, na.omit(test))
svm.perf <- table(test$mortality, svm.pred, dnn = c('Actual', 'Predicted'))
svm.perf

#决策树
library(rpart)
set.seed(1234)
dtree <- rpart(class ~ ., data = train, method = 'class',
               parms = list(split = "information"))


#热图
data_heatmap <- gene_79962[c('SMU1','NNMT','TMEM184C','TBC1D22B','DECR1','VPS45'),]
data_heatmap <- as.data.frame(scale(data_heatmap))


#lasso回归
library(glmnet)
x <- as.matrix(d7[,c(1:105)])
y <- as.numeric(d7[,106])
lasso_fit <- glmnet(x,y, alpha = 1, family = "gaussian")
plot(lasso_fit, xvar = "lambda", label = TRUE)
coef(lasso_fit, s = lasso_fit$lambda[30])
print(lasso_fit)
h <- coef(lasso_fit, s = lasso_fit$lambda[26])
h <- as.matrix(h)
h <- as.data.frame(h)
h$gene <- rownames(h)
h1 <- h[h$s1 !=0,]
rownames(h1)
cvfit <- cv.glmnet(x,y,alpha = 1, family = "gaussian")
plot(cvfit)
cvfit$lambda.1se
cvfit$lambda.min
coef1 <- coef(cvfit$glmnet.fit, s = 0.0127231, exact = F)
coef1 


genefromforest <- read.csv("genefromforest_new.csv")
intersect(genefromforest$x, rownames(h1))
