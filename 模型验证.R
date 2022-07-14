gene <- gene_54514[,c("HACD1" , "DECR1", 
                      "TNFRSF8","DAAM2", "V13630", "V13630" )]
View(gene[rownames(test_54514_new),])
gene <- gene[rownames(test_54514_new),]
gene$mortality <- test_54514_new$group
gene$age <- test_54514_new$age
gene$appachII <- test_54514_new$appachII
gene[,1:4] <- sapply(gene[,1:4], as.numeric)
fit.logit <- glm(mortality~ HACD1+DECR1+TNFRSF8+DAAM2, data = gene, family = binomial())
summary(fit.logit)
prob <- predict(fit.logit, gene, type = "response")
logit.pred <- factor(prob > .5, levels = c(FALSE, TRUE),
                     labels = c("survive", "death"))
logit.perf <- table(gene$mortality, logit.pred,
                    dnn = c("Actual", "Predicted"))
logit.perf
library(pROC)
library(plotROC)
roc(gene$mortality, gene$age, plot=TRUE, print.thres=TRUE, print.auc=TRUE)

#多个曲线
age <- pd[rownames(test),]
age <- age$age.ch1
age <- as.data.frame(age)
rownames(age) <- rownames(test)
test_roc_mortality<- cbind(gene$mortality,gene$age,gene$appachII,prob)
test_roc_mortality$gene_model <-  prob
colnames(test_roc_mortality) <- c("mortality","age","appachII","gene_model")
longtest <- melt_roc(test_roc_mortality , "mortality", c("age","appachII","gene_model"))
head(longtest)
ggplot(longtest, aes(d = D, m = M, color = name)) + geom_roc() + style_roc()


#nomograph
#nomograph
library(rms)
fit1 <- lrm(mortality ~HACD1+DECR1+TNFRSF8+DAAM2+age+appachII, data = gene,
            x = T, y = T)
fit1
dd <- datadist(gene)
options(datadist = "dd")
nom1 <- nomogram(fit1, fun=plogis,fun.at=c(.001, .01, .05, seq(.1,.9, by=.1), .95, .99, .999),lp=F, funlabel="Risk")  
plot(nom1)
cal1 <- calibrate(fit1, method='boot', B=1000) 
plot(cal1,xlim=c(0,1.0),ylim=c(0,1.0))


install.packages("regplot")
library(regplot)
regplot(fit1, 
        observation=pbc[1,], 
        odds=TRUE, 
        interval="confidence")










