ibrary(dplyr)
library(tidyr)
library(tibble)
library(genefilter)
library(ggplot2)
library(plotly)
​
# load AUC dataframe (CellLine, Drug, sgRNA, Rep, AUC)
auc <- read_csv("Downloads/AUC_LongFormat.csv")
​
# remove drugs from well B7
auc <-auc %>% filter(Drug != "Ruxolitinib", Drug != "STF083010")
​
# remove ZR-75-1 sgKEAP1 &  ZR-75-1 sgHSF1
auc <- auc %>% filter(!(CellLine == "ZR-75-1" & (sgRNA=="sgKEAP1" | sgRNA=="sgHSF1")))
​
  
# reshape >> CellLine, Drug, sgNT.1 (AUC), sgNT.2 (AUC) 
sgNT.df <- auc %>% 
  filter(sgRNA=="sgNT") %>%
  unite(col="sgNT",sgRNA, Rep, sep=".") %>% 
  spread(sgNT, AUC)
​
# reshape >> CellLine, Drug, sgRNA, AUC.1, AUC.2
sgRNA.df <- auc %>% filter(sgRNA !="sgNT") %>% spread(Rep, AUC)
colnames(sgRNA.df) <- c("CellLine", "Drug", "sgRNA", "AUC.1", "AUC.2")
​
# reshape >> CellLine, Drug, sgRNA, AUC.1, AUC.2, sgNT.1, sgNT.2 
auc.reshaped <- merge(sgRNA.df, sgNT.df, by=c("CellLine","Drug"))
​
# reshape >> CellLine.Drug.sgRNA, AUC.1, AUC.2, sgNT.1, sgNT.2
auc.reshaped <- auc.reshaped %>% unite(col="CellLine.Drug.sgRNA", CellLine, Drug, sgRNA, sep=".")
​
# run t-tests (AUC vs sgNT)
auc.test <- auc.reshaped %>% column_to_rownames("CellLine.Drug.sgRNA") %>% as.matrix
g <- factor( c(0,0,1,1) ) 
test <- rowttests(auc.test, g) 
test[,"CellLine.Drug.sgRNA"] <- rownames(test)
​
# combine t-test results with reshaped dataframe
auc.reshaped <- merge(auc.reshaped, test, by=c("CellLine.Drug.sgRNA"))
​
# make volcano plot
auc.reshaped <- auc.reshaped %>% mutate(AUC.mean = (AUC.1+AUC.2)/2,
                                        sgNT.mean = (sgNT.1+sgNT.2)/2)
​
# fold changes (not log2)
auc.reshaped <- auc.reshaped %>% mutate(foldchange.div = AUC.mean / sgNT.mean,
                                        foldchange.sub = AUC.mean - sgNT.mean,
                                        foldchange.sub.div = (AUC.mean - sgNT.mean)/sgNT.mean,
                                        p.neglog10 = -1 * log10(p.value))
​
# (MDA-MB-231.Ruxolitinib MDA-MB-231.STF083010)
​
p1 <- ggplot(auc.reshaped, aes(x=foldchange.div, y=p.neglog10)) + geom_point(aes(text=CellLine.Drug.sgRNA), alpha=0.3, size=0.8)
p2 <- ggplot(auc.reshaped, aes(x=foldchange.sub, y=p.neglog10) ) + geom_point(aes(text=CellLine.Drug.sgRNA), alpha=0.3, size=0.8)
p3 <- ggplot(auc.reshaped, aes(x=foldchange.sub.div, y=p.neglog10)) + geom_point(aes(text=CellLine.Drug.sgRNA), alpha=0.3, size=0.8)
ggplotly(p1)
​
# fold changes log2
auc.reshaped <- auc.reshaped %>% mutate(foldchange.div = log((AUC.mean / sgNT.mean)+2),
                                        foldchange.sub = log((AUC.mean - sgNT.mean)+2),
                                        foldchange.sub.div = log(((AUC.mean - sgNT.mean)/sgNT.mean)+2),
                                        p.neglog10 = -1 * log10(p.value))
​
p1 <- ggplot(auc.reshaped, aes(x=foldchange.div, y=p.neglog10)) + geom_point(aes(text=CellLine.Drug.sgRNA), alpha=0.3, size=0.8)
p2 <- ggplot(auc.reshaped, aes(x=foldchange.sub, y=p.neglog10) ) + geom_point(aes(text=CellLine.Drug.sgRNA), alpha=0.3, size=0.8)
p3 <- ggplot(auc.reshaped, aes(x=foldchange.sub.div, y=p.neglog10)) + geom_point(aes(text=CellLine.Drug.sgRNA), alpha=0.3, size=0.8)
ggplotly(p1)
​
# remove positive control (outlier)
auc.reshaped.noPC <- auc.reshaped %>% filter(CellLine.Drug.sgRNA != "MDA-MB-231.YM155.sgSLC35F2")
