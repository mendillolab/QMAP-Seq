require(dplyr)
require(tidyr)
require(tibble)
require(genefilter)
require(ggplot2)
require(readr)
require(plotly)

# parameters
dose4.sgNT.mean.cutoff = 90
sgNT.mean.totDiff.cutoff = 25

# initialize variables needed to generate exclude list (list of CLxDrugxGene's to remove)
CL.drug.sgRNA.exclude <- c()
sgRNAs <- c("sgATF3", "sgATF4", "sgATF6", 
            "sgATG7", "sgERN1", "sgHSF1", "sgHSF2", 
            "sgKEAP1", "sgNFE2L2", "sgSLC35F2", "sgXBP1") 
n_sgRNAs <- length(sgRNAs)

# load QMAPP dataframes
auc <- read_csv("Downloads/AUC_LongFormat.csv")
qmapp.long <- read_csv("QMAPP/QMAPP.concatenated.long.csv")

# identify combinations where highest dose is ineffective against non-targeting
CL.drug.exclude <- qmapp.long %>% 
  filter(dose=="Dose4", sgRNA=="sgNT0.4") %>% 
  unite(col="sgNT",sgRNA, rep, sep=".") %>% 
  spread(sgNT, survival) %>% 
  unite("CL.drug", cell_line, drug, sep=".") %>% 
  mutate(sgNT.mean = (sgNT0.4.Rep1+sgNT0.4.Rep2)/2) %>%
  filter(sgNT.mean >= dose4.sgNT.mean.cutoff) %>% 
  select(CL.drug) %>% unlist 

# Add all genes for these CLxDrug combos to exclude list 
sgRNAs.rep <- rep(sgRNAs, length(CL.drug.exclude))
CL.drug.exclude <- rep(CL.drug.exclude, each=n_sgRNAs)
CL.drug.sgRNA.exclude <- c(CL.drug.sgRNA.exclude, paste(CL.drug.exclude, sgRNAs.rep, sep="."))

# identify combinations where overall change is small for non-targeting
sgNT.dose1 <- qmapp.long %>%
  filter(dose=="Dose1", sgRNA=="sgNT0.4") %>% 
  unite(col="sgNT.dose.rep", sgRNA, dose, rep, sep=".") %>%
  spread(sgNT.dose.rep, survival) %>% 
  select(-Concentration)

sgNT.dose4 <- qmapp.long %>%
  filter(dose=="Dose4", sgRNA=="sgNT0.4") %>% 
  unite(col="sgNT.dose.rep", sgRNA, dose, rep, sep=".") %>%
  spread(sgNT.dose.rep, survival) %>% 
  select(-Concentration)

CL.drug.exclude <- merge(sgNT.dose1, sgNT.dose4, by=c("cell_line", "drug")) %>%
  mutate(sgNT.dose1.mean = (sgNT0.4.Dose1.Rep1 + sgNT0.4.Dose1.Rep2)/2, 
         sgNT.dose4.mean = (sgNT0.4.Dose4.Rep1 + sgNT0.4.Dose4.Rep2)/2,
         sgNT.totDiff = sgNT.dose1.mean - sgNT.dose4.mean) %>%
  filter(sgNT.totDiff < sgNT.mean.totDiff.cutoff) %>% 
  unite(col="CL.drug", cell_line, drug, sep=".") %>% 
  select(CL.drug) %>% unlist

sgRNAs.rep <- rep(sgRNAs, length(CL.drug.exclude))
CL.drug.exclude <- rep(CL.drug.exclude, each=n_sgRNAs)
CL.drug.sgRNA.exclude <- c(CL.drug.sgRNA.exclude, paste(CL.drug.exclude, sgRNAs.rep, sep=".")) %>% unique
        
# remove drugs from well B7
auc <-auc %>% filter(Drug != "Ruxolitinib", Drug != "STF083010")

# remove ZR-75-1 sgKEAP1 &  ZR-75-1 sgHSF1
auc <- auc %>% filter(!(CellLine == "ZR-75-1" & (sgRNA=="sgKEAP1" | sgRNA=="sgHSF1")))

# reshape >> CellLine, Drug, sgNT.1 (AUC), sgNT.2 (AUC) 
sgNT.df <- auc %>% 
  filter(sgRNA=="sgNT") %>%
  unite(col="sgNT",sgRNA, Rep, sep=".") %>% 
  spread(sgNT, AUC)

# reshape >> CellLine, Drug, sgRNA, AUC.1, AUC.2
sgRNA.df <- auc %>% filter(sgRNA !="sgNT") %>% spread(Rep, AUC)
colnames(sgRNA.df) <- c("CellLine", "Drug", "sgRNA", "AUC.1", "AUC.2")

# reshape >> CellLine, Drug, sgRNA, AUC.1, AUC.2, sgNT.1, sgNT.2 
auc.reshaped <- merge(sgRNA.df, sgNT.df, by=c("CellLine","Drug"))

# reshape >> CellLine.Drug.sgRNA, AUC.1, AUC.2, sgNT.1, sgNT.2
auc.reshaped <- auc.reshaped %>% unite(col="CellLine.Drug.sgRNA", CellLine, Drug, sgRNA, sep=".")

# filter out combinations on the exclude list
auc.reshaped <- auc.reshaped %>% filter(!(CellLine.Drug.sgRNA %in% CL.drug.sgRNA.exclude))

# run t-tests (AUC vs sgNT)
auc.test <- auc.reshaped %>% column_to_rownames("CellLine.Drug.sgRNA") %>% as.matrix
g <- factor( c(0,0,1,1) ) 
test <- rowttests(auc.test, g) 
test[,"CellLine.Drug.sgRNA"] <- rownames(test)

# combine t-test results with reshaped dataframe
auc.reshaped <- merge(auc.reshaped, test, by=c("CellLine.Drug.sgRNA"))

# make volcano plot
auc.reshaped <- auc.reshaped %>% mutate(AUC.mean = (AUC.1+AUC.2)/2,
                                        sgNT.mean = (sgNT.1+sgNT.2)/2)

# fold changes (not log2)
auc.reshaped <- auc.reshaped %>% mutate(foldchange.div = AUC.mean / sgNT.mean,
                                        foldchange.sub = AUC.mean - sgNT.mean,
                                        foldchange.sub.div = (AUC.mean - sgNT.mean)/sgNT.mean,
                                        p.neglog10 = -1 * log10(p.value))

# save 
write_csv(auc.reshaped, path="QMAPP.filtered.tests.csv")
# (MDA-MB-231.Ruxolitinib MDA-MB-231.STF083010)

# p1 <- ggplot(auc.reshaped, aes(x=foldchange.div, y=p.neglog10)) + geom_point(aes(text=CellLine.Drug.sgRNA), alpha=0.3, size=0.8)
# p2 <- ggplot(auc.reshaped, aes(x=foldchange.sub, y=p.neglog10) ) + geom_point(aes(text=CellLine.Drug.sgRNA), alpha=0.3, size=0.8)
# p3 <- ggplot(auc.reshaped, aes(x=foldchange.sub.div, y=p.neglog10)) + geom_point(aes(text=CellLine.Drug.sgRNA), alpha=0.3, size=0.8)
# ggplotly(p2)

#fold changes log2
# auc.reshaped <- auc.reshaped %>% mutate(foldchange.div = log((AUC.mean / sgNT.mean)+2),
#                                         foldchange.sub = log((AUC.mean - sgNT.mean)+2),
#                                         foldchange.sub.div = log(((AUC.mean - sgNT.mean)/sgNT.mean)+2),
#                                         p.neglog10 = -1 * log10(p.value))
# 
# p1 <- ggplot(auc.reshaped, aes(x=foldchange.div, y=p.neglog10)) + geom_point(aes(text=CellLine.Drug.sgRNA), alpha=0.3, size=0.8)
# p2 <- ggplot(auc.reshaped, aes(x=foldchange.sub, y=p.neglog10) ) + geom_point(aes(text=CellLine.Drug.sgRNA), alpha=0.3, size=0.8)
# p3 <- ggplot(auc.reshaped, aes(x=foldchange.sub.div, y=p.neglog10)) + geom_point(aes(text=CellLine.Drug.sgRNA), alpha=0.3, size=0.8)
# ggplotly(p2)

# remove positive control (outlier)
#auc.reshaped.noPC <- auc.reshaped %>% filter(CellLine.Drug.sgRNA != "MDA-MB-231.YM155.sgSLC35F2")
