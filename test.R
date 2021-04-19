# script: test.R
# filter AUCs and run t-tests
# input: path for qmap.AUC and qmap.long
# output: path to store test folder [[can be same folder]]

require(dplyr)
require(tidyr)
require(tibble)
require(genefilter)
require(readr)
# require(plotly)
# require(ggplot2)

#### handle file paths #####
args <- commandArgs()
input <- sub('--input=', '', args[grep('--input=', args)])
output <- sub('--output=', '', args[grep('--output=', args)])
#input <- "../results"
#output <- "../results"
qmap.auc <- read_csv(paste0(output,"/QMAP.concatenated.AUC.csv"))
qmap.long <- read_csv(paste0(output,"/QMAP.concatenated.long.csv"))

# parameters
dose4.sgNT.mean.cutoff = 90
sgNT.mean.totDiff.cutoff = 25

# initialize variables needed to generate exclude list (list of CLxDrugxGene's to remove)
CL.drug.sgRNA.exclude <- c() 
sgRNAs <- c("sgATF3", "sgATF4", "sgATF6", 
            "sgATG7", "sgERN1", "sgHSF1", "sgHSF2", 
            "sgKEAP1", "sgNFE2L2", "sgSLC35F2", "sgXBP1" ) 
n_sgRNAs <- length(sgRNAs)

# identify combinations where highest dose is ineffective against non-targeting
CL.drug.exclude <- qmap.long %>% 
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
sgNT.dose1 <- qmap.long %>%
  filter(dose=="Dose1", sgRNA=="sgNT0.4") %>% 
  unite(col="sgNT.dose.rep", sgRNA, dose, rep, sep=".") %>%
  spread(sgNT.dose.rep, survival) 

sgNT.dose4 <- qmap.long %>%
  filter(dose=="Dose4", sgRNA=="sgNT0.4") %>% 
  unite(col="sgNT.dose.rep", sgRNA, dose, rep, sep=".") %>%
  spread(sgNT.dose.rep, survival)

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
qmap.auc <- qmap.auc %>% filter(drug != "Ruxolitinib", drug != "STF083010")

# remove ZR-75-1 sgKEAP1 &  ZR-75-1 sgHSF1 #! colnames
qmap.auc <- qmap.auc %>% filter(!(cell_line == "ZR-75-1" & (sgRNA=="sgKEAP1" | sgRNA=="sgHSF1")))

# reshape into test dataframe [cell_line.drug.sgRNA] [AUC.1] [AUC.2] [sgNT.1] [sgNT.2]
sgNT.df <- qmap.auc %>% filter(sgRNA=="sgNT0.4")
colnames(sgNT.df) <- c("cell_line", "drug", "sgRNA", "sgNT.1", "sgNT.2")
sgRNA.df <- qmap.auc %>% filter(sgRNA != "sgNT0.4")
qmap.forTest <- merge(sgRNA.df, sgNT.df, by=c("cell_line","drug")) %>% select(-sgRNA.y)
qmap.forTest <- qmap.forTest %>% unite(col="cell_line.drug.sgRNA", cell_line, drug, sgRNA.x, sep=".")

# filter out combinations on the exclude list
qmap.forTest <- qmap.forTest %>% filter(!(cell_line.drug.sgRNA %in% CL.drug.sgRNA.exclude))

# run t-tests (AUC vs sgNT)
qmap.test <- qmap.forTest %>% column_to_rownames("cell_line.drug.sgRNA") %>% as.matrix
g <- factor( c(0,0,1,1) ) 
qmap.test <- rowttests(qmap.test, g) 
qmap.test[,"cell_line.drug.sgRNA"] <- rownames(qmap.test)

# combine t-test results with reshaped dataframe
qmap.test <- merge(qmap.forTest, qmap.test, by=c("cell_line.drug.sgRNA"))

# make volcano plot
qmap.test <- qmap.test %>% mutate(AUC.mean = (AUC.1+AUC.2)/2,
                                        sgNT.mean = (sgNT.1+sgNT.2)/2)

# fold changes (not log2)
qmap.test <- qmap.test %>% mutate(foldchange.div = AUC.mean / sgNT.mean,
                                        foldchange.sub = AUC.mean - sgNT.mean,
                                        foldchange.sub.div = (AUC.mean - sgNT.mean)/sgNT.mean,
                                        p.neglog10 = -1 * log10(p.value))

# save 
write_csv(qmap.test, path=paste0(output,"/QMAP.filtered.tests.csv"))

print("Finished test.R")
