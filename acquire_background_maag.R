library(limma)
library(dplyr)
library(ENmix)
#setwd("C:/GRA/CodingProjects/PlasticityInducedEpiChanges/Data/Maag/Maag_raw_data_all_samples/")

## Read in data 
sdrf <- read.delim("sample_data_relationship.txt")
treatment <- sdrf[, "Characteristics.treatment."]
treatment <- factor(treatment, levels = c("control", "HFS_30m", "HFS_2h", "HFS_5h"))

all_raw_data <-  c("Ramaciotti_257549210001_S01_CGH_107_Sep09_1_1.txt", 
                   "Ramaciotti_257549210001_S01_CGH_107_Sep09_1_2.txt",
                   "Ramaciotti_257549210001_S01_CGH_107_Sep09_1_3.txt",
                   "Ramaciotti_257549210001_S01_CGH_107_Sep09_1_4.txt",
                   "Ramaciotti_257549210002_S01_CGH_107_Sep09_1_1.txt",
                   "Ramaciotti_257549210002_S01_CGH_107_Sep09_1_2.txt",
                   "Ramaciotti_257549210002_S01_CGH_107_Sep09_1_3.txt",
                   "Ramaciotti_257549210002_S01_CGH_107_Sep09_1_4.txt")


all_raw_data <- read.maimages(files = all_raw_data, source = 'agilent', other.columns = "gIsWellAboveBG")


## Filter unwanted probes
control_probes <- all_raw_data$genes$ControlType==1
neg_control_probes <- all_raw_data$genes$ControlType==-1
unexpr_probes <- rowSums(all_raw_data$other$gIsWellAboveBG > 0) >= 4 # if use this then get no diff expr genes?
data_filt <- all_raw_data[!control_probes & !neg_control_probes,]


## Normalise within arrays & between arrays
raw_data_normalised <- normalizeWithinArrays(data_filt, method = "loess")
data_norm_filt <- normalizeBetweenArrays(raw_data_normalised, method = "quantile")


## Differential expression
design <- model.matrix(~treatment)
fit <- lmFit(data_norm_filt, design = design)
fit <- eBayes(fit, trend=T, robust = T)
summary(decideTests(fit[,-1]))

results <- topTable(fit = fit, coef = 4, n = length(as.vector(fit$genes$ProbeUID)))


## Gene ontology analysis
g <- goana(fit, coef=4, species="Rn", geneid = "EntrezID")
GO_results <- topGO(g, n=20, truncate="50")