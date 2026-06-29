rm(list=ls())

###############################################
# Installing required packages
###############################################
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ChAMP")
BiocManager::install("methylGSA")

library("ChAMP")
library("methylGSA")

###############################################
# Setting path
# Directory must contain *.idat files and *.csv file with phenotype
###############################################
path <- "path/to/directory/with/files"
setwd(path)

arraytype <- "EPICv2"
detPcut <- 0.01

###############################################
# Import and filtration
###############################################
myLoad <- champ.load(
  directory = path,
  arraytype = arraytype,
  method = "minfi",
  methValue = "B",
  autoimpute = TRUE,
  filterDetP = TRUE,
  ProbeCutoff = 0.1,
  SampleCutoff = 0.1,
  detPcut = detPcut,
  filterBeads = TRUE,
  beadCutoff = 0.05,
  filterNoCG = TRUE,
  filterSNPs = TRUE,
  filterMultiHit = TRUE,
  filterXY = FALSE,
  force = TRUE
)
pd <- as.data.frame(myLoad$pd)

###############################################
# Functional normalization
###############################################
myNorm <- getBeta(preprocessFunnorm(myLoad$rgSet))
cpgs <- intersect(rownames(myLoad$beta), rownames(myNorm))
myNorm <- myNorm[cpgs, ]
rownames(myNorm) <- cpgs
colnames(myNorm) <- colnames(myLoad$beta)

###############################################
# Combat correction
###############################################
pd <- myLoad$pd
pd$Status <- as.factor(pd$Status)
pd$Slide <- as.factor(pd$Slide)
pd$Array <- as.factor(pd$Array)
myNorm <- champ.runCombat(
  beta = myNorm,
  pd = pd,
  variablename = c("Age", "Sex", "Status"),
  batchname = c("Slide", "Array"),
  logitTrans = TRUE
)

###############################################
# Save
###############################################
write.csv(myNorm, file = "betas.csv")
