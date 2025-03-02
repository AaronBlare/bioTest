rm(list=ls())

###############################################
# Installing packages
###############################################
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ChAMP")
BiocManager::install("methylGSA")
BiocManager::install("preprocessCore")
library("ChAMP")
library("preprocessCore")
library(readxl)

###############################################
# Setting variables
###############################################
arraytype <- '450K'
dataset <- 'GSE87571'

###############################################
# Setting path
###############################################
path_data <- "E:/YandexDisk/pydnameth/datasets/GPL13534/GSE87571/raw/idat"
path_save <- "D:/bioTest/attack/data"
setwd(path_data)

###############################################
# Import and filtration
###############################################
myLoad <- champ.load(
  directory = path_data,
  arraytype = arraytype,
  method = "minfi",
  methValue = "B",
  autoimpute = TRUE,
  filterDetP = TRUE,
  ProbeCutoff = 0.1,
  SampleCutoff = 0.1,
  detPcut = 0.01,
  filterBeads = FALSE,
  beadCutoff = 0.05,
  filterNoCG = FALSE,
  filterSNPs = FALSE,
  filterMultiHit = FALSE,
  filterXY = FALSE,
  force = TRUE
)
pd <- as.data.frame(myLoad$pd)

myNormBMIQ <- champ.norm(
    beta=myLoad$beta,
    rgSet=myLoad$rgSet,
    mset=myLoad$mset,
    resultsDir=path_save,
    method="BMIQ",
    plotBMIQ=FALSE,
    arraytype=arraytype,
    cores=3)
myNormFun <- champ.norm(
    beta=myLoad$beta,
    rgSet=myLoad$rgSet,
    mset=myLoad$mset,
    resultsDir=path_save,
    method="FunctionalNormalization",
    plotBMIQ=FALSE,
    arraytype=arraytype,
    cores=3)
myNormSWAN <- champ.norm(
    beta=myLoad$beta,
    rgSet=myLoad$rgSet,
    mset=myLoad$mset,
    resultsDir=path_save,
    method="SWAN",
    plotBMIQ=FALSE,
    arraytype=arraytype,
    cores=3)
myNormPBC <- champ.norm(
    beta=myLoad$beta,
    rgSet=myLoad$rgSet,
    mset=myLoad$mset,
    resultsDir=path_save,
    method="PBC",
    plotBMIQ=FALSE,
    arraytype=arraytype,
    cores=3)

###############################################
# Save modified pheno
###############################################
write.csv(pd, file = "D:/bioTest/attack/data/pheno.csv")

###############################################
# Save DNAm data
###############################################
write.csv(myNormBMIQ, file = "D:/bioTest/attack/data/betasBMIQ.csv")
write.csv(myNormFun, file = "D:/bioTest/attack/data/betasFun.csv")
write.csv(myNormSWAN, file = "D:/bioTest/attack/data/betasSWAN.csv")
write.csv(myNormPBC, file = "D:/bioTest/attack/data/betasPBC.csv")
