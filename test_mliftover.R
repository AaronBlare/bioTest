rm(list=ls())

###############################################
# Installing packages
###############################################
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
install.packages("devtools")
devtools::install_url("https://cran.r-project.org/src/contrib/Archive/kpmt/kpmt_0.1.0.tar.gz")
BiocManager::install("ChAMP")
library("ChAMP")

###############################################
# Setting variables
###############################################
dataset_450 <- 'GSE87571'
arraytype_450 <- '450K'

dataset_EPIC <- 'GSE234461'
arraytype_EPIC <- 'EPIC'

###############################################
# Setting path
# Directory must contain *.idat files and *.csv file with phenotype
###############################################
path_data_450 <- "E:/YandexDisk/DNAm draft/GEO/GPL13534/24_GSE87571/raw"
path_data_EPIC <- "E:/YandexDisk/DNAm draft/GEO/GPL21145/93_GSE234461/raw"

path_work <- path_data_450
setwd(path_work)

###############################################
# Import and filtration
###############################################
myLoad_450 <- champ.load(
  directory = path_data_450,
  arraytype = arraytype_450,
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
pd_450 <- as.data.frame(myLoad_450$pd)

###############################################
# Functional normalization
###############################################
betas_450 <- getBeta(preprocessFunnorm(myLoad_450$rgSet))

###############################################
path_work <- path_data_EPIC
setwd(path_work)

###############################################
# Import and filtration
###############################################
myLoad_EPIC <- champ.load(
  directory = path_data_EPIC,
  arraytype = arraytype_EPIC,
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
pd_EPIC <-  as.data.frame(myLoad_EPIC$pd)

###############################################
# Functional normalization
###############################################
betas_EPIC <- getBeta(preprocessFunnorm(myLoad_EPIC$rgSet))

betas_450k_to_epic = mLiftOver(betas_450, "EPIC", impute=FALSE)
length(betas_450k_to_epic)