rm(list=ls())

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DSS")
BiocManager::install("minfi")
BiocManager::install("ChAMPdata")
BiocManager::install("wateRmelon")
BiocManager::install("ChAMP")
BiocManager::install("methylGSA")
BiocManager::install("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
install.packages("reticulate")
install.packages("readxl")
install.packages("stringr")

library(devtools)
devtools::install_github("YuanTian1991/ChAMP")
devtools::install_github("YuanTian1991/ChAMPData")

library("ChAMP")
library("methylGSA")
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(readxl)
library(stringr)

### Python test
install.packages("reticulate")
Sys.setenv(RETICULATE_PYTHON = "C:/Users/alena/anaconda3/envs/py312/python.exe")
library("reticulate")
py_config()
Sys.which('python')
use_condaenv('py312')
py_run_string('print(1+1)')
rm(list=ls())
pd <- import("pandas")
path_load <- "E:/YandexDisk/pydnameth/datasets/GPL21145/GSEUNN/special/043_yakutia_EWAS/00_all_region/data_for_R"
pheno_py <- pd$read_pickle(paste(path_load, "/pheno_R_all_region.pkl", sep=''))
pheno_py$Region <- as.factor(pheno_py$Region)
betas_py <- pd$read_pickle(paste(path_load, "/betas_R_all_region.pkl", sep=''))
###

dmp_pval <- 1
dmr_pval <- 0.05
dmr_min_probes <- 10
gsea_pval <- 0.05
methylglm_minsize <- 10
methylglm_maxsize <- 1000

path <- "E:/YandexDisk/bbd/fmba/dnam/processed/special_63/swan"
setwd(path)

pheno <- read_excel("pheno_swan.xlsx")
pheno <- as.data.frame(pheno)
names(pheno) <- str_replace_all(names(pheno), c(" " = ".", "," = ""))
pheno$Special.Status <- as.factor(pheno$Special.Status)
colnames(pheno)[colnames(pheno) == '...1'] <- 'ID'
rownames(pheno) <- pheno[,1]
pheno <- pheno[,c("Age","Sex","Special.Status")]

betas <- read.csv("betas_swan.csv")
rownames(betas) <- betas[,1]
betas[,1] <- NULL
colnames(betas) <- gsub("^X", "", colnames(betas))

### DMP function test
dmp <- champ.DMP(
  beta = betas,
  pheno = pheno$Special.Status,
  compare.group = c("Control", "Case"),
  adjPVal = dmp_pval,
  adjust.method = "BH",
  arraytype = "EPICv2"
)
write.csv(dmp$Control_to_Case, file = "DMP_orgn_champ.csv")

cpgs_fltr <- read.csv("cpgs_fltd.csv")
cpgs_fltr <- as.character(cpgs_fltr[,1])
betas_fltr <- betas[row.names(betas) %in% cpgs_fltr,]
dmp_fltr <- champ.DMP(
  beta = betas_fltr,
  pheno = pheno$Special.Status,
  compare.group = c("Control", "Case"),
  adjPVal = dmp_pval,
  adjust.method = "BH",
  arraytype = "EPICv2"
)
write.csv(dmp_fltr$Control_to_Case, file = "DMP_fltr_champ.csv")
###
### DMR function test
dmr <- champ.DMR(
  beta = data.matrix(betas),
  pheno = pheno$Special.Status,
  compare.group = c("Control", "Case"),
  arraytype = "EPICv2",
  method = "Bumphunter", # "Bumphunter" "ProbeLasso" "DMRcate"
  minProbes = dmr_min_probes,
  adjPvalDmr = dmr_pval,
  cores = 8,
  ## following parameters are specifically for Bumphunter method.
  maxGap = 300,
  cutoff = NULL,
  pickCutoff = TRUE,
  smooth = TRUE,
  smoothFunction = loessByCluster,
  useWeights = FALSE,
  permutations = NULL,
  B = 250,
  nullMethod = "bootstrap"
)
write.csv(dmr$BumphunterDMR, file = "DMR_orgn_champ.csv")

dmr_fltr <- champ.DMR(
  beta = data.matrix(betas_fltr),
  pheno = pheno$Special.Status,
  compare.group = c("Control", "Case"),
  arraytype = "EPICv2",
  method = "Bumphunter", # "Bumphunter" "ProbeLasso" "DMRcate"
  minProbes = dmr_min_probes,
  adjPvalDmr = dmr_pval,
  cores = 8,
  ## following parameters are specifically for Bumphunter method.
  maxGap = 300,
  cutoff = NULL,
  pickCutoff = TRUE,
  smooth = TRUE,
  smoothFunction = loessByCluster,
  useWeights = FALSE,
  permutations = NULL,
  B = 250,
  nullMethod = "bootstrap"
)
write.csv(dmr_fltr$BumphunterDMR, file = "DMR_fltr_champ.csv")

RSobject <- RatioSet(betas, annotation = c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38"))
RSanno <- getAnnotation(RSobject)[, c("chr", "pos", "Name", "UCSC_RefGene_Name")]
loi.lv <- list()
cpg.idx <- unique(unlist(apply(dmr[[1]], 1, function(x) rownames(RSanno)[which(RSanno$chr == x[1] & RSanno$pos >= as.numeric(x[2]) & RSanno$pos <= as.numeric(x[3]))])))
loi.lv[["DMR"]] <- unique(unlist(sapply(RSanno[cpg.idx, "UCSC_RefGene_Name"], function(x) strsplit(x, split = ";")[[1]])))
write.csv(data.frame(loi.lv$DMR), file = "DMR_genes_orgn_champ.csv", row.names=FALSE)

RSobject <- RatioSet(betas_fltr, annotation = c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38"))
RSanno <- getAnnotation(RSobject)[, c("chr", "pos", "Name", "UCSC_RefGene_Name")]
loi.lv <- list()
cpg.idx <- unique(unlist(apply(dmr[[1]], 1, function(x) rownames(RSanno)[which(RSanno$chr == x[1] & RSanno$pos >= as.numeric(x[2]) & RSanno$pos <= as.numeric(x[3]))])))
loi.lv[["DMR"]] <- unique(unlist(sapply(RSanno[cpg.idx, "UCSC_RefGene_Name"], function(x) strsplit(x, split = ";")[[1]])))
write.csv(data.frame(loi.lv$DMR), file = "DMR_genes_fltr_champ.csv", row.names=FALSE)

gsea <- champ.GSEA(beta=betas, 
  DMP=dmp[[1]],
  DMR=dmr,
  CpGlist=NULL,
  Genelist=NULL,
  pheno=pheno$Special.Status,
  method="fisher",
  arraytype="EPICv2",
  Rplot=TRUE,
  adjPval=gsea_pval,
  cores=8)
if (!all(is.na(gsea$DMP))) {
  gtResult_DMP_orgn <- data.frame(row.names(gsea$DMP), gsea$DMP)
  write.csv(gtResult_DMP_orgn, file = "GSEA(fisher)_DMP_orgn_champ.csv", row.names=TRUE)
}
if (!all(is.na(gsea$DMR))) {
  gtResult_DMR_orgn <- data.frame(row.names(gsea$DMR), gsea$DMR)
  write.csv(gtResult_DMR_orgn, file = "GSEA(fisher)_DMR_orgn_champ.csv", row.names=TRUE)
}

gsea_fltr <- champ.GSEA(beta=betas_fltr, 
  DMP=dmp_fltr[[1]],
  DMR=dmr_fltr,
  CpGlist=NULL,
  Genelist=NULL,
  pheno=pheno$Special.Status,
  method="fisher",
  arraytype="EPICv2",
  Rplot=TRUE,
  adjPval=gsea_pval,
  cores=8)
if (!all(is.na(gsea_fltr$DMP))) {
  gtResult_DMP_fltr <- data.frame(row.names(gsea_fltr$DMP), gsea_fltr$DMP)
  write.csv(gtResult_DMP_fltr, file = "GSEA(fisher)_DMP_fltr_champ.csv", row.names=TRUE)
}
if (!all(is.na(gsea_fltr$DMR))) {
  gtResult_DMR_fltr <- data.frame(row.names(gsea_fltr$DMR), gsea_fltr$DMR)
  write.csv(gtResult_DMR_fltr, file = "GSEA(fisher)_DMR_fltr_champ.csv", row.names=TRUE)
}
###

### Seems that GSEA function does not work with EPICv2
gsea <- champ.GSEA(
  beta = betas,
  DMP = NULL,
  DMR = NULL,
  CpGlist = NULL,
  Genelist = NULL,
  pheno = pheno$Special.Status,
  method = "ebayes",
  arraytype = "EPICv2",
  Rplot = FALSE,
  adjPval = gsea_pval,
  cores = 8
)
gtResult <- data.frame(row.names(gsea[[3]]), gsea[[3]])
colnames(gtResult)[1] <- "ID"
write.csv(gtResult, file = "GSEA(ebayes)_gtResult_orgn.csv", row.names=TRUE)
write.csv(gsea$GSEA[[1]], file = "GSEA(ebayes)_Rank(P)_orgn.csv", row.names=TRUE)
###
### Seems that methylglm function does not work with EPICv2
dmp_df <- data.frame(row.names(dmp$Control_to_Case), dmp$Control_to_Case)
colnames(dmp_df)[1] <- "CpG"
cpg_pval <- setNames(dmp_df$adj.P.Val, dmp_df$CpG)
GSEA_methylglm <- methylglm(
  cpg.pval = cpg_pval,
  array.type = "EPICv2",
  group = "all",
  GS.idtype = "SYMBOL",
  GS.type = "GO",
  minsize = methylglm_minsize,
  maxsize = methylglm_maxsize,
  parallel = TRUE
)
write.csv(GSEA_methylglm, file = "GSEA(methylglm)_GO_orgn.csv", row.names=FALSE)
GSEA_methylglm <- methylglm(
  cpg.pval = cpg_pval,
  array.type = "EPIC",
  group = "all",
  GS.idtype = "SYMBOL",
  GS.type = "KEGG",
  minsize = methylglm_minsize,
  maxsize = methylglm_maxsize,
  parallel = TRUE
)
write.csv(GSEA_methylglm, file = "GSEA(methylglm)_KEGG_orgn.csv", row.names=FALSE)
GSEA_methylglm <- methylglm(
  cpg.pval = cpg_pval,
  array.type = "EPIC",
  group = "all",
  GS.idtype = "SYMBOL",
  GS.type = "Reactome",
  minsize = methylglm_minsize,
  maxsize = methylglm_maxsize,
  parallel = TRUE
)
write.csv(GSEA_methylglm, file = "GSEA(methylglm)_Reactome_orgn.csv", row.names=FALSE)
###
### Seems that GSEA function with gometh method does not work with EPICv2
GSEA_gometh <- champ.GSEA(
  beta = betas,
  DMP = NULL,
  DMR = dmr,
  CpGlist = NULL,
  Genelist = NULL,
  pheno = pheno$Special.Status,
  method = "gometh",
  arraytype = "EPICv2",
  Rplot = TRUE,
  adjPval = dmr_pval,
  cores = 8
)
write.csv(data.frame(GSEA_gometh$DMR), file = "DMR_GSEA_gometh_orgn.csv", row.names=FALSE)
###