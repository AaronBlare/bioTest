rm(list=ls())

###############################################
# Installing required packages
###############################################
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ChAMP")
BiocManager::install("methylGSA")
BiocManager::install("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
BiocManager::install("limma")
BiocManager::install("missMethyl")

library("ChAMP")
library("methylGSA")
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(readxl)
library(stringr)
library(limma)
library(missMethyl)
library(openxlsx)

###############################################
# Setting path
# Directory must contain file with beta-values and *.xlsx file with phenotype
###############################################
path <- "path/to/directory/with/files"
setwd(path)

###############################################
# Read data
###############################################
pheno <- read_excel("data.xlsx")
pheno <- as.data.frame(pheno)
names(pheno) <- str_replace_all(names(pheno), c(" " = ".", "," = ""))
pheno$Special.Status <- as.factor(pheno$Special.Status)
colnames(pheno)[colnames(pheno) == '...1'] <- 'ID'
rownames(pheno) <- pheno[,1]
pheno <- pheno[,c("Age","Sex","Special.Status")]

betas <- read.csv("betas.csv")
rownames(betas) <- betas[,1]
betas[,1] <- NULL
colnames(betas) <- gsub("^X", "", colnames(betas))
betas <- betas[grepl("^cg", rownames(betas)),]

###############################################
# Get limma results for GSEA
###############################################
group <- factor(pheno$Special.Status, levels=c("Control","Case"))
design <- model.matrix(~group)
row.names(design) <- row.names(pheno)

fit <- lmFit(betas, design)
fit <- eBayes(fit, proportion=0.01, robust=TRUE)
top <- topTable(fit, adjust="BH", sort.by="B", number=nrow(fit))
write.xlsx(top, file = "limma.xlsx", row.names=TRUE)

###############################################
# Get required data and annotations for EWAS
###############################################
cpgs <- read.xlsx("limma.xlsx")
manifest <- read.csv("path/to/manifest/file")
manifest <- manifest[, c("IlmnID", "UCSC_RefGene_Name")]

cpg_pval <- setNames(cpgs$adj.P.Val, cpgs$X)
anno <- prepareAnnot(manifest)

###############################################
# Get methylglm results for GSEA
###############################################
GSEA_methylglm_go <- methylglm(
  cpg.pval = cpg_pval,
  array.type = NULL,
  FullAnnot = anno,
  group = "all",
  GS.idtype = "SYMBOL",
  GS.type = "GO",
  minsize = 5,
  maxsize = 1000,
  parallel = FALSE
)
write.csv(GSEA_methylglm_go, file = "GSEA_methylglm_GO.csv", row.names=FALSE)

GSEA_methylglm_kegg <- methylglm(
  cpg.pval = cpg_pval,
  array.type = NULL,
  FullAnnot = anno,
  group = "all",
  GS.idtype = "SYMBOL",
  GS.type = "KEGG",
  minsize = 5,
  maxsize = 1000,
  parallel = FALSE
)
write.csv(GSEA_methylglm_kegg, file = "GSEA_methylglm_KEGG.csv", row.names=FALSE)

GSEA_methylglm_react <- methylglm(
  cpg.pval = cpg_pval,
  array.type = NULL,
  FullAnnot = anno,
  group = "all",
  GS.idtype = "SYMBOL",
  GS.type = "Reactome",
  minsize = 5,
  maxsize = 1000,
  parallel = FALSE
)
write.csv(GSEA_methylglm_react, file = "GSEA_methylglm_Reactome.csv", row.names=FALSE)
