rm(list=ls())

library("ChAMP")
library("methylGSA")
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(readxl)
library(stringr)
library(limma)
library(DMRcate)

path <- "E:/YandexDisk/bbd/fmba/dnam/processed/special_63"
setwd(path)

pheno <- read_excel("pheno.xlsx")
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

group <- factor(pheno$Special.Status, levels=c("Control","Case"))
age <- pheno$Age
design <- model.matrix(~group)
row.names(design) <- row.names(pheno)

annotation_diff_var <- cpg.annotate(
  datatype="array", 
  object=na.omit(data.matrix(betas)), what="Beta", 
  arraytype="EPICv2", 
  analysis.type="diffVar", design=design, 
  contrasts = FALSE, cont.matrix = NULL, 
  fdr=0.05, varFitcoef=2) 
diff_var_DMRs <- dmrcate(annotation_diff_var, lambda=1000, C=2)
if (diff_var_DMRs) {
  results.ranges.diff.var <- extractRanges(diff_var_DMRs)
  results.ranges.diff.var.sign <- results.ranges.diff.var[results.ranges.diff.var$HMFDR<=0.05,]
  if (!all(is.na(results.ranges.diff.var.sign))) {
    write.csv(results.ranges.diff.var.sign, file = "DMRcate_diff_var_group_orgn.csv", row.names=TRUE)
  }
}