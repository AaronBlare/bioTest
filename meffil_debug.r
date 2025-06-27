rm(list=ls())

library(meffil)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(readxl)
library(stringr)

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

group <- pheno$Special.Status
age <- data.frame(Age = pheno$Age)
rownames(age) <- rownames(pheno)

beta.nodup <- meffil.collapse.dups(data.matrix(betas))

set.seed(1337)
ewas.ret.cont <- meffil.ewas(beta.nodup, variable=group, covariates=age, isva=F) 