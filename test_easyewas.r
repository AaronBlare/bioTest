rm(list=ls())

library(easyEWAS)
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
rownames(pheno) <- as.character(rownames(pheno))
pheno[,1] <- as.character(pheno[,1])
pheno <- pheno[,c("ID", "Age","Sex","Special.Status")]

betas <- read.csv("betas.csv")
rownames(betas) <- betas[,1]
colnames(betas) <- gsub("^X", "", colnames(betas))
colnames(betas) <- as.character(colnames(betas))

# prepare the data file ------
res <- initEWAS(outpath = path)
res <- loadEWAS(input = res,
                ExpoData = pheno,
                MethyData = betas)