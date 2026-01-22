rm(list=ls())

library(devtools)
devtools::install_github("YuanTian1991/ChAMP")
devtools::install_github("YuanTian1991/ChAMPData")

library("ChAMP")
library("methylGSA")
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(readxl)
library(stringr)

path <- "D:/Yandex.Disk/DNAm draft/Lesnoy_CVD/GSE220622/generalized"
setwd(path)

pheno <- read_excel("pheno_funnorm.xlsx")
pheno <- as.data.frame(pheno)
names(pheno) <- str_replace_all(names(pheno), c(" " = ".", "," = ""))
pheno$Status <- as.factor(pheno$Status)
colnames(pheno)[colnames(pheno) == '...1'] <- 'ID'
rownames(pheno) <- pheno[,1]
pheno <- pheno[,c("Age","Sex","Status")]

betas <- read.csv("betas_funnorm.csv")
rownames(betas) <- betas[,1]
betas[,1] <- NULL
colnames(betas) <- gsub("^X", "", colnames(betas))

cpgs <- read.csv("GSEA(ebayes)_group_orgn_limma.csv")
#cpgs <- cpgs[cpgs$adj.P.Val <= 0.05, ]

genes <- read.csv("GSEA(ebayes)_group_genes_orgn_limma.csv")

manifest <- read.csv("infinium-methylationepic-v-1-0-b5-manifest-file.csv")
manifest <- manifest[, c("IlmnID", "UCSC_RefGene_Name")]

####################################################################
### methylglm function test
####################################################################

cpg_pval <- setNames(cpgs$adj.P.Val, cpgs$X)
anno <- prepareAnnot(manifest)

GSEA_methylglm_go <- methylglm(
  cpg.pval = cpg_pval,
  array.type = 'EPIC',
  group = "all",
  GS.idtype = "SYMBOL",
  GS.type = "GO",
  minsize = 100,
  maxsize = 500,
  parallel = FALSE
)
write.csv(GSEA_methylglm_go, file = "GSEA(methylglm)_GO_orgn_100_500.csv", row.names=FALSE)

GSEA_methylglm_kegg <- methylglm(
  cpg.pval = cpg_pval,
  array.type = 'EPIC',
  group = "all",
  GS.idtype = "SYMBOL",
  GS.type = "KEGG",
  minsize = 100,
  maxsize = 500,
  parallel = FALSE
)
write.csv(GSEA_methylglm_kegg, file = "GSEA(methylglm)_KEGG_orgn_100_500.csv", row.names=FALSE)

GSEA_methylglm_react <- methylglm(
  cpg.pval = cpg_pval,
  array.type = 'EPIC',
  group = "all",
  GS.idtype = "SYMBOL",
  GS.type = "Reactome",
  minsize = 100,
  maxsize = 500,
  parallel = FALSE
)
write.csv(GSEA_methylglm_react, file = "GSEA(methylglm)_Reactome_orgn_100_500.csv", row.names=FALSE)

GSEA_methylrra_go <- methylRRA(
    cpg.pval = cpg_pval,
    array.type = 'EPIC',
    group = "all",
    GS.idtype = "SYMBOL",
    GS.type = "GO",
    minsize = 100,
    maxsize = 500
)
write.csv(GSEA_methylrra_go, file = "GSEA(methylrra)_GO_orgn_100_500.csv", row.names=FALSE)

GSEA_methylrra_kegg <- methylRRA(
    cpg.pval = cpg_pval,
    array.type = 'EPIC',
    group = "all",
    GS.idtype = "SYMBOL",
    GS.type = "KEGG",
    minsize = 100,
    maxsize = 500
)
write.csv(GSEA_methylrra_kegg, file = "GSEA(methylrra)_KEGG_orgn_100_500.csv", row.names=FALSE)

GSEA_methylrra_react <- methylRRA(
    cpg.pval = cpg_pval,
    array.type = 'EPIC',
    group = "all",
    GS.idtype = "SYMBOL",
    GS.type = "Reactome",
    minsize = 100,
    maxsize = 500
)
write.csv(GSEA_methylrra_react, file = "GSEA(methylrra)_Reactome_orgn_100_500.csv", row.names=FALSE)

####################################################################
### limma gometh function test
####################################################################

library(limma)
library(missMethyl)

cpgs_orgn <- read.csv("cpgs_orgn.csv")
cpgs_orgn <- as.character(cpgs_orgn[,1])

cpgs_top <- cpgs[cpgs$adj.P.Val <= 0.05, ]
gsea_res_all <- gometh(
  cpgs_top$X,
  all.cpg = cpgs_orgn,
  collection = c("GO", "KEGG"),
  array.type = "EPIC")

gsea_res_adj <- gsea_res_all[gsea_res_all$FDR<=0.05,]
if (!all(is.na(gsea_res_adj))) {
  write.csv(gsea_res_adj, file = "GSEA(gometh)_group_orgn.csv", row.names=TRUE)
}

####################################################################
### dmGSEA test
####################################################################

library(dmGsea)
library(data.table)

colnames(cpgs)[colnames(cpgs) == 'X'] <- 'Name'
colnames(cpgs)[colnames(cpgs) == 'adj.P.Val'] <- 'p'
cpgs <- cpgs[,c("Name","p")]

gsProbe(
  cpgs,
  FDRthre=0.05,
  arrayType='EPIC',
  gSetName='GO',
  species="Human",
  outfile="gsProbe_go",
  ncore=1)

gsGene(
  cpgs,
  method="Threshold",
  FDRthre=0.05,
  arrayType='EPIC',
  gSetName="GO",
  species="Human",
  outfile="gsGene_go",
  ncore=1)

gsGene(
  cpgs,
  method="Ranking",
  arrayType='EPIC',
  gSetName="GO",
  species="Human",
  outfile="gsGene_go_rank",
  ncore=1)

gsProbe(
  cpgs,
  FDRthre=0.05,
  arrayType='EPIC',
  gSetName='KEGG',
  species="Human",
  outfile="gsProbe_kegg",
  ncore=1)

gsGene(
  cpgs,
  method="Threshold",
  FDRthre=0.05,
  arrayType='EPIC',
  gSetName="KEGG",
  species="Human",
  outfile="gsGene_kegg",
  ncore=1)

gsGene(
  cpgs,
  method="Ranking",
  arrayType='EPIC',
  gSetName="KEGG",
  species="Human",
  outfile="gsGene_kegg_rank",
  ncore=1)

gsProbe(
  cpgs,
  FDRthre=0.05,
  arrayType='EPIC',
  gSetName='MSigDB',
  species="Human",
  outfile="gsProbe_msig",
  ncore=1)

gsGene(
  cpgs,
  method="Threshold",
  FDRthre=0.05,
  arrayType='EPIC',
  gSetName="MSigDB",
  species="Human",
  outfile="gsGene_msig",
  ncore=1)

gsGene(
  cpgs,
  method="Ranking",
  arrayType='EPIC',
  gSetName="MSigDB",
  species="Human",
  outfile="gsGene_msig_rank",
  ncore=1)
