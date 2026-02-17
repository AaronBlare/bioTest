rm(list=ls())

library(devtools)
devtools::install_github("YuanTian1991/ChAMP")
devtools::install_github("YuanTian1991/ChAMPData")

library("ChAMP")
library("methylGSA")
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(readxl)
library(stringr)

path <- "D:/Yandex.Disk/DNAm draft/Lesnoy_CVD/120_1"
setwd(path)

pheno <- read_excel("pheno_funnorm.xlsx")
pheno <- as.data.frame(pheno)
names(pheno) <- str_replace_all(names(pheno), c(" " = ".", "," = ""))
pheno$Special.Status <- as.factor(pheno$Special.Status)
colnames(pheno)[colnames(pheno) == '...1'] <- 'ID'
rownames(pheno) <- pheno[,1]
pheno <- pheno[,c("Age","Sex","Special.Status")]

betas <- read.csv("betas_funnorm.csv")
rownames(betas) <- betas[,1]
betas[,1] <- NULL
colnames(betas) <- gsub("^X", "", colnames(betas))

cpgs <- read.csv("GSEA(ebayes)_group_orgn_limma.csv")
#cpgs <- cpgs[cpgs$adj.P.Val <= 0.05, ]

genes <- read.csv("GSEA(ebayes)_group_genes_orgn_limma.csv")

manifest <- read.csv("EPIC-8v2-0_A2.csv")
manifest <- manifest[, c("IlmnID", "UCSC_RefGene_Name")]

####################################################################
### ChAMP GSEA function test
####################################################################

gsea <- champ.GSEA(beta=betas,
 DMP=NULL,
 DMR=NULL,
 CpGlist=cpgs$X,
 Genelist=NULL,
 pheno=pheno$Special.Status,
 method="gometh",
 arraytype="EPIC_V2",
 Rplot=TRUE,
 adjPval=0.05,
 cores=8)
capture.output(print(PathwayList), file = "GSEA(gometh)_cpg_orgn_champ.txt")

gsea <- champ.GSEA(beta=betas,
 DMP=NULL,
 DMR=NULL,
 CpGlist=NULL,
 Genelist=genes$loi.lv.CpG,
 pheno=pheno$Special.Status,
 method="gometh",
 arraytype="EPICv2",
 Rplot=TRUE,
 adjPval=0.05,
 cores=8)

####################################################################
### methylglm function test
####################################################################

cpg_pval <- setNames(cpgs$adj.P.Val, cpgs$X)
anno <- prepareAnnot(manifest)

GSEA_methylglm_go <- methylglm(
  cpg.pval = cpg_pval,
  array.type = NULL,
  FullAnnot = anno,
  group = "all",
  GS.idtype = "SYMBOL",
  GS.type = "GO",
  minsize = 5,
  maxsize = 10000,
  parallel = FALSE
)
write.csv(GSEA_methylglm_go, file = "GSEA(methylglm)_GO_orgn_5_10000.csv", row.names=FALSE)

GSEA_methylglm_kegg <- methylglm(
  cpg.pval = cpg_pval,
  array.type = NULL,
  FullAnnot = anno,
  group = "all",
  GS.idtype = "SYMBOL",
  GS.type = "KEGG",
  minsize = 5,
  maxsize = 10000,
  parallel = FALSE
)
write.csv(GSEA_methylglm_kegg, file = "GSEA(methylglm)_KEGG_orgn_5_10000.csv", row.names=FALSE)

GSEA_methylglm_react <- methylglm(
  cpg.pval = cpg_pval,
  array.type = NULL,
  FullAnnot = anno,
  group = "all",
  GS.idtype = "SYMBOL",
  GS.type = "Reactome",
  minsize = 5,
  maxsize = 10000,
  parallel = FALSE
)
write.csv(GSEA_methylglm_react, file = "GSEA(methylglm)_Reactome_orgn_5_10000.csv", row.names=FALSE)

GSEA_methylrra_go <- methylRRA(
    cpg.pval = cpg_pval,
    array.type = NULL,
    FullAnnot = anno,
    group = "all",
    GS.idtype = "SYMBOL",
    GS.type = "GO",
    minsize = 5,
    maxsize = 1000
)
write.csv(GSEA_methylrra_go, file = "GSEA(methylrra)_GO_orgn_5_1000.csv", row.names=FALSE)

GSEA_methylrra_kegg <- methylRRA(
    cpg.pval = cpg_pval,
    array.type = NULL,
    FullAnnot = anno,
    group = "all",
    GS.idtype = "SYMBOL",
    GS.type = "KEGG",
    minsize = 5,
    maxsize = 1000
)
write.csv(GSEA_methylrra_kegg, file = "GSEA(methylrra)_KEGG_orgn_5_1000.csv", row.names=FALSE)

GSEA_methylrra_react <- methylRRA(
    cpg.pval = cpg_pval,
    array.type = NULL,
    FullAnnot = anno,
    group = "all",
    GS.idtype = "SYMBOL",
    GS.type = "Reactome",
    minsize = 5,
    maxsize = 1000
)
write.csv(GSEA_methylrra_react, file = "GSEA(methylrra)_Reactome_orgn_5_1000.csv", row.names=FALSE)

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
  array.type = "EPIC_V2")

gsea_res_adj <- gsea_res_all[gsea_res_all$FDR<=0.05,]
if (!all(is.na(gsea_res_adj))) {
  write.csv(gsea_res_adj, file = "GSEA(gometh)_group_orgn.csv", row.names=TRUE)
}

####################################################################
### clusterProfiler testing
####################################################################

devtools::install_github(repo = "YuLab-SMU/clusterProfiler", host = "https://api.github.com")
library(clusterProfiler)

enrich_res_go <- enrichGO(
    genes$loi.lv.CpG,
    'org.Hs.eg.db',
    keyType = "SYMBOL",
    ont = "ALL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 5,
    maxGSSize = 1000,
)

gene_uniprot <- bitr(genes$loi.lv.CpG, fromType="SYMBOL", toType="UNIPROT", OrgDb="org.Hs.eg.db")
enrich_res_kegg <- enrichKEGG(
    gene_uniprot$UNIPROT,
    organism = "hsa",
    keyType = "UNIPROT",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 1000,
    qvalueCutoff = 0.2
)
enrich_res_mkegg <- enrichMKEGG(
    gene_uniprot$UNIPROT,
    organism = "hsa",
    keyType = "UNIPROT",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 1000,
    qvalueCutoff = 0.2
)

enrich_res_pc <- enrichPC(genes$loi.lv.CpG)

gene_entrez <- bitr(genes$loi.lv.CpG, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
enrich_res_wp <- enrichWP(gene_entrez$ENTREZID, "Homo sapiens")

####################################################################
### dmGSEA test
####################################################################

library(dmGsea)
library(data.table)
library(org.Hs.eg.db)

colnames(cpgs)[colnames(cpgs) == 'X'] <- 'Name'
colnames(cpgs)[colnames(cpgs) == 'adj.P.Val'] <- 'p'
cpgs <- cpgs[,c("Name","p")]

all_symbols <- unique(unlist(strsplit(na.omit(manifest$UCSC_RefGene_Name), ";")))
all_symbols <- unique(trimws(all_symbols))
symbol_to_entrez <- mapIds(org.Hs.eg.db,
                          keys = all_symbols,
                          column = "ENTREZID",
                          keytype = "SYMBOL",
                          multiVals = "first")

convert_symbols_fast <- function(symbol_string) {
  if (is.na(symbol_string) || symbol_string == "") return("")
  
  symbols <- trimws(unlist(strsplit(symbol_string, ";")))
  entrez_ids <- na.omit(unique(symbol_to_entrez[symbols]))
  
  if (length(entrez_ids) == 0) return(NA)
  paste(entrez_ids, collapse = ";")
}

if (is.data.table(manifest)) {
  manifest[, entrezid := sapply(UCSC_RefGene_Name, convert_symbols_fast)]
} else {
  manifest$entrezid <- sapply(manifest$UCSC_RefGene_Name, convert_symbols_fast)
}

colnames(manifest)[colnames(manifest) == 'IlmnID'] <- 'Name'
manifest_entrez <- manifest[,c("Name","entrezid")]

gsProbe(
  cpgs,
  FDRthre=0.05,
  GeneProbeTable=manifest_entrez,
  arrayType=NULL,
  gSetName='GO',
  species="Human",
  outfile="gsProbe_go",
  ncore=1)

gsGene(
  cpgs,
  method="Threshold",
  FDRthre=0.05,
  GeneProbeTable=manifest_entrez,
  gSetName="GO",
  species="Human",
  outfile="gsGene_go",
  ncore=1)

gsGene(
  cpgs,
  method="Ranking",
  GeneProbeTable=manifest_entrez,
  gSetName="GO",
  species="Human",
  outfile="gsGene_go_rank",
  ncore=1)

gsProbe(
  cpgs,
  FDRthre=0.05,
  GeneProbeTable=manifest_entrez,
  arrayType=NULL,
  gSetName='KEGG',
  species="Human",
  outfile="gsProbe_kegg",
  ncore=1)

gsGene(
  cpgs,
  method="Threshold",
  FDRthre=0.05,
  GeneProbeTable=manifest_entrez,
  gSetName="KEGG",
  species="Human",
  outfile="gsGene_kegg",
  ncore=1)

gsGene(
  cpgs,
  method="Ranking",
  GeneProbeTable=manifest_entrez,
  gSetName="KEGG",
  species="Human",
  outfile="gsGene_kegg_rank",
  ncore=1)

gsProbe(
  cpgs,
  FDRthre=0.05,
  GeneProbeTable=manifest_entrez,
  arrayType=NULL,
  gSetName='MSigDB',
  species="Human",
  outfile="gsProbe_msig",
  ncore=1)

gsGene(
  cpgs,
  method="Threshold",
  FDRthre=0.05,
  GeneProbeTable=manifest_entrez,
  gSetName="MSigDB",
  species="Human",
  outfile="gsGene_msig",
  ncore=1)

gsGene(
  cpgs,
  method="Ranking",
  GeneProbeTable=manifest_entrez,
  gSetName="MSigDB",
  species="Human",
  outfile="gsGene_msig_rank",
  ncore=1)