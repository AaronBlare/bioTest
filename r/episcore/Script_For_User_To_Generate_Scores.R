rm(list = ls())

## FOR THE USER - PLACE INPUT FILES HERE - CAN CHANGE readRDS or read.csv if you so wish 

setwd("C:/Git/bioTest/r/episcore")

## Start to Process Files 
data=read.csv("data_Les_120.csv", row.names='cpg')
sexageinfo=read.csv("age_sex_Les_120.csv")
sexageinfo$ID = paste0('X', sexageinfo$ID)

message("1. Loading data") 

message("1.1 Loading Methylation data - rows to be CpGs and columns to be individuals") 

cpgs <- read.csv("./Predictors_Shiny_by_Groups.csv", header = T) 

## Read in separate files for Bernabeu age prediction as it requires two-step process 
coefficients <- read.delim("cage_coefficients_linear.tsv", sep = "\t")
coefficients_log <- read.delim("cage_coefficients_log.tsv", sep = "\t")
row.names(coefficients) <- coefficients[,"CpG_Site"]
row.names(coefficients_log) <- coefficients_log[,"CpG_Site"]

## Prepare CpGs for Bernabeu age predictor 
# Coefficients for log models
coef_log_2 <- coefficients_log[rownames(coefficients_log)[grep('_2', rownames(coefficients_log))],,drop=FALSE]
coef_log_2_simp <- gsub('_2', '', rownames(coef_log_2))
coef_log <- coefficients_log[which(!(rownames(coefficients_log) %in% rownames(coef_log_2))),,drop=FALSE]


# Coefficients for non-log models
coef_2 <- coefficients[rownames(coefficients)[grep('_2', rownames(coefficients))],,drop=FALSE]
coef_2_simp <- gsub('_2', '', rownames(coef_2))
coef1 <- coefficients[rownames(coefficients)[which(!(rownames(coefficients) %in% rownames(coef_2)))],,drop=FALSE]

# Beta means - for Bernabeu predictor (processed separately from other models)
means <- read.delim("cpg_meanbeta_gs20k.tsv")
row.names(means) <- means[,"cpg"]

# Total CpGs
cpgs_linear <- union(rownames(coef1), rownames(coef_log))
cpgs_squared <- union(coef_2_simp, coef_log_2_simp)
all_cpgs <- union(cpgs_linear, cpgs_squared)

## Check if Data needs to be Transposed

message("2. Quality Control and data Preparation") 

message("2.1 Checking if Row Names are CpG Sites") 

if(ncol(data) > nrow(data)){
  message("It seems that individuals are rows - data will be transposed!")
  data<-t(data) 
}

message("2.2 Subsetting CpG sites to those required for Predictor Calculation") 

## Subset CpG sites to those present on list for predictors 

coef=data[intersect(rownames(data), cpgs$CpG_Site),]

# Process separately for Bernabeu predictor 
coef.cage=data[intersect(rownames(data), all_cpgs),]


message("2.3 Checking if Beta Values are Present") 

## Conversion to Beta Values if M Values are likely present  

m_to_beta <- function (val) 
{
  beta <- 2^val/(2^val + 1)
  return(beta)
}

coef<-if((range(coef,na.rm=T)> 1)[[2]] == "TRUE") { message("Suspect that M Values are present. Converting to Beta Values");m_to_beta(coef) } else { message("Suspect that Beta Values are present");coef}
coef.cage<-if((range(coef.cage,na.rm=T)> 1)[[2]] == "TRUE") { message("Suspect that M Values are present. Converting to Beta Values");m_to_beta(coef.cage) } else { message("Suspect that Beta Values are present");coef.cage}


message("2.4 Convert NA Values to Mean for each Probe") 

## Convert NAs to Mean Value for all individuals across each probe 

na_to_mean <-function(methyl){
  methyl[is.na(methyl)]<-mean(methyl,na.rm=T)
  return(methyl)
}

coef <- t(apply(coef,1,function(x) na_to_mean(x)))
coef.cage <- t(apply(coef.cage,1,function(x) na_to_mean(x)))

message("2.5 Find CpGs not present in uploaded file, add these with mean Beta Value for CpG site from Training Sample") 

## Identify CpGs missing from input dataframe, include them and provide values as mean methylation value at that site

coef <- if(nrow(coef) == length(unique(cpgs$CpG_Site))) { message("All sites present"); coef } else if(nrow(coef)==0){ 
  message("There Are No Necessary CpGs in The dataset - All Individuals Would Have Same Values For Predictors. Analysis Is Not Informative!")
} else { 
  missing_cpgs = cpgs[-which(cpgs$CpG_Site %in% rownames(coef)),c("CpG_Site","Mean_Beta_Value")]
  message(paste(length(unique(missing_cpgs$CpG_Site)), "unique sites are missing - add to dataset with mean Beta Value from Training Sample", sep = " "))
  mat = matrix(nrow=length(unique(missing_cpgs$CpG_Site)),ncol = ncol(coef))
  row.names(mat) <- unique(missing_cpgs$CpG_Site)
  colnames(mat) <- colnames(coef) 
  mat[is.na(mat)] <- 1
  missing_cpgs1 <- if(length(which(duplicated(missing_cpgs$CpG_Site))) > 1) { 
    missing_cpgs[-which(duplicated(missing_cpgs$CpG_Site)),]
  } else {missing_cpgs
  }  
  ids = unique(row.names(mat))
  missing_cpgs1 = missing_cpgs1[match(ids,missing_cpgs1$CpG_Site),]
  mat=mat*missing_cpgs1$Mean_Beta_Value
  coef=rbind(coef,mat)
} 


coef.cage <- if (nrow(coef.cage)==length(all_cpgs)) { 
  message("All sites present")
  coef.cage
} else if (nrow(coef.cage)==0) { 
  message("There Are No Necessary CpGs in The dataset - All Individuals Would Have Same Values For Predictors. Analysis Is Not Informative!")
} else {
  missing_cpgs <- all_cpgs[which(!(all_cpgs %in% rownames(coef.cage)))]
  message(paste(length(missing_cpgs), "unique sites are missing - adding to dataset with mean Beta Value from GS training sample (N = 18,413)", sep = " "))
  mat = matrix(nrow=length(missing_cpgs), ncol = ncol(coef.cage))
  row.names(mat) <- missing_cpgs
  colnames(mat) <- colnames(coef.cage)
  mat[is.na(mat)] <- 1
  ids <- unique(row.names(mat))
  missing_cpg_means <- means[ids,]
  mat <- mat*missing_cpg_means[,"mean"]
  coef.cage <- rbind(coef.cage,mat)
} 


message("3. Calculating the Predictors") 
loop=unique(cpgs$Predictor)
out <- data.frame()
for(i in loop){ 
  tmp=coef[intersect(row.names(coef),cpgs[cpgs$Predictor %in% i,"CpG_Site"]),]
  tmp_coef = cpgs[cpgs$Predictor %in% i, ]
  if(nrow(tmp_coef) > 1) { 
    tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG_Site),]
    out[colnames(coef),i]=colSums(tmp_coef$Coefficient*tmp)
  } else {
    tmp2 = as.matrix(tmp)*tmp_coef$Coefficient 
    out[colnames(coef),i] = colSums(tmp2)
  }
} 
out$'Epigenetic Age (Zhang)' <- out$'Epigenetic Age (Zhang)' + 65.79295
out$ID <- row.names(out) 
out <- out[,c(ncol(out),1:(ncol(out)-1))] 


message("4. Obtaining cAge predictions - Bernabeu et al.") 
message("4.1. Preparing data for prediction") 

  scores <- coef.cage[rownames(coef1),]
  scores_quadratic <- coef.cage[coef_2_simp,]**2
  rownames(scores_quadratic) <- paste0(rownames(scores_quadratic), "_2")
  scores_linear <- rbind(scores, scores_quadratic)
  
  ## Prep for log predictor
  scores_log <- coef.cage[rownames(coef_log),]
  scores_quadratic_log <- coef.cage[coef_log_2_simp,]**2
  rownames(scores_quadratic_log) <- paste0(rownames(scores_quadratic_log), "_2")
  scores_log <- rbind(scores_log, scores_quadratic_log)
  
  ## Calculate cAge with linear model
  coefficients <- coefficients[rownames(scores_linear),]
  pred_linear <- scores_linear * coefficients[,"Coefficient"]
  pred_linear_pp <- colSums(pred_linear,na.rm=T)
  pred_linear_pp <- pred_linear_pp + 30.7873138968084
  
  ## Identify any individuals predicted as under 20s, and re-run with model trained on log(age)
  over20s <- names(pred_linear_pp[pred_linear_pp > 20])
  pred_linear_pp1 <- pred_linear_pp[over20s]
  under20s <- names(pred_linear_pp[pred_linear_pp < 20])
  
  ## Now re-run model for those
  coefficients_log <- coefficients_log[rownames(scores_log),]
  pred_log <- scores_log * coefficients_log[,"Coefficient"]
  pred_log_pp <- colSums(pred_log,na.rm=T)
  pred_log_pp <- pred_log_pp + 3.87249881475691
  pred_log_pp <- exp(pred_log_pp[under20s])
  predictions <- c(pred_log_pp, pred_linear_pp1)
  results <- data.frame("ID" = names(predictions), "Epigenetic Age (Bernabeu)" = predictions,check.names=F)

  out=merge(out,results,by="ID")
  out=out[,c(1,2,ncol(out),3:(ncol(out)-1))]
## combine Sex and Age information
if(!exists("sexageinfo")){
  out$'True Age' <- NA 
  out$Sex <- NA 
} else { 
  ids = out$ID
  sexageinfo = sexageinfo[match(ids, sexageinfo$ID),] 
  out$'True Age' <- sexageinfo$Age
  out$Sex <- sexageinfo$Sex
  message("4. Sex and Age Info Added")
}

## Correct order of output file 

out <- out[,c(1,ncol(out),(ncol(out)-1), 2:c(ncol(out)-2))]



## combine covariates  information
if(!exists("covariates")){
NULL 
  } else { 
  ids = out$ID
  covariates = covariates[match(ids, covariates$ID),] 
  out <- merge(out, covariates, by = "ID")
  for(i in names(out)[5:10]){ 
    out[,i] <- resid(lm(as.formula(paste(paste(bt(i),"~"),paste(bt(paste(c(names(covariates)[-which(names(covariates) %in% "ID")]))), collapse = "+"))), na.action = na.exclude, data = out))
  }
  out <- out[,c(1:10)]
  message("5. Covariates")
}

write.csv(out, "episcores_Les_120.csv", row.names = F)

## Save File and Finish Up 
message("Analysis Finished! Thank you for using our application. Output File is called \"out\"") 
