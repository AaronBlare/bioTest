
## Installation of Packages

my_packages <- c("shiny", "shinyWidgets", "shinythemes", "data.table", "shinyalert") 
## Find the non-installed ones 
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]
## install them
if(length(not_installed)) install.packages(not_installed)                             

## Load packages 
library(shiny)
library(shinyWidgets)
library(shinythemes)
library(data.table)
library(shinyalert)


cpgs <- read.csv("./Predictors_Shiny_by_Groups.csv", header = T) 
bt <- function(x) sprintf("`%s`", x)

## Read in separate files for Bernabeu age prediction as it requires two-step process 
coefficients <- read.delim("./cage_coefficients_linear.tsv", sep = "\t")
coefficients_log <- read.delim("./cage_coefficients_log.tsv", sep = "\t")
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
coefs <- coefficients[rownames(coefficients)[which(!(rownames(coefficients) %in% rownames(coef_2)))],,drop=FALSE]

# Beta means - for Bernabeu predictor (processed separately from other models)
means <- read.delim("cpg_meanbeta_gs20k.tsv")
row.names(means) <- means[,"cpg"]


timeoutSeconds <- 300

inactivity <- sprintf("function idleTimer() {
                      var t = setTimeout(logout, %s);
                      window.onmousemove = resetTimer; // catches mouse movements
                      window.onmousedown = resetTimer; // catches mouse movements
                      window.onclick = resetTimer;     // catches mouse clicks
                      window.onscroll = resetTimer;    // catches scrolling
                      window.onkeypress = resetTimer;  //catches keyboard actions
                      
                      function logout() {
                      Shiny.setInputValue('timeOut', '%ss')
                      }
                      
                      function resetTimer() {
                      clearTimeout(t);
                      t = setTimeout(logout, %s);  // time is in milliseconds (1000 is 1 second)
                      }
                      }
                      idleTimer();", timeoutSeconds*1000, timeoutSeconds, timeoutSeconds*1000)
## User Interface Portion 
ui <- fluidPage(theme = shinytheme("united"),
                tags$head(
                  tags$style(      ".title {margin-left: 370px; width: 1000px}",
                    HTML(".shiny-notification {
                         position:fixed;
                         top: calc(50%);
                         left: calc(50%);
                         }
                         "
                         
                    ),
                    ".modal-body {padding: 20px}
                    .modal-header {background-color: #8ae9ff}
                    .modal { text-align: center}"
                  )
                    ),    tags$div(class="title", titlePanel("MethylDetectR - Calculate Your Scores")), 
                tags$script(inactivity),
                sidebarLayout(
                  sidebarPanel(
                    useShinyalert(),
                    fileInput("file1", "Upload RDS File - Methylation",
                              multiple = TRUE,
                              accept = c(".rds")),
                    
                    fileInput("file2", "Upload CSV File - Sex and Age Information (Optional)",
                              multiple = TRUE,
                              accept = c("text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv")),
                    
              
                    
                    actionButton("go", "Run the Analysis"),
                    
                    br(),
                    br(),
                    
                    uiOutput("downloadData"),
                    
                    
                    br(),
                    actionButton("reset", "Press Here to Reset - Please Reset If Running New Dataset"),
                    br(),
                    br(),
                    actionButton("version", "Please Cite The Following Papers - Press Here"),
                    br(), 
                    br(),
                    actionButton("press", "Press Here For General Information and Useful Links"),
                    br(),
                    br(),
                    actionButton("format", "Press Here For File Formats and Useful Links")),
                  mainPanel(tabsetPanel(tabPanel("Calculating the Scores", verbatimTextOutput("text"), verbatimTextOutput("text1"),
                                                 verbatimTextOutput("text2"), verbatimTextOutput("text_m2beta"), verbatimTextOutput("text3"), 
                                                 verbatimTextOutput("text4"), verbatimTextOutput("text4b"),verbatimTextOutput("text7"), verbatimTextOutput("text8"),
                                                 verbatimTextOutput("text9")
                  )))
                ))

## Server Portion 
server <- function(input, output,session) {
  options(shiny.maxRequestSize = 5000*1024^2)
  
  
  dat <- reactive ({ req(input$file1)
    df1 <- readRDS(input$file1$datapath)
    output$text <- renderText("Finished Reading Methylation File")
    return(df1) 
  })
  
  
  
  sexageinfo <-  reactive({ 
    if (is.null(input$file2)) return(NULL)
    df2 <- read.csv(input$file2$datapath)
    return(df2)
  })
  covariates <-  reactive({ 
    if (is.null(input$file3)) return(NULL)
    df3 <- read.csv(input$file3$datapath)
    return(df3)
  })
  
  
  observeEvent(input$reset, { 
    session$reload()    
  } ) 
  
  dat3 <- reactive({ if(ncol(dat()) <= nrow(dat())){
    data1 <- dat()
    return(data1)} else{ 
      data2 <- t(dat())
      return(data2)}  }) 
  
  
  
  
  url <- a("'MethylDetectR'", href="https://shiny.igmm.ed.ac.uk/MethylDetectR/", target = "_blank")
  url_demo <- a("'MethylDetectR' Demo", href="https://shiny.igmm.ed.ac.uk/MethylDetectR_Demo/", target = "_blank")
  url2 <- a("'MethylDetectR' Website", href="https://www.ed.ac.uk/centre-genomic-medicine/research-groups/marioni-group/methyldetectr", target = "_blank")
  url3 <- a("Example Input and Output Files in Zenodo", href="https://zenodo.org/record/4291880#.YFzuUq_7TIU", target = "_blank")
  url4 <- a("'MethylDetectR Paper'", href="https://wellcomeopenresearch.org/articles/5-283", target = "_blank")
 
  p1 <- a("McCartney & Hillary - Lifestyle/Biochemical", href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1514-1", target = "_blank")
  p2 <- a("Gadd - Proteins", href="https://elifesciences.org/articles/71802", target = "_blank")
  p3 <- a("Zhang - Epigenetic Age", href="https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0667-1", target = "_blank")
  p4 <- a("Bernabeu - Epigenetic Age", href="https://www.biorxiv.org/content/10.1101/2022.09.08.507115v1", target = "_blank")
  p5 <- a("Hillary & Marioni - MethylDetectR", href="https://wellcomeopenresearch.org/articles/5-283", target = "_blank")
  
  
  observeEvent(input$version,{ 
    # show pop-up ...
    showModal(modalDialog(
      title = "Please Cite These Papers If You Use This App",
      div(HTML(paste0("1. Lifestyle and Biochemical Trait Predictors: McCartney and Hillary ",  em("et al."), " (2018). PMID: 30257690")),style="font-size:90%"),
      tagList(p1),
      br(),
      br(),
      div(HTML(paste0("2. Protein Trait Predictors: Gadd ",  em("et al."), " (2022). PMID: 35023833")),style="font-size:90%"),
      tagList(p2),
      br(),
      br(),
      div(HTML(paste0("3. Epigenetic Age: Zhang ",  em("et al."), " (2019). PMID: 30257690")),style="font-size:90%"),
      tagList(p3),
      br(),
      br(),
      div(HTML(paste0("4. Epigenetic Age: Bernabeu ",  em("et al."), " (2022). PMID: pending")),style="font-size:90%"),
      tagList(p4),
      br(),
      br(),
      div(HTML(paste0("5. MethylDetectR: Hillary and Marioni (2021). PMID: 33969230")),style="font-size:90%"),
      tagList(p5),
      br(),
      br(),
      easyClose = TRUE,
      footer = NULL))
    
  })
  
  
  
  observeEvent(input$press,{ 
    # show pop-up ...
    showModal(modalDialog(
      title = "Thank you for using 'MethylDetectR - Calculate Your Scores'",
      div("DNA methylation is an important biological process that helps to determine whether genes are switched on or off. It occurs at CpG sites across the genome and is influenced by our genetics and environment.", style="font-size:90%"),
      br(),
      br(),
      div("To develop the predictors in this app, we can for instance look at patterns of DNA methylation in smokers and non-smokers. We can identify general differences between these groups and make algorithms using these differences to estimate smoking status of other individuals.", style="font-size:90%"),
      br(),
      br(),
      div("These predictors work well at the level of the population by looking at general differences. However, they can make inaccurate predictions at an individual level. They will improve with more refined algorithms and bigger studies.", style="font-size:90%"),
      br(),
      br(),
      div("Applications will time out after three minutes of inactivity. Important information on data privacy and protection can be accessed in our website using the link below.", style="font-size:90%"),
      br(),
      br(),
      div("Please use the following links to access the other parts of our platform.", style="font-size:90%"),
      br(),
      tagList(url),
      br(),
      tagList(url_demo),
      br(),
      tagList(url2),
      br(),
      tagList(url3),
      br(),
      tagList(url4),
      
      easyClose = TRUE,
      footer = NULL))
    
  })
  
  observeEvent(input$format,{ 
    # show pop-up ...
    showModal(modalDialog(
      title = "How to Format Your Files",
      div("DNA methylation data should be uploaded as an .rds file. Uploads of greater than 500 Mb are discouraged. To make your file smaller, please truncate CpG sites to those in the 'Truncate_to_these_CpGs.csv' file in Zenodo (link at bottom).", style="font-size:90%"),
      br(),
      div("An R script to generate methylation-based scores without the app is also made available in the Zenodo repository.", style="font-size:90%"),
      br(),
      div("We recommend to have individuals as columns and CpG sites as rows but either version is accepted. Beta or M values are accepted, as are missing values.", style="font-size:90%"),
      br(),
      br(),
      div("An optional 'SexAgeInfo' file can be uploaded. This is a .csv file and should have three columns: ID's of individuals in the methylation file ('ID' column), a 'Sex' column for these individuals written as 'Male' or 'Female' or 'NA' and an 'Age' column detailing chronological ages of the individuals.", style="font-size:90%"),
      br(),
      br(),
      div("Once you have downloaded the DNAm-based scores, you may wish to adjust them for covariates, such as cell type counts or proportions. Importantly, when the predictors were created, they were trained on phenotypic data that were not adjusted for cell-type heterogeneity. The need for covariate adjustments will be specific to the aims of each study.", style="font-size:90%"),
      br(),
      br(),
      div("Applications will time out after three minutes of inactivity. Important information on data privacy and protection can be accessed in our website using the link below.", style="font-size:90%"),
      br(),
      div("Please use the following links to access the other parts of our platform.", style="font-size:90%"),
      br(),
      tagList(url),
      br(),
      tagList(url_demo),
      br(),
      tagList(url2),
      br(),
      tagList(url3),
      br(),
      tagList(url4),
      easyClose = TRUE,
      footer = NULL))
    
  })
  
  
  observeEvent(input$go, { 
    
    if(is.null(input$file1))
    {
      # show pop-up ...
      showModal(modalDialog(
        title = "Methylation File Not Yet Uploaded",
        paste0("Please Allow For Upload. To Dismiss Click Anywhere Else On Page"),
        easyClose = TRUE,
        footer = NULL
      ))} 
    else if(nrow(dat()) == 0 | ncol(dat()) ==0){
      shinyalert("Sorry, An Error Has Occurred","There Were No Individuals Detected in the Dataset",type = "warning",  showCancelButton = TRUE, inputId = "foo4")
      observeEvent(input$foo4, {
        session$reload()
      })
    }
    
    else if(withProgress(message = "Performing Quality Control Checks", {
      t<- vector()
      for(i in 1:ncol(dat3())){
        t[i] <- is.numeric(dat3()[,i])
        incProgress(1/ncol(dat3()), detail = paste((i/ncol(dat3()))*100, "% complete"))
        Sys.sleep(0.2)
      } 
      length(which(as.data.frame(t)[,1] %in% "FALSE")) >  0}) ){
      shinyalert("Sorry, An Error Has Occurred","Non-Numeric Values Have Been Detected - Please Check You Have The Correct Input Format",type = "warning",  showCancelButton = TRUE, inputId = "foo1")
      observeEvent(input$foo1, {
        session$reload()
      })
    }  else { 
    
      showModal(modalDialog("Calculating Your Scores... Please wait...",footer =NULL))
      ## Set up function to convert NAs to sample mean 
      
      na_to_mean <-function(x){
        x[is.na(x)]<-mean(x,na.rm=T)
        return(x)
      }
      
      ## Set up function to convert M values to beta values 
      
      m_to_beta <- function (val) 
      { 
        beta <- 2^val/(2^val + 1)
        return(beta)
      }
      
      
      
      ## Transposing step 
      
      dat1 <- reactive({ 
        
        if(ncol(dat()) > nrow(dat())){
          output$text2 <- renderText("Data Is Being Transposed")
          data1 <-t(dat())
          return(data1)} else{ 
            output$text2 <- renderText("Data Does Not Need To Be Transposed") 
            data2 <- dat()
            return(data2)}
      }) 
      
      dat.cage <- reactive({ 
        if(ncol(dat()) > nrow(dat())){
          data1.1 <-t(dat())
          return(data1.1)} else{ 
            data2.1 <- dat()
            return(data2.1)}
      }) 
      
    
      # Overlap of CpGs with those for predictors + conversion of M to beta values if needs be
      ## All predictors less Bernabeu 
      coef <- reactive({ 
        coef1.1 <- dat1()[intersect(rownames(dat1()), cpgs$CpG_Site),]
        output$text3 <- renderText("Subsetting CpG Sites To Those Needed by Predictors...") 
        if((range(coef1.1,na.rm=T)> 1)[[2]] == "TRUE"){
        output$text_m2beta <- renderText("M Values Detected - Converting to Beta Values") 
        fix <- m_to_beta(coef1.1)
        return(fix)
        } else { 
          output$text_m2beta <- renderText("Beta Values Detected") 
          safe <- coef1.1 
          return(safe)
          }
        })
      
      
      ## Separate dataframe for Bernabeu predictor 
      coef.cage <- reactive({ 
        coefcage1.1 <- dat.cage()[intersect(rownames(dat.cage()),all_cpgs),]
        if((range(coefcage1.1,na.rm=T)> 1)[[2]] == "TRUE"){
          fix1 <- m_to_beta(coefcage1.1)
          return(fix1)
        } else { 
          safe1 <- coefcage1.1
          return(safe1)
        }
      })
      
      
      
      ## All predictors less Bernabeu 
      coef1 <- reactive({ if(nrow(coef()) == length(unique(cpgs$CpG_Site))) {  
        output$text4 <- renderText("No Missing CpG Sites") 
        x <- coef()
        return(x) 
      } else if(nrow(coef()) == 0) { 
        shinyalert("Sorry, There Are No Neccessary CpG Sites Present", "This Means All Individuals Will Have Same Values For Predictors. Analysis Will Not Be Informative",type = "warning",  showCancelButton = TRUE, inputId = "foo2")
        observeEvent(input$foo2, {
          session$reload()
        })
      } 
        
        else {
          missing_cpgs = cpgs[-which(cpgs$CpG_Site %in% rownames(coef())),c("CpG_Site","Mean_Beta_Value")]
          message(paste(length(unique(missing_cpgs$CpG_Site)), "unique sites are missing - add to dataset with mean Beta Value from Training Sample", sep = " "))
          mat = matrix(nrow=length(unique(missing_cpgs$CpG_Site)),ncol = ncol(coef()))
          row.names(mat) <- unique(missing_cpgs$CpG_Site)
          colnames(mat) <- colnames(coef()) 
          mat[is.na(mat)] <- 1
          missing_cpgs1 <- if(length(which(duplicated(missing_cpgs$CpG_Site))) > 1) { 
            missing_cpgs[-which(duplicated(missing_cpgs$CpG_Site)),]
          } else {missing_cpgs
          }  
          ids = unique(row.names(mat))
          missing_cpgs1 = missing_cpgs1[match(ids,missing_cpgs1$CpG_Site),]
          mat=mat*missing_cpgs1$Mean_Beta_Value
          coef.2=rbind(coef(),mat)
          output$text4 <- renderText(paste(paste("Substituting in", length(unique(missing_cpgs$CpG_Site)), sep = " "), "Missing CpG Sites With Mean Beta Values from Training Sample"))
          return(coef.2)} 
      })
      
      ## Bernabeu predictor
      coef.cage1 <- reactive({ if(nrow(coef.cage()) == length(unique(means$cpg))) {  
        xcage <- coef.cage()
        output$text4b <- renderText("For Bernabeu Age Predictor - No Missing CpG Sites") 
        return(xcage) 
      } else if(nrow(coef.cage()) == 0) { 
        print(message("There Are No Neccessary CpG Sites for Bernabeu age predictor"))
      } 
        else {
          missing_cpgs.cage = means[-which(means$cpg %in% rownames(coef.cage())),c("cpg","mean")]
          print(message(paste(length(unique(missing_cpgs.cage$cpg)), "unique sites are missing for Bernabeu predictor - adding mean betas from Training Dataset", sep = " ")))
          mat.cage = matrix(nrow=length(unique(missing_cpgs.cage$cpg)),ncol = ncol(coef.cage()))
          row.names(mat.cage) <- unique(missing_cpgs.cage$cpg)
          colnames(mat.cage) <- colnames(coef.cage()) 
          mat.cage[is.na(mat.cage)] <- 1
          missing_cpgscage1 <- if(length(which(duplicated(missing_cpgs.cage$cpg))) > 1) { 
            missing_cpgs.cage[-which(duplicated(missing_cpgs.cage$cpg)),]
          } else {missing_cpgs.cage
          }  
          ids.cage = unique(row.names(mat.cage))
          missing_cpgscage1 = missing_cpgscage1[match(ids.cage,missing_cpgscage1$cpg),]
          mat.cage=mat.cage*missing_cpgscage1$mean
          coef.cage2=rbind(coef.cage(),mat.cage)
          output$text4b <- renderText(paste(paste("For Bernabeu Age Predictor - Substituting in", length(unique(missing_cpgscage1$cpg)), sep = " "), "Missing CpG Sites With Mean Beta Values from Generation Scotland"))
          return(coef.cage2)} 
      })
      
  
    
      
      coef3 <- reactive({ coef2.1 <- t(apply(coef1(),1,function(x) na_to_mean(x)))
      return(coef2.1)
      }) 
      
      coef.cage3 <- reactive({ coef2.2 <- t(apply(coef.cage1(),1,function(x) na_to_mean(x)))
      
      return(coef2.2)
      }) 

      out1 <- reactive({
        loop = unique(cpgs$Predictor)
        out <- data.frame()
        for(i in loop){ 
          tmp=coef3()[intersect(row.names(coef3()),cpgs[cpgs$Predictor %in% i,"CpG_Site"]),]
          tmp_coef = cpgs[cpgs$Predictor %in% i, ]
          if(nrow(tmp_coef)>1){ 
            tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG_Site),]
            out[colnames(coef3()),i]=colSums(tmp_coef$Coefficient*tmp)
            output$text7 <- renderText("Calculating the Predictors...") } else { 
              tmp2 = as.matrix(tmp)*tmp_coef$Coefficient 
              out[colnames(coef3()),i] = tmp2[,1]
              output$text7 <- renderText("Calculating the Predictors...")
            }
        }
        
        out$'Epigenetic Age (Zhang)' <- out$'Epigenetic Age (Zhang)' + 65.79295
        out$ID <- row.names(out) 
        out <- out[,c(ncol(out),1:(ncol(out)-1))] 
        
        ## Calculating age predictor - Bernabeu 
      
        scores <- coef.cage3()[which(row.names(coef.cage3()) %in% rownames(coefs)),]
        scores_quadratic <- coef.cage3()[which(row.names(coef.cage3()) %in% coef_2_simp),]**2
        rownames(scores_quadratic) <- paste0(rownames(scores_quadratic), "_2")
        scores_linear <- rbind(scores, scores_quadratic)
        
        ## Prep for log predictor
        scores_log <- coef.cage3()[which(row.names(coef.cage3()) %in% rownames(coef_log)),]
        scores_quadratic_log <- coef.cage3()[which(row.names(coef.cage3()) %in% coef_log_2_simp),]**2
        rownames(scores_quadratic_log) <- paste0(rownames(scores_quadratic_log), "_2")
        scores_log <- rbind(scores_log, scores_quadratic_log)
        
        ## Calculate cAge with linear model
        coefficients <- coefficients[which(row.names(coefficients) %in% rownames(scores_linear)),]
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
        
        if(is.null(sexageinfo())){
          out$'True Age' <- NA
          out$Sex <- NA 
          out3 <- out[,c(1,ncol(out),(ncol(out)-1), 2:c(ncol(out)-2))]
          return(out3)
        } else { 
          ids = out$ID
          linker = sexageinfo()[match(ids, sexageinfo()$ID),] 
          out$'True Age' <- linker$Age
          out$Sex <- linker$Sex 
          out3 <- out[,c(1,ncol(out),(ncol(out)-1), 2:c(ncol(out)-2))]
          return(out3)} 
      })
        
    
      output$text8 <- renderText({ 
        
        
        if(is.null(input$file2)){
        
        out.true.age.sex <- "No Sex/Age Info File - Sex and 'True Age' Listed As 'NA'"
        } 
        
        if(!is.null(input$file2)){ 
          
          out.true.age.sex <- "Sex/Age Info File Included - Sex and 'True Age' Incorporated Into Output File"
  
          }
        
        out.true.age.sex
        
      }) 
      
      
       
        
      output$text9 <-if(nrow(out1()) >0){ 
        renderText("Analysis is Finished. Thank You For Using Our Application") } 
      else { NULL 
       }

      
     
      output$text1 <- renderText({
        
        if(!is.null(input$file2)) {
          out.cov.sexage <-  "Finished Reading Sex and Age Info File"
        } else { 
          out.cov.sexage <- NULL 
          }
        
        out.cov.sexage
      })
      
      
      
     
      ## Save File and Finish Up 
   
      
      
      output$downloadData <- renderUI({
          downloadButton("download", "Download Your Results")

        
      })
      
      
     
       
        
     

      output$download <-  downloadHandler(
        filename = function(){paste("Epigenetic_Scores", ".csv", sep = "")}, 
        content = function(file){
          fwrite(out1(), file, row.names = F,dateTimeAs = "write.csv")
          
        } ) 
        
        
   
      removeModal()
    }
  })
  
  observeEvent(input$timeOut, { 
    print(paste0("Session (", session$token, ") timed out at: ", Sys.time()))
    showModal(modalDialog(
      title = "Timeout",
      paste("Session Time Due To", input$timeOut, "Inactivity -", Sys.time()),
      footer = NULL
    ))
    session$close()
    stopApp()
  })
  
}

shinyApp(ui=ui,server=server)
