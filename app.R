# COVID-19 HLA DATABASE PORTAL - 
# v0.6.2
# By: Livia Tran 
# October 20, 2020
# Date this version deployed: October 20, 2020


library(shiny)
library(HIBAG)
library(shinyjs)
library(httr)
library(jsonlite)
library(rmarkdown)
library(rjson)
library(shinybusy)
library(dplyr)

options(warn=2)
options(shiny.maxRequestSize = 100*1024^2)

#create log file if not already in existence 
if(!file.exists("/srv/shiny-server/HLA-Imputation-Portal/log_file.log")){
    file.create("/srv/shiny-server/HLA-Imputation-Portal/log_file.log")
}

#writeLog
#singularizes parameters for log message appending 
writeLog<-function(mess){
    write.table(paste(Sys.time(), sep=" ", mess), file="/srv/shiny-server/HLA-Imputation-Portal/log_file.log", append = T, row.names=F, col.names = F)
}

HIBAGgerVANCE <- function(ethnicity, genoMethod, SNPdata, useremail, project){
    
    
    h <- c("hibagger@hlacovid19.org", "712ZTRhZBKPzZmmGxD56", "application/json")
    
    names(h) <- c("X-User-Email", "X-User-Token", "Content-type")
    
    content <-fromJSON(content(POST("https://database-hlacovid19.org/query/hibag_preflight",
                                    body = paste('{"email":"',useremail,'","project_name":"',project,'","origin_identifiers":', toJSON(SNPdata$sample.id),'}', sep = ""),
                                    add_headers(.headers = h), encode='json'), "text"))
    
    
    if(content$user_approved == TRUE & content$project_found == TRUE){
        projValidity<-TRUE
        
        #useremail="livia.tran@ucsf.edu"
        #project= "sm_test"
        writeLog(paste0("The e-mail ", useremail, " and ", project, " project are valid."))
        
    }
    
    if(content$user_approved == TRUE & content$project_found == FALSE){
        projValidity<-paste0("The ",  project, " project has not been registered in the COVID-19|HLA Database. A registered project name is required for imputation.")
        
        writeLog(paste0("The e-mail ", useremail, " is valid, but the ", project, " project is not valid."))
    }
    
    
    if(content$user_approved == FALSE){
        projValidity<-paste0("The ",  useremail, " email address has not been registered in the COVID-19|HLA Database. A registered e-mail is required for imputation.")
        
        writeLog(paste0("The e-mail ", useremail, " is not valid."))
    }
    
    
    if(projValidity!=TRUE){
        return(projValidity)
    }
    
    #turn origin_identifiers response list into a dataframe 
    subjectDF<-data.frame(validity=unlist(content$origin_identifiers))
    subjectDF$OIDs<-rownames(subjectDF)
    rownames(subjectDF)<-NULL
    
    if(all((subjectDF$validity %in% "valid"))==FALSE){
        
        writeLog("No origin identifiers present in the database.")
        return("None of the origin identifiers are present in the database.")}
    
    
    if(any((subjectDF$validity %in% "valid")==FALSE & all((subjectDF$validity %in% "valid"))==FALSE)){
        
        writeLog(paste(subjectDF$OIDs[which(!subjectDF$validity %in% "valid")],collapse=", ")," origin identifiers are not present in the database.")
        
        return(paste0("The following origin identifiers [ " ,paste(subjectDF$OIDs[which(!subjectDF$validity %in% "valid")],collapse=", "), sep = "] ", "are not found in the database."))
    }
    
    if(all((subjectDF$validity %in% "valid"))==TRUE){
        subjectValidity<-TRUE
        
        writeLog("All origin identifers are present in the database.")
        
        
        if(projValidity == TRUE & subjectValidity == TRUE){
            
            hla.id<-c("A", "B", "C", "DRB1", "DQA1", "DQB1", "DPB1")
            
            predictions<-imputations<-stats<-modelConvert<-sapply(hla.id, function(x) NULL)
            
            #convert saved hlaObjects to hla Models
            if(file.exists(paste("Models/", paste(paste(genoMethod, ethnicity, sep="_"), ".rda", sep=""), sep=""))){
                load(paste("Models/", paste(paste(genoMethod, ethnicity, sep="_"), ".rda", sep=""), sep=""))    
                
                writeLog(paste0(genoMethod, sep="_", ethnicity, " model loaded."))
            }
            
            else{
                return(paste0("An imputation model for the ", ethnicity, " group is not available for the ", genoMethod, " platform."))
                writeLog(paste0("An imputation model for the ",ethnicity, " group is not available for the ", genoMethod, " platform."))
            }
            
            
            for(i in 1:length(hla.id)){
                
                #convert saved hlaObjects to hla Models
                modelConvert[[i]]<-hlaModelFromObj(get(paste(genoMethod, ethnicity, sep="_"))[[hla.id[[i]]]])
                
                
                #writeLog("HLA model converted from object.")
                
                #if ethnicity for model entered does not exist, return error message 
                if(length(modelConvert[[i]])==1){
                    
                    writeLog(paste0("An imputation model for the ", ethnicity, " group is not available for the ",  genoMethod, " platform."))
                    
                    return(paste(paste(paste("An imputation model for the", ethnicity), "group is not available for the", genoMethod), "platform."))
                }
                
                
                predictions[[i]]<-tryCatch({predict(modelConvert[[hla.id[[i]]]], SNPdata)}, error =function(error) {geterrmessage()})
                
                
                #sub out conversion from warning to error message 
                #sub out conversion from warning to error message 
                if(length(predictions[[i]])==1){
                    if(grepl("warning", predictions[[i]])){
                        
                        writeLog(paste0("ERROR: ", predictions[[i]]))
                        
                        #add error/warning message after this line
                        return(trimws(gsub("[(]converted from warning[)]", "", predictions[[i]])))
                    }
                }
                
                
                #if length is one, error message is present 
                if(length(predictions[[i]]) == 1){
                    
                    writeLog(paste0("ERROR: ", predictions[[i]]))
                    
                    return(predictions[[i]])
                }
                
                writeLog(paste0("Imputation completed for HLA-", hla.id[[i]], "."))
                
                
                names(predictions[[i]]$value)[names(predictions[[i]]$value) == 'allele1'] <- paste(hla.id[[i]], "_1", sep="")
                names(predictions[[i]]$value)[names(predictions[[i]]$value) == 'allele2'] <- paste(hla.id[[i]], "_2", sep="")
                names(predictions[[i]]$value)[names(predictions[[i]]$value) == 'prob'] <- paste(hla.id[[i]], "prob", sep="_")
                names(predictions[[i]]$value)[names(predictions[[i]]$value) == 'matching'] <- paste(hla.id[[i]], "matching", sep="_")
                
                
                imputations[[i]] <- predictions[[i]]$value[,c(1:3)]
                imputations[[i]][,2]<-paste(hla.id[[i]], imputations[[i]][,2], sep="*")
                imputations[[i]][,3]<-paste(hla.id[[i]], imputations[[i]][,3], sep="*")
                
                
                stats[[i]] <-predictions[[i]]$value[,c(1,4:5)]
                
            }
            
            imputationsDF<-Reduce(function(...) merge(..., by="sample.id", all=TRUE), imputations)
            statsDF <- Reduce(function(...) merge(..., by="sample.id", all=TRUE), stats)
            
            imputationsDF<-dplyr::rename(imputationsDF, origin_identifier = "sample.id")
            statsDF<-dplyr::rename(statsDF, origin_identifier = "sample.id")
            
            imputationsDF<-tibble::add_column(imputationsDF, typing_method_name = "HIBAG")
            statsDF<-tibble::add_column(statsDF, typing_method_name = "HIBAG")
            
            imputationsDF<-tibble::add_column(imputationsDF, typing_method_version = "3.1.1.")
            statsDF<-tibble::add_column(statsDF, typing_method_version = "3.1.1.")
            
            imputationsDF<-tibble::add_column(imputationsDF, pop = ethnicity)
            statsDF<-tibble::add_column(statsDF, pop = ethnicity)
            
            write.csv(imputationsDF, file=paste("/var/hibag_output/", paste(paste(paste(paste("__", gsub(".", "_dot_", gsub("@", "_at_", useremail), fixed=T), sep=""), "__", sep=""), paste(project, "__", sep=""), sep=""), paste(ethnicity,"__imputations__", Sys.Date(), ".csv", sep=""), sep=""), sep=""), row.names = F)
            
            writeLog(paste("/var/hibag_output/", paste(paste(paste(paste("__", gsub(".", "_dot_", gsub("@", "_at_", useremail), fixed=T), sep=""), "__", sep=""), paste(project, "__", sep=""), sep=""), paste(ethnicity,"__imputations__", Sys.Date(), ".csv", sep=""), sep=""), sep=""))
            
            write.csv(statsDF, file=paste("/var/hibag_output/", paste(paste(paste(paste("__", gsub(".", "_dot_", gsub("@", "_at_", useremail), fixed=T), sep=""), "__", sep=""), paste(project, "__", sep=""), sep=""), paste(ethnicity,"__stats__", Sys.Date(), ".csv", sep=""), sep=""), sep=""), row.names = F)
            
            writeLog(paste("/var/hibag_output/", paste(paste(paste(paste("__", gsub(".", "_dot_", gsub("@", "_at_", useremail), fixed=T), sep=""), "__", sep=""), paste(project, "__", sep=""), sep=""), paste(ethnicity,"__stats__", Sys.Date(), ".csv", sep=""), sep=""), sep=""))
            
            write.table(paste0(Sys.time(), sep=" ", "Imputations for all loci are complete."), file="/srv/shiny-server/HLA-Imputation-Portal/log_file", append=T, row.names=F, col.names=F)
            
            return("Imputation completed. HLA genotypes and statistics have been loaded to the database.")
            
            
        }
    }
}



ui = fluidPage(
    add_busy_spinner(spin = "fading-circle"),
    useShinyjs(),
    titlePanel("COVID-19|HLA & Immunogenetics Consortium HLA Imputation Portal"),
    sidebarPanel(tags$head(
        tags$style(type="text/css", "select { max-width: 240px; }"),
        tags$style(type="text/css", ".span4 { max-width: 290px; }"),
        tags$style(type="text/css", ".well { max-width: 280px; }")
    ),
    textInput("email", h4("Submitter's e-mail")),
    textInput("projID", h4("Project Name")),
    selectInput("ethnicity", h4("Ethnic population"), 
                choices = list("Subsaharan African/African American" = "African", "East Asian/Asian American" = "Asian", "European/European American" = "European", "Hispanic/Hispanic American" = "Hispanic", "Multi-ethnic/Other" = "Broad"), multiple=F),
    selectInput("gmethod", h4("Genotyping Method"),
                choices = list("Affymetrix Genome-Wide Human SNP Array 5.0" = "Affy500K", "Affymetrix Axiom® Genome-Wide ASI 1 Array Plate" = "AffyAxiom_ASI", "Affymetrix Axiom® Human Origins 1 Array" = "AffyAxiom_Human_Origin_1", "Affymetrix Axiom® UK Biobank Array" = "AffyAxiomUKB", "Affymetrix Axiom® World Arrays Genome-Wide EAS 1 Array Plate" = "AffyAxiom_WA_EAS", "Affymetrix Axiom® Genome-Wide CHB 1 Array Plate-1" = "AffyAxiomCHB-1", "Affymetrix Axiom® Exome 319 Array Plate" = "AffyAxiomExome319", "Affymetrix Axiom® World Arrays Genome-Wide EUR 1 Array Plate" = "AffyAxiom_WA_EUR", "Affymetrix Genome-Wide Human SNP Array 6.0" = "AffySNPv6",
                               "Axiom™ Precision Medicine Research Array (PMRA)" = "Axiom_PMRA", "Illumina Human610-Quad v1.0" = "Human610", "Illumina Human660W-Quad v1.0" = "Human660W", "Illumina HumanCoreExome" = "HumanCoreExome", "Illumina HumanExome" = "HumanExome12", "Illumina HumanOmni2.5" = "HumanOmni2.5", "Illumina HumanOmni5" = "HumanOmni5", "Illumina Infinium HumanOmni5Exome" = "HumanOmni5Exome", "Illumina HumanOmniExpress " = "HumanOmniExpress", 
                               "Illumina Sentrix HumanHap300" = "IlluminaHumanHap300", "Illumina ImmunoChip" = "Immunochip", "Illumina Infinium Global Screening Array v2.0" = "InfiniumGlobal", "Illumina Infinium Multi-Ethnic Global BeadChip" = "InfiniumMEGA", "Illumina Infinium OmniExpress-24 v1.2" = "InfiniumOmniExpress24", "Illumina MHC Exon-centric Golden Gate panel" = "MHC_Exon_Golden_Gate",
                               "Illumina HumanOmni1-Quad" = "Omni1_Quad", "Illumina Infinium Omni2.5Exome-8 BeadChip Kit" = "Omni2.5Exome", "Illumina Human Omni Zhonghua-8 Beadchip" = "OmniZhonghua", "Illumina Infinium OncoArray-500K BeadChip" = "OncoArray500K", "Illumina PsychArray-B" = "PsychArray-B"), multiple = F),
    tags$head(tags$style('h5 {color:red;}')),
    
    
    h5("REQUIRED files for imputation"),
    fileInput("file1", ".bed",
              accept = ".bed", multiple=F),
    fileInput("file2", ".bim",
              accept = ".bim", multiple=F),
    fileInput("file3", ".fam",
              accept = ".fam", multiple=F),
    HTML("<br>"),
    
    actionButton("run", "Impute!"),
    HTML("<br>"),
    HTML("<br>"),
    actionButton("reset", "Clear"),
    
    ),
    mainPanel(
        tabsetPanel(type="tabs",
                    
                    tabPanel("Main", 
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             span(span(textOutput("endMess"), style="font-size:20px"), style="font-weight:bold"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),
                             
                             fluidRow(column(12,img(src=base64enc::dataURI(file="logo.png", mime="image/png"), align="right", width="300"), height = 200)),
                             
                    ),
                    
                    tabPanel("About", uiOutput("aboutText"),     
                             HTML("<br>"), HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"), HTML("<br>"), HTML("<br>"), HTML("<br>"),
                             fluidRow(column(12,img(src=base64enc::dataURI(file="logo.png", mime="image/png"), align="right", width="300"), height = 200)),)),
    ))


server = function(input, output, session) {
    
    bedFile<-reactiveValues(bed=NULL)
    bimFile<-reactiveValues(bim=NULL)
    famFile<-reactiveValues(fam=NULL)
    
    fileClear <- reactiveValues(
        data = NULL,
        clear = FALSE
    )
    
    disable("run")
    
    isolate(observeEvent(input$file1, {
        bedUpload<-input$file1
        bedFile$bed<-bedUpload$datapath
        
    }))
    
    isolate(observeEvent(input$file2, {
        bimUpload<-input$file2
        bimFile$bim<-bimUpload$datapath
    }))
    
    isolate(observeEvent(input$file3, {
        famUpload<-input$file3
        famFile$fam<-famUpload$datapath
    }))
    
    #only enable action button if all files are uploaded
    isolate(observe({
        if(length(bedFile$bed!=0) & length(bimFile$bim!=0) & length(famFile$fam!=0)){
            enable("run")
        }
    }))
    
    endMess <- reactiveValues(functionOut=NULL)
    ethnicity<-reactiveValues(ethnic=NULL)
    genotypeMethod<- reactiveValues(gMethod = NULL)
    email<-reactiveValues(userEmail=NULL)
    project<-reactiveValues(projectID=NULL)
    
    
    #ethnicity
    observeEvent(input$ethnicity, {
        ethnicity$ethnic<-input$ethnicity
    })
    
    #gmethod
    observeEvent(input$gmethod, {
        genotypeMethod$gmethod<-input$gmethod
    })
    
    #email
    observeEvent(input$email, {
        email$userEmail<-input$email
    })
    
    #project ID
    observeEvent(input$projID, {
        project$projectID<-input$projID
    })
    
    #counter <- reactiveValues(countervalue = 0) # Defining & initializing the reactiveValues object
    
    
    observeEvent(input$run, { 
        showModal(modalDialog("Your imputation is running...", footer=NULL))
        
        SNPdata<-hlaBED2Geno(bedFile$bed, famFile$fam, bimFile$bim)                  
        
        endMess$functionOut<-tryCatch({HIBAGgerVANCE(ethnicity$ethnic, genotypeMethod$gmethod, SNPdata, useremail= email$userEmail, project=  project$projectID)}, error = function(error) {"Unknown error. Imputation terminated."})
        
        if(grepl("overlapping", endMess$functionOut)==TRUE){
            output$endMess<-renderText({"There are no overlapping SNPs between the selected training dataset and your SNP data."})
        }
        
        else{
            output$endMess<-renderText(endMess$functionOut)
        }
        
        if(grepl("50%", endMess$functionOut)==TRUE){
            output$endMess<-renderText({paste("More than 50% of the SNPs in the", paste(genotypeMethod$gmethod, ethnicity$ethnic, sep ="_"), "model are not present in your SNP data.", sep =" ")})
        }
        removeModal()
        
    })
    
    
    #clear button logic
    isolate(observeEvent(input$reset, {
        updateSelectInput(session, "ethnicity", selected = c("African"))
        updateSelectInput(session, "gmethod", selected = c("Affy500K"))
        SNPdata<-NULL
        fileClear$data <- NULL
        fileClear$clear <- TRUE
        updateTextInput(session, "email", value = "")
        updateTextInput(session, "projID", value = "")
        output$endMess<-NULL
        reset('file1')
        reset('file2')
        reset('file3')
        disable("run")
    }))
    
    
    #target = '_blank' forces HIBAG models hyperlink to automatically be opened in a new window 
    output$aboutText<-renderUI({
        HTML("<b>Version 0.6.2</b>
    <br>
    <br>
    <font size='4'><b><u>Authors:</b></font></u><br>
  <br>
  Livia Tran - <i>livia.tran@ucsf.edu</i> 
  <br>
  Steven Mack - <i>steven.mack@ucsf.edu </i>
  <br>
  <br>
  This portal applies the <b>HIBAG</b> R package, described by <i>Zheng et al. (2014)</i>, to impute genotypes for the HLA-A, -B, -C, -DRB1, -DQA1, -DQB1, -DPA1, and -DPB1 loci. The portal supports HLA imputation from data generated using 30 SNP genotyping platforms for individuals of African, Asian, European, Hispanic, or Multi-ethnic origin, using 129  
  models, which can be found <a href='https://zhengxwen.github.io/HIBAG/platforms.html' target='_blank'> here.</a> 
  <br>
  <br>
  After 'Impute!' is clicked, a busy spinner on the top right corner will appear to indicate that the application is running.
  <br>
  <br>
  <b>REQUIRED INFORMATION </b>
  <br>
  - Submitter's e-mail address : Please provide the email address associated with the submitter's database account. 
  <br>
  - Project Name: Please provide the database Project Name associated with the origin identifiers (defined in the 'origin_identifier' field of the data dictionary) for the individuals being imputed.
 <br>
 <br>
 <b> REQUIRED FILES</b>
 <br>
  - .bed, .bim, and .fam files for the subjects to be imputed. Origin identifiers included in the .fam file must have already been submitted to the database and associated with the Project Name.
  <br>
  <br>
  <b>NOTE</b>
  <br>
  If there are no SNPs in common between your data and the selected training model, the following warning message will be displayed: <i><b>'There are no overlapping SNPs between the selected training dataset and your SNP data'. </b> </i>
  <br>
  If more than 50% of the SNPs in your data are not included in the selected training model, the following warning message will be displayed: <i><b>'More than 50% of the SNPs in the [selected training model] are not present in your SNP data'. </b> </i>
     <br>
      <br>
    The COVID-19 | HLA & Immunogenetics Consortium HLA Imputation Portal has been developed with the support of National Institute of Allergy and Infectious Disease Grant R01 AI128775-07S1.")})
}



shinyApp(ui = ui, server = server)

