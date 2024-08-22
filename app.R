#caMeL - a chemometrics applied machine learning interface
#Copyright (C) 2024  Agnieszka Konkolewska

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.


list.of.packages <- c("shiny", "glmnet", "caret", "signal", "mdatools", "chemometrics","EMSC","plotly", "shinyFiles", "tidyverse", "pls",
                      "shinycssloaders", "DT", "rlang", "reticulate", "hyperSpec")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(plotly)
library(shiny)
library(shinyFiles)
library(shinycssloaders)
library(DT)
library(tidyverse)
library(rlang)
library(reticulate)


# source("read.spc.R")
# source("split.string.R")

#Load model
models <- read.table("<folderpath>/caMeL/models/modelSettings.txt") #hit enter at the end of the file
if (file.exists("<folderpath>/caMeL/models/Options.RDS")){
  ops <- readRDS("<folderpath>/caMeL/models/Options.RDS")
  mainDir <- ops$path
  setwd(mainDir)
}else{
  mainDir <- getwd()
}

predictSpectra <- function(spectraPred, modelscript){
  #takes two objects available in the environment: spectraPred and modelscript
  source(paste0("<folderpath>/caMeL",modelscript))
  predicted <<- predicted
  predMD <<- predMD
}


options(readr.show_col_types = FALSE)

convertToCSV <- function(each){
  spectrum <- hyperSpec::read.spc(each)
  spectra <- spectrum$spc
  
  newSpec.spc <<- spectra
  
  newName  <- basename(substr(each, 1, nchar(each)-4))
  newName <<- paste0(newName, ".csv")
  
}

convertToCSV2 <- function(each){
  spectrum <- hyperSpec::read.spc(each)
  meanSP <- apply(spectrum$spc,2, mean)
  
  #Create new mean spectra
  newSpec <- empty(spectrum, nrow = 1, spc=meanSP) #get features from the original hyperSpec object
  newSpec@data$filename <- basename(each)
  
  #Convert to data.frame
  newSpec.spc <- newSpec@data
  newSpec.spc <<- newSpec$spc 
  #matplot(newSpec@wavelength, t(newSpec.spc), type="l", main = paste0(basename(each)), ylab = "Reflectance", xlab = "Wavelength [nm]", lty=1)
  newName <-  substr(each, 1, nchar(each)-4)
  newName <<- paste0(substr(newName,3, nchar(newName)), ".csv")
  
}

if (!file.exists("predictions.csv")){
  #old_files <- read.table("fileList.txt")
  print("No prediction file")
  n <- length(models$ModelPrefix)
  
  predictTrait <- data.frame(matrix(ncol=(2+n),nrow=0, dimnames=list(NULL, c("Sample", paste0("New",1:(n+1))))))
  colnames(predictTrait)[-1] <- c(models$ModelPrefix, "Timestamp")
  predictTrait <<- predictTrait
  write.csv(predictTrait, "predictions.csv", quote=FALSE, row.names = (FALSE))
}else{
  predictTrait <- as.data.frame(read_csv("predictions.csv"), col_types = cols(), show_col_types = FALSE)
  predictTrait$Timestamp <- as.character(predictTrait$Timestamp)
  predictTrait <<- predictTrait
}
if (!file.exists("fileList.txt")){
  
  old_files <- c("Start")
  write.table(old_files, "fileList.txt", row.names=FALSE, col.names=FALSE, quote=FALSE) #create empty list
}else{
  
  old_files <- read.table("fileList.txt")
  if(length(old_files)==0){
    old_files <- c("Start")
    write.table(old_files, "fileList.txt", row.names=FALSE, col.names=FALSE, quote=FALSE) #create empty list
  }
  
}

predictTrait <- as.data.frame(read_csv("predictions.csv"), col_types = cols(), show_col_types = FALSE)



# Define UI for application
ui <- navbarPage("caMeL",
                 tabPanel("Home",
                          
                          # Sidebar
                          sidebarLayout(
                            sidebarPanel(
                              
                              uiOutput("buttonPanel"),
                              uiOutput("traitPanel")
                              
                            ),
                            mainPanel(
                              tabsetPanel(
                                #textOutput("dir"),
                                tabPanel(title="Results", textOutput("chosen"),
                                         plotlyOutput("predsPlot"),
                                         textOutput("info")%>% withSpinner(color="#0dc5c1"),
                                         plotOutput("spectra")
                                ),
                                tabPanel(title="Dataset",textOutput("info2"),
                                         DT::dataTableOutput("mytable")
                                ),
                                tabPanel(title="Sample",
                                         DT::dataTableOutput("distTable")
                                )
                              )
                            )
                          )
                 ),
                 tabPanel("Settings",
                          sidebarPanel(
                            shinyDirButton('folder', 'Select a folder', 'Please select a folder with .spc files', FALSE),
                            uiOutput("optionsPanel")
                          )
                 )
)

# Define server logic
server <- function(input, output,session) {
  
  volumes <- getVolumes()

  observe({
    shinyDirChoose(input, 'folder', roots=volumes(), filetypes=c('', 'txt'))
    print(input$folder)
    
  })
  
  newFolder <- reactive({
    
    output$buttonPanel <- renderUI ({
      actionButton("button", "Start scanning")
    })
    
    validate({
      need(!is.atomic(input$folder), label = "Proper folder ")
    } )
    
    rootFol <- gsub("\\)", "", gsub("\\(", "", stringr::str_extract(input$folder$root, "\\((.*?)\\)")))
    text <- paste0(rootFol,"/",paste(sapply(input$folder$path[-1], unlist), collapse = "/"))
    text <- gsub("^NA", "", text)
    pathToFolder$path <<- paste0(text,"/")
    print(pathToFolder$path)
    setwd(pathToFolder$path)
    
    mainDir <<- pathToFolder$path
    text
    
  })
  
  pathToFolder <- reactiveValues(path=mainDir)
  
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$folder
               },
               handlerExpr = {
                 
                   output$optionsPanel <- renderUI ({
                     
                     
                     checkboxGroupInput("options", label = h3("Options"), choices = list("Convert (.spc to .csv)" = 1, "Predict from NIRS" = 2, "Calculate mean spectra from a multispec file" = 3, "Merge files" = 4),
                                        selected = c(1,2)) # ops$inputOptions
                     
                   })
                 

                 
                 output$chosen = renderText({ newFolder()})
                 
                 print(pathToFolder$path)
               })
  
  
  #text_to_display2 <-mainDir
  
  output$chosen <- renderText({  mainDir })
  
  autoInvalidate <- reactiveTimer(30000)
  
  observe({ #todo make sure the input$options is selected (e.g. do not run this, if it is empty)
    # Invalidate and re-execute this reactive expression every time the
    # timer fires.
    autoInvalidate()
    
    if(clicked()==1){

      savedPath <- pathToFolder$path
      options <- list(path= savedPath, inputOptions = input$options)
      saveRDS(options, "<folderpath>/caMeL/models/Options.RDS")
      setwd(pathToFolder$path)

      if (!file.exists("fileList.txt")){
        old_files <- c("Start")
        write.table(old_files, "fileList.txt", row.names=FALSE, col.names=FALSE, quote=FALSE) #create empty list
      }else{
        old_files <- read.table("fileList.txt")

      }
      
      subDir <<- "converted"
      
      if (file.exists(file.path(subDir))){
        print("A folder for converted spc files exists already at this folder.")
      }else {
        print("Creating new folder for converted data.")
        dir.create(file.path(subDir))
      }

      
      if (!file.exists("distanceMs.txt")){
        n <- length(models$ModelPrefix)
        parts <- c(" MD", " CN")
        library(stringi)

        distance <- data.frame(matrix(ncol=(1+2*n),nrow=0, dimnames=list(NULL, c("Sample", rep(models$ModelPrefix, each=2) %s+% rep(parts, times=n)))))
        write.table(distance, "distanceMs.txt", quote = FALSE, row.names=FALSE, sep = "\t")
        }
      if (!file.exists("predictions.csv")){
        #old_files <- read.table("fileList.txt")
        n <- length(models$ModelPrefix)
        
        predictTrait <- data.frame(matrix(ncol=(2+n),nrow=0, dimnames=list(NULL, c("Sample", paste0("New",1:(n+1))))))
        colnames(predictTrait)[-1] <- c(models$ModelPrefix, "Timestamp")
        predictTrait <<-  predictTrait
        write.csv(predictTrait, "predictions.csv", quote=FALSE, row.names = (FALSE))
      }else{
        predictTrait <- as.data.frame(read_csv("predictions.csv"), col_types = cols(), show_col_types = FALSE)
        predictTrait$Timestamp <- as.character(predictTrait$Timestamp)
        predictTrait <<-  predictTrait
      }

      output$predsPlot <- renderPlotly({
        print(input$traits)
        dat <- as.data.frame(read_csv(paste0(mainDir, "/", "predictions.csv"), col_types = cols(), show_col_types = FALSE))
        
        if (nrow(dat) != 0) {
          dat$No <- 1:nrow(dat)
          Trait <-dat[[input$traits]]
          dat$Timestamp <- reorder(dat$Timestamp, dat$No)
          p <- ggplot(data=dat, aes(y=Trait, x=No))+ geom_blank()+theme_classic()+ labs(y = input$traits)
          
          #p <- p+geom_point(aes(y=Trait, x=Timestamp, text=Sample))
          
          #p <- ggplot(data=dat, aes(y=Trait, x=reorder(Timestamp, No)))+ geom_blank()+theme_classic()+ labs(y = input$traits)
          
          p <- p+geom_point(aes(text=Sample))+
            scale_x_discrete("Timestamp", limits = dat$No, labels = dat$Timestamp)+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(axis.title.x = element_blank())
          #p
          p
        }else{
          p <- ggplot(data=dat, aes(y=dat[[input$traits]], x=Sample))+ geom_blank()+theme_classic()+ labs(y = input$traits)
          
          print("Empty dataframe")
          p
        }
        ggplotly(p, source = "predsPlot")
        
      })
      
      output$mytable = DT::renderDataTable({
        dat <- as.data.frame(read_csv(paste0(mainDir, "/", "predictions.csv")), col_types = cols(), show_col_types = FALSE)
      })
      
      
      if(4 %in% input$options){
        
        files <- list.files(pattern = "\\.csv$", recursive = TRUE)
        
        if(file.exists("predictions.csv")) {
          
          files <-files[-grep("predictions.csv", files)]
        }
        if(file.exists("merged.csv")) {
          
          files <-files[-grep("merged.csv", files)]
          showNotification("There already exists a merged file at this location. Stopping.", duration = 10)
          clicked(0)
          output$buttonPanel <- renderUI ({
            actionButton("button", "Start scanning")
          })
          #
        }else{
          
          output$info <- renderText({
            d <- event_data("plotly_click", source = "predsPlot")
            if (is.null(d)) {
              "Click events appear here (double-click to clear)"
            } else {
              x <- round(d$x, 2)
              y <- round(d$y, 2)
              #cat("[", x, ", ", y, "]", sep = "")
              chosen(x)
              print(paste0("You have selected sample number: ", chosen()))
            }
          })
          
          iter <- 1
          for (file in files){
            datBind <- as.data.frame(read_csv(file), show_col_types = FALSE)
            datBind$Names <- file
            datBind <- datBind[,c(ncol(datBind),1:(ncol(datBind)-1))]
            if (iter == 1){
              write.table(datBind, "merged.csv", quote=FALSE,  sep=",", row.names=FALSE)
              iter <- 2
            }else{
            write.table(datBind, "merged.csv", quote=FALSE, col.names=FALSE, append=TRUE, sep=",", row.names=FALSE)
            }
          }
          
          print("Merging files")
          showNotification(paste("Files have been merged. See file: merged.csv"), duration = 300)

        }
      }else{
        
        old_files <- read.table("fileList.txt")

        old_files <- tools::file_path_sans_ext(old_files$V1)
        base_old <- basename(old_files)
        if(1 %in% input$options){
          new_files <- setdiff(list.files(path=getwd(),pattern = "\\.spc$", recursive = TRUE, full.names = TRUE), old_files)
          new_files <- new_files[!sapply(new_files, function(x) any(str_detect(x, base_old)))]
          
        }else{
          new_files <- setdiff(list.files(path=getwd(),pattern = "\\.csv$", recursive = TRUE, full.names = TRUE), old_files)
          new_files <- new_files[!sapply(new_files, function(x) any(str_detect(x, base_old)))]
          
          if(file.exists("predictions.csv")) {
            
            new_files <-new_files[-grep("predictions.csv", new_files)]
          }
          
          if(file.exists("merged.csv")) {
            
            new_files <-new_files[-grep("merged.csv", new_files)] 
          }
          #new_files <-new_files[-grep("result", new_files)]
        }
        
        sapply(new_files, function(x) {
          if(1 %in% input$options){

            output$info <- renderText({
              d <- event_data("plotly_click", source = "predsPlot")
              if (is.null(d)) {
                "Click events appear here (double-click to clear)"
              } else {
                x <- round(d$x, 2)
                y <- round(d$y, 2)
                #cat("[", x, ", ", y, "]", sep = "")
                chosen(x)
                print(paste0("You have selected sample number: ", chosen()))
              }
            })
            
            
            print("Conversion selected")
            if(3 %in% input$options){
              print("Conversion with mean spectra selected")
              convertToCSV2(x)
            }else{
              convertToCSV(x)
            }
            
            pathNew <-paste0(file.path(mainDir, subDir), "/", newName)
            #pathNew <- here("converted", newName)
            if(!dir.exists(dirname(pathNew))) {
              dir.create(dirname(pathNew), recursive=TRUE)
              print("Created new folder")}
            write.csv(newSpec.spc, pathNew, row.names = FALSE)
            print("Converting new file...")
            
          }#end of conversion
          

    
          if(2 %in% input$options){ #if prediction was chosen
            #Load spectra from csv file
            #where x is a path to file
            print(x)
            
            output$info <- renderText({
              d <- event_data("plotly_click", source = "predsPlot")
              if (is.null(d)) {
                "Click events appear here (double-click to clear)"
              } else {
                x <- round(d$x, 2)
                y <- round(d$y, 2)
                #cat("[", x, ", ", y, "]", sep = "")
                chosen(x)
                print(paste0("You have selected sample number: ", chosen()))
              }
            })
            
            if(1 %in% input$options){
              #for pipepc conversion we need to re-create a path to .csv
              
              newName <- basename(x)
              newName <- paste0(substr(newName,1, nchar(newName)-4), ".csv")
              x <- newName
              newSpec.spc <- read_csv(paste0("./converted/",newName), col_types = cols(), show_col_types = FALSE)
              
            }else{
              newSpec.spc <- read_csv(x, col_types = cols(), show_col_types = FALSE)
              #If we want calculation on mean spectra from csv
              if(3 %in% input$options){
                print("Prediction with mean spectra selected")
                meanSP <- apply(newSpec.spc,2, mean)
                newSpec.spc <-  as.data.frame(t(meanSP))
              }
            }
            
            if (dim(newSpec.spc)[1]>1){ #Detect multispec files
              
              # output$info <- renderText({
              #   print("Did you mean to calculate mean spectra?")
              # })

              showNotification(paste("Multispectral files detected. Did you mean to calculate mean spectra?"), duration = 3)
              
            }else{
              allPreds <- c()
              allMD <- c()
              for (model in 1:length(models$Model)){
                
                spectraPred <- newSpec.spc
                #Do stuff with spectra. Return predictions (predicted object) and predicted Mahalanobis distance (predMD, it will be set to 0 if nonthing is calculated is calculated)
                predictSpectra(spectraPred, modelscript=models$Script[model])
                
                if (!exists(quote(predMD))){
                  predMD <- 0
                }
                pred <- predicted
                
                allPreds <- c(allPreds,pred)
                allMD <- c(allMD, predMD)
                
              }

              timePred <- as.character(Sys.time())
              predictTrait[1,] <- c(x, allPreds,  timePred)
              predictTrait[1,]$Timestamp <- timePred
              write.table(predictTrait[1,], "predictions.csv", quote=FALSE, col.names=FALSE, append=TRUE, sep=",", row.names=FALSE)
 
              
              res <- c(x, allMD) #Add name
              write.table(t(res), paste0("distanceMs.txt"), append=TRUE, quote = FALSE, row.names=FALSE, col.names = FALSE, sep = "\t")
              
              #print(res)
              }
            #print(head(newSpec.spc))
          }
          
          
        })
        
        old_files = c(old_files, new_files)
        if(length(old_files)>1){
          write.table(old_files, "fileList.txt", row.names=FALSE, col.names=FALSE, quote=TRUE)
          
        }
        

        # write.table(new_files, "fileList.txt", row.names=FALSE, col.names=FALSE, quote=TRUE)
        
        
        output$mytable = DT::renderDataTable({
          dat <- as.data.frame(read_csv(paste0(mainDir, "/", "predictions.csv")), col_types = cols(), show_col_types = FALSE)
        })
        
        
        output$distTable = DT::renderDataTable({
          dat <- as.data.frame(read.table(paste0(mainDir, "/", "distanceMs.txt"), header = TRUE, sep = "\t"))
        })
        
        output$predsPlot <- renderPlotly({
          #print(input$traits)
          dat <- as.data.frame(read_csv(paste0(mainDir, "/", "predictions.csv"), col_types = cols(), show_col_types = FALSE))
          
          if (nrow(dat) != 0) {
            datMs <- as.data.frame(read.table(paste0(mainDir, "/", "distanceMs.txt"), header = TRUE, sep = "\t"))
            
            datMs <- cbind(datMs$Sample, datMs[ , grepl(paste0(input$traits,"\\."), colnames(datMs) ) ])
            
            # datMs <- datMs %>%
            #   dplyr::select(matches(c("Sample", input$traits)))
            
            colnames(datMs) <- c("Sample", "MD", "LN")
             datMs$MS <- ifelse(datMs$LN >0.1 & datMs$MD > 2.5, 3, 1)
            datMs$MS <- ifelse(datMs$LN >0.05 & datMs$MD <= 2.5, 2, datMs$MS)
            datMs$MS <- factor(datMs$MS, levels = c(3, 2, 1))
            
            dat <- dat %>%
              dplyr::select(all_of(c("Sample", "Timestamp",input$traits)))
            
          
            dat <- left_join(dat, datMs)
            dat$No <- 1:nrow(dat)
            Trait <-dat[[input$traits]]
            dat$Timestamp <- reorder(dat$Timestamp, dat$No)
            cp <- c(`1` = "black", `2` = "blue", `3` = "red")
            dat$MS <- ordered(dat$MS, levels=c(1,2,3))
            p <- ggplot(data=dat, aes(y=Trait, x=No))+ geom_blank()+theme_classic()+ labs(y = input$traits)

            p <- p+geom_point(aes(colour = MS))+
              scale_colour_manual(values = cp)+
              geom_point(aes(colour = MS, shape = MS, size =  MS, text=
                               paste("No.", No, "<br>Name: ", Sample, "<br>Trait: ",
                                     Trait, "<br>MD: ", MD, "<br>LN: ", LN)),  tooltip = 'text')+
              scale_shape_manual(values=c(16, 1, 1))+
              scale_size_manual(values=c(1,4, 4))+
              theme(legend.position="none")+
              #guides(color = FALSE, size = FALSE, shape = FALSE)+
              scale_x_discrete("Timestamp", limits = dat$No, labels = dat$Timestamp)+
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(axis.title.x = element_blank())
            p
          }else{
            p <- ggplot(data=dat, aes(y=dat[[input$traits]], x=Sample))+ geom_blank()+theme_classic()+ labs(y = input$traits)
            
            print("Empty dataframe")
            p
          }
          ggplotly(p, source = "predsPlot", tooltip = "text")
          
        })
        
      }} 
    closeAllConnections()}
    
    )
  
  output$info <- renderText({
    d <- event_data("plotly_click", source = "predsPlot")
    if (is.null(d)) {
      "Click events appear here (double-click to clear)"
    } else {
      x <- round(d$x, 2)
      y <- round(d$y, 2)
      #cat("[", x, ", ", y, "]", sep = "")
      chosen(x)
      print(paste0("You have selected sample number: ", chosen()))
    }
  })
  
  output$spectra <- renderPlot({
    
    spec <- chosen()
    if(spec != 0 ){
      #print(spec)
      dat <- as.data.frame(read_csv(paste0(mainDir, "/", "predictions.csv"), col_types = cols(), show_col_types = FALSE))
      checked <- ( dat[spec,]$Sample)
      if(file.exists(checked)){
        #if(1 %in% input$options){
        spectra <- read_csv(paste0(checked), col_types = cols(), show_col_types = FALSE)
      }else{
        spectra <- read_csv(paste0(mainDir, "/converted/", checked), col_types = cols(), show_col_types = FALSE)
      }
      waves <- parse_number(colnames(spectra))
      matplot(waves, t(spectra), main = "Spectra", type = "l", lty=1,pch=".",xlab="Wavelength [nm]",ylab="Reflectance", col = "#0284db", bty = "l")
    }
  })
  
  chosen <- reactiveVal(0) 
  
  clicked <- reactiveVal(0)  
  
  observeEvent(input$button, {
    print(clicked)
    if(clicked()==0){
      clicked(clicked()+1)
      output$buttonPanel <- renderUI ({
        actionButton("button", "Stop scanning")
      })
    }else{
      clicked(clicked()-1)
      output$buttonPanel <- renderUI ({
        actionButton("button", "Start scanning")
      })
    }
    
  })
  
  output$buttonPanel <- renderUI ({
    actionButton("button", "Start scanning")
  })
  
  
  output$traitPanel <- renderUI ({
    
    selectInput("traits", choices = colnames(predictTrait)[-c(1, ncol(predictTrait))], label = "Select trait")
  })
  
  observeEvent(input$options, {
    clicked(0)
    output$buttonPanel <- renderUI ({
      actionButton("button", "Start scanning")
    })
    
    if(4 %in% input$options){
      
      output$optionsPanel <- renderUI ({
        
        checkboxGroupInput("options", label = h3("Options"), choices = list("Convert (.spc to .csv)" = 1, "Predict from NIRS" = 2, "Calculate mean spectra from a multispec file" = 3, "Merge files" = 4),
                           selected = c(4))
      })
    }
    savedPath <- pathToFolder$path
    options <- list(path= savedPath, inputOptions = input$options)
    saveRDS(options, "<folderpath>/caMeL/models/Options.RDS")
  })
  
  
  output$optionsPanel <- renderUI ({
    if (file.exists("options.RDS")){
      ops <- readRDS("<folderpath>/caMeL/models/Options.RDS")
      checkboxGroupInput("options", label = h3("Options"), choices = list("Convert (.spc to .csv)" = 1, "Predict from NIRS" = 2, "Calculate mean spectra from a multispec file" = 3, "Merge files" = 4),
                         selected = ops$inputOptions) # ops$inputOptions
    }else{
      checkboxGroupInput("options", label = h3("Options"), choices = list("Convert (.spc to .csv)" = 1, "Predict from NIRS" = 2, "Calculate mean spectra from a multispec file" = 3, "Merge files" = 4),
                         selected = c(1,2)) # ops$inputOptions
      
    }
    
  })

  session$onSessionEnded(function() {
    stopApp()
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
