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

library(EMSC)

cvPlot <- function(data_long, type = "RMSE", trait_vector, where = "bottom"){
  data_long <- data_long %>%
    dplyr::filter(Stats == type)
  labs <- c("None", bquote("SGolay" ~ 1^st ~"ord. derivative"), bquote("SGolay" ~ 2^nd ~"ord. derivative"))
  
  allWins <- attributes(data_long)["best"]$best
  allWins <- allWins[allWins$Trait %in% trait_vector,]
  data_long <- data_long[data_long$Trait %in% trait_vector,]
  model.labs <- c("PLS", "Elastic net")
  names(model.labs) <- c("pls", "elasnet")
  trait.labs <-  as.character(levels(data_long$Trait))
  names(trait.labs) <-  as.character(levels(data_long$Trait))
  s <-ggplot(data_long, aes(y=Value, x=Treatment, fill = Derivative))+
    
    geom_boxplot() +
    #scale_alpha_manual(values=c(0.2,1)[as.factor(winMode)])+
    facet_grid(Trait ~ Model, scales="free", 
               labeller = labeller(Model = model.labs, Trait = trait.labs))+
    xlab("Spectral pre-treatment") + ylab(type)+
    #scale_fill_manual(values=wes_palette(n=3, name="Darjeeling1"))+
    #scale_fill_manual(values =c(alpha("#0284db", 0.7), alpha("#ffd300", 0.7), alpha("#f00e11", 0.7)))+
    #scale_fill_manual(values=cols)+scale_color_manual(values=cols)
    scale_fill_brewer(palette="RdBu",name="Mathematical\ntreatment",
                      breaks=levels(data_long$Derivative),
                      labels=labs)+
    theme_classic()+
    
    theme(axis.text.x = element_text(size = 10, angle = 45,hjust=1),
          axis.text.y = element_text(size = 10),
          axis.title=element_text(size=15))
  
  #strip.background = element_rect(colour = "white", fill="lightgrey"), )
  
  get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  legendCV <<- get_legend(s)
  
  u= c(1:length(spTreat))[allWins$Treatment]#-0.24
  u <- u+c(-0.25, 0, 0.25)[allWins$Derivative]
  v <- s+
    theme(legend.position="none")
  #+ geom_rect(data = data.frame(Trait = allWins$Trait, Model = allWins$Model), aes(ymin = -Inf, ymax = Inf, xmin = u-0.12, xmax = u+0.13), alpha = 0.05, fill="blue", inherit.aes = FALSE)
  if (where == "top"){
    v+theme(axis.title.x = element_blank(), axis.text.x=element_blank())+ labs(y="RMSECV")
  }else{
    v+labs(y="RMSEP")
  }  
  
  
}
getLatexTable <- function(a, param){
  base <- a[,c(param, "wl")] 
  colnames(base)[1] <- "Y"
  base <- base %>%
    dplyr::filter(abs(Y) > 0.5)%>%
    arrange(desc(Y))%>%
    mutate_at(2, round, 0) %>%
    mutate_at(1, round, 2) %>%
    relocate("wl")
  colnames(base) <-c("Lambda", "Coefficient")
  base <- t(base)
  base <- rbind(formatC(as.numeric(base[1,]),format="f",digits=0),
                formatC(as.numeric(base[2,]),format="f",digits=2))
  row.names(base) <-c("Lambda", "Coefficient")
  
  cropWidth <-function(x){
    if(ncol(x)==20 ){
      if (ncol(x)==20){
        return(x)
      }
    }else if (ncol(x) >20){
      partA <- as.data.frame(x)[,1:20]
      partB <- as.data.frame(x)[,20:ncol(x)]
      colnames(partB) <- paste0("V", c(1:ncol(partB)))
      return(rbind(partA, cropWidth(partB)))
    }else if (ncol(x) <20){
      x <- cbind(x,c("0", "0"))
      colnames(x) <- paste0("V", c(1:ncol(x)))
      return(cropWidth(x))
      
    }
  }
  base <- cropWidth(base)
  base <- rbind(rep(0, 20),base)
  rownames(base)[1] <- param
  
  return(base)
}



library(chemometrics)

#custom diagPlot
diagPlot <- function (res, title = "", col = "black", pch = 16, labelOut = TRUE, id = 3) 
{
  
  labelDDMY <- function (x, y, 
                         id.n.sd = 3, id.n.od = 3, off = 0.03) 
  {
    xrange <- par("usr")
    yrange <- par("usr")[4]-par("usr")[3]
    xrange <- xrange[2] - xrange[1]
    if (id.n.sd > 0 && id.n.od > 0) {
      n <- length(x)
      ind.sd <- sort(x, index.return = TRUE)$ix
      ind.sd <- ind.sd[(n - id.n.sd + 1):n]
      ind.od <- sort(y, index.return = TRUE)$ix
      ind.od <- ind.od[(n - id.n.od + 1):n]
      sd.out <- as.integer(names(x[x>co.sd])) #leverage points
      od.out <- as.integer(names(y[y>co.od])) # orthogonal outliers
      ind.sdod <-intersect(sd.out,od.out) #bad leverage point
      sd.out <- sd.out[!(sd.out %in% ind.sdod)]
      od.out <- od.out[!(od.out %in% ind.sdod)]
      norm.out <- as.integer(names(x))[-c(sd.out, od.out, ind.sdod)]
      points(x[od.out], y[od.out], pch = 19, col = alpha("#ffd300", 0.7)) #yellow
      points(x[sd.out], y[sd.out], pch = 19, col = alpha("#878787", 0.7)) # grey
      points(x[norm.out], y[norm.out], pch = 19, col = alpha("#0284db", 0.7)) # blue
      points(x[ind.sdod], y[ind.sdod], pch = 19, col = alpha("#f00e11", 0.7)) #red #f00e11
      lab <- ind.od
      if (is.character(names(y))) {
        lab <- names(y[ind.od])
      }
      text(x[ind.od] + off * xrange, y[ind.od], lab)
      lab <- ind.sd
      if (is.character(names(x))) {
        lab <- names(x[ind.sd])
      }
      text(x[ind.sd] + off * xrange, y[ind.sd], lab)
      text(x = 0, y = co.od + 0.08, expression("C"[OD]), col = "darkgrey")
      text(x = co.sd + 0.03* xrange, y = 0, expression("C"[SD]), col = "darkgrey")
      
      legend(x=xrange*0.6, y = yrange*0.9, 
             legend = c("Regular observation", "Orthogonal point", "Good leverage point", "Bad leverage point"), 
             col = c(alpha("#0284db", 0.7), alpha("#ffd300", 0.7), alpha("#878787", 0.7), alpha("#f00e11", 0.7)),
             pch = 19, 
             cex = 0.9,
             inset = c(-0.45, 0),
             bty = "n", 
             y.intersp=1, #interline in the legend
             text.col = "black")
    }
  }
  sd <- res$sd
  od <- res$od
  co.sd <- res$cutoff.sd # add label
  co.od <- res$cutoff.od #add label
  if (is.null(pch)) {
    pch <- ifelse(is.null(col), 1, 19)
  }
  if (is.null(col)) {
    col <- "black"
  }
  #par(mar = c(5, 5, 5, 5))
  plot(sd, od, xlab = "Score distance", ylab = "Orthogonal distance", 
       main = title, pch = pch, col = "white", xlim = c(0, max(sd) * 
                                                          1.4), ylim = c(0, max(od) * 1.1), axes=F)
  box(bty="l")
  axis(2) #alpha("#0284db", 0.7)
  axis(1)
  abline(v = co.sd, col = "darkgrey", lty = "dashed") 
  abline(h = co.od, col = "darkgrey", lty = "dashed")
  if (labelOut) {
    labelDDMY(sd, od, id.n.sd = id, id.n.od = id)
  }
  
}
SNV<-function(spectra){                                             
  spectra<-as.matrix(spectra)
  
  spectrat<-t(spectra)
  spectrat_snv<-scale(spectrat,center=TRUE,scale=TRUE)
  spectra_snv<-t(spectrat_snv)
  return(spectra_snv)}

BASELINE <- function(spectra){
  #This pre-treatment uses hyperSpec package
  #Convert data to a hyperSpec S4 object:
  newHS<-new("hyperSpec", spc = spectra)
  #Compute baselines using order 2 polynomials:
  baseline<-spc.fit.poly.below(fit.to = newHS, poly.order = 2)
  #Baseline removal:
  newspectra<-newHS@data$spc-baseline@data$spc
  newspectra<-as.matrix(newspectra)
  return(newspectra)
}


DETREND <- function(spectra){
  newspectra<-pracma::detrend(t(spectra),tt="linear") # detrend
  newspectra <- as.matrix(t(newspectra))
  return(newspectra)
}

MSC <- function(spectra){
  newspectra <- pls::msc(as.matrix(spectra))
  refSpectra <- attributes(newspectra)["reference"]
  refSpectra <<- as.matrix(refSpectra$reference)
  return(newspectra)
}

EMSCmy <- function(spectra){
  refEMSC <<- EMSC(spectra, degree = 6)
  newspectra <- refEMSC$corrected
  return(newspectra)
}


addMSC <- function(spectra){
  newspectra <- mdatools::prep.msc(as.matrix(spectra), mspectrum = refSpectra)
  #print("Applying mean spectra reference")
  return(newspectra)
}

addEMSC <- function(spectra){
  newspectra <- predict(refEMSC, spectra)
  newspectra <- newspectra$corrected
  #print("Applying mean spectra reference")
  return(newspectra)
}

RMSEPmy<-function(predicted,observed,na.rm=FALSE)
{if(!na.rm)
  if(!(length(which(is.na(predicted)))+length(which(is.na(observed))))==0)
    stop("Datasets include NA values. To override this error, set na.rm=TRUE.")
  sqrt(sum((observed-predicted)^2,na.rm=TRUE)/
         length(which(!is.na(observed-predicted))))} #from chillR

RPDmy<-function(predicted,observed,na.rm=FALSE)
{if(!na.rm)
  if(!(length(which(is.na(predicted)))+length(which(is.na(observed))))==0)
    stop("Datasets include NA values. This may indicate a serious prediction problem. To override this error, set na.rm=TRUE.")
  sd(observed,na.rm=TRUE)/
    sqrt(sum((observed-predicted)^2,na.rm=TRUE)/
           length(which(!is.na(observed-predicted))))} #from chillR

MAE <-function(observed, predicted){ #mae - mean absolute error
  return(sum(abs(observed - predicted))/length(observed))   #tested against (ModelMetrics::mae(observed, predicted))
}

BIAS <- function(observed, predicted){
  return(sum(predicted-observed)/length(observed))
}

SEP <- function(observed, predicted) {  #standard error of prediction (SEP)
  b <- (predicted - observed - BIAS(observed, predicted))^2
  return(sqrt(sum(b)/(length(observed)-1)))
}

RPD.2 <- function(observed, predicted){
  return(sd(observed)/SEP(observed, predicted))
}

RER <- function(observed, predicted){ #RER, Range Error Ratio 
  range <- range(observed)[2]-range(observed)[1]
  return(range/SEP(observed, predicted))
}

#MEAN ABSOLUTE PERCENTAGE ERROR (MAPE)
MAPE = function(y_actual,y_predict){
  mean(abs((y_actual-y_predict)/y_actual))*100
}
rsq <- function(y_actual,y_predict){
  cor(y_actual,y_predict)^2
}


trainFolds <- function(dataSet, calSet, K, hyperGrid){
  
  foldsAll <<- splitFolds(K, dataSet, calSet)
  outTer <- buildModel(dataSet, calSet, K, hyperGrid)
  return(outTer)
}

getVal <- function(outTer){
  #Data from test validations
  df.V <- data.frame()
  for (elem in 1:K) { #change t to K after tests
    df.V<- rbind(df.V, outTer[[elem]]$val)}
  
  df.V <- df.V %>%
    gather(Stats, Value, R2:MAPE, factor_key=TRUE)
  
  # df.V <- df.V %>%
  #   mutate(ID = 1:nrow(df.V))
  # 
  return(df.V)
}

applySP <- function(grid, dataSet, mod = 1){
  
  library(pls)
  
  #Apply spectral pretreatments
  switch(as.character(grid),
         "none" = {
           #print("no spectral pretreatment")
           dataSet <- as.matrix(dataSet)
         },
         "msc" = {
           #print("msc")
           if (mod==1) {
             dataSet <- MSC(dataSet) # apply spectral transformations
           }else{
             dataSet <- addMSC(dataSet)
           }
           
         },
         "msc+detrend" = {
           #print("msc+detrend")
           if (mod==1){
             dataSet <- DETREND(MSC(dataSet)) # apply spectral transformations
           }else{
             dataSet <- DETREND(addMSC(dataSet))
           }
           #spectraPlot(dataSet, calibQual.fold[[param]])
         },
         "baseline" = {
           #print("baseline")
           dataSet <- BASELINE(dataSet)
           #matplot(t(dataSet))
         },
         "msc+baseline" = {
           #print("msc+baseline")
           if (mod==1){
             dataSet <- BASELINE(MSC(dataSet)) # apply spectral transformations
           }else{
             dataSet <- BASELINE(addMSC(dataSet))
           }
         },
         "snv" = {
           # #SNV
           #print("snv")
           dataSet <-SNV(dataSet) # Standard Normal Variate (SNV)
           #spectraPlot(dataSet, calibQual[[1]])
         },
         "snv+detrend" = {
           # #SNV
           #print("snv+detrend")
           dataSet <-DETREND(SNV(dataSet))
         },
         "snv+baseline" = {
           #print("snv+baseline")
           dataSet <-BASELINE(SNV(dataSet)) 
         },
         "EMSC" = {
           #print("EMSC")
           if (mod==1) {
             dataSet <- EMSCmy(dataSet)
           }else{
             dataSet <- addEMSC(dataSet)
           }
           
           #spectraPlot(dataSet, calibQual.fold[[param]])
         }
         
  )
  return(dataSet)
}


applyMT <- function(grid, dataSet){
  
  #Apply spectral pretreatments
  switch(as.character(grid),
         "none" = {
           #print("no spectral pretreatment")
           dataSet <- as.matrix(dataSet)
         },
         "2.7.2" = {
           #print("2.7")
           dataSet<-t(apply(dataSet, 1, FUN=signal::sgolayfilt, p = 2, n = 7, m = 2, ts = 1))
           #spectraPlot(dataSet, calibQual.fold[[param]])
         },
         "1.7.2" = {
           #dataSet <-calibSet.foldRem
           #print("1.7")
           dataSet<-t(apply(dataSet, 1, FUN=signal::sgolayfilt, p = 2, n = 7, m = 1, ts = 1))
           #spectraPlot(dataSet, calibQual.fold[[param]])
         },
         "0.7.2" = {
           #dataSet <-calibSet.foldRem
           #print("1.7")
           dataSet<-t(apply(dataSet, 1, FUN=signal::sgolayfilt, p = 2, n = 7, m = 0, ts = 1))
           #spectraPlot(dataSet, calibQual.fold[[param]])
         }
  )
  return(dataSet)
}

applyModel <- function(grid, df_train, mode = 1, pcs = 1, k = 1, nlv = 1, alpha = 0, lambda = 0){
  
  #Select model
  switch(as.character(grid),
         "pls" = {
           #print("pls")
           #PLS
           if (mode == 1){
             comps <- pcs
             fitModel <- pls::plsr(Y ~ X, data = df_train, validation = "CV", scale = TRUE) #val CV
            comps <- selectNcomp(fitModel, method = "onesigma")
             fitModel <- pls::plsr(Y ~ X, ncomp = comps, data = df_train, validation = "CV", scale = TRUE)
           
          }else{
             comps <- pcs
             fitModel <- pls::plsr(Y ~ X, ncomp = comps, data = df_train, validation = "none", scale = TRUE)
           }
           
           attributes(fitModel)$msc.ref<-refSpectra
           attributes(fitModel)$emsc.ref<-refEMSC
           attributes(fitModel)$tuning<-c(comps,0)
           
           obs <- df_train$Y
           
           preds <- fitModel$fitted.values[,,comps]
           
           r2 <- rsq(obs, preds)
           #rmsec <- sd(fitModel$residuals[,,comps])
          rmsec <- RMSEPmy(preds, obs)
           rpd <- RPDmy(obs,preds)
           rer <- RER(obs,preds)
           mape <- MAPE(obs,preds)
           descriptive <- c(r2, rpd, rmsec, rer, mape)
           names(descriptive) <- c("R2", "RPD", "RMSE", "RER", "MAPE")
           
         },
         "lwplsr" = {
           trait <- df_train$Y
           data <- df_train$X
           nlvdis <- 0
           h <- 2
           fitModel <- rchemo::lwplsr(data, trait, nlvdis = nlvdis, h = h, k = k, nlv = nlv, diss  = "mahal")
           attributes(fitModel)$msc.ref<-refSpectra
           attributes(fitModel)$emsc.ref<-refEMSC
           

           obs <- df_train$Y
           res <- predict(fitModel, data)
           preds <- res$pred

           r2 <- rsq(obs, preds)
           #rmsec <- sd(fitModel$residuals[,,comps])
           rmsec <- RMSEPmy(preds, obs)
           rpd <- RPDmy(obs,preds)
           rer <- RER(obs,preds)
           mape <- MAPE(obs,preds)
           descriptive <- c(r2, rpd, rmsec, rer, mape)
           names(descriptive) <- c("R2", "RPD", "RMSE", "RER", "MAPE")
         },
         "elasnet" = {
           #print("Elastic net")
           #Elastic net
           if (mode == 1){
             repCV10 <- trainControl(method = "repeatedcv", 
                                     number = 10, #cv
                                     returnResamp = "all", #all to see all results
                                     savePredictions = "final", 
                                     allowParallel = T, 
                                     verboseIter = F)
             
             
             fitModel <- train(Y ~ ., data = df_train,
                               method = 'glmnet', 
                               preProc = c("center", "scale"),
                               tuneGrid=expand.grid(
                                 .alpha=seq(0, 1, by = 0.1), #alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.
                                 .lambda=seq(0, 1, by = 0.05)),  
                               trControl = repCV10) 
             
           }else{
             repCV10 <- trainControl(returnResamp = "all", #all to see all results
                                     savePredictions = "final", 
                                     allowParallel = F, 
                                     verboseIter = F)
             #check if the order is correct
             fitModel <- train(Y ~ ., data = df_train,
                               method = 'glmnet', 
                               preProc = c("center", "scale"),
                               tuneGrid=expand.grid(
                                 .alpha=alpha, #alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.
                                 .lambda=lambda),  
                               trControl = repCV10) 
             
           }
           
           
           #plot(fitModel)
           #fitModel$bestTune
           attributes(fitModel)$msc.ref<-refSpectra
           attributes(fitModel)$emsc.ref<-refEMSC
           
           best(fitModel$results, metric = 'RMSE', maximize = F) # indicate the best model
           attributes(fitModel)$tuning<-unlist(fitModel$bestTune)
           obs <- fitModel$pred$obs
           preds <- fitModel$pred$pred
           r2 <- rsq(obs, preds)
           rmsec <- RMSEPmy(preds, obs)
           rpd <- RPDmy(obs, preds)
           rer <- RER(obs,preds)
           mape <- MAPE(obs,preds)
           descriptive <- c(r2, rpd, rmsec, rer, mape)
           names(descriptive) <- c("R2", "RPD", "RMSE", "RER", "MAPE")
         }
         
  )#end of model switch
  return(list(model = fitModel, stats = descriptive, preds = preds, obs = obs))
}

buildModel <- function(dataSet, calSet, t = K, hyperGrid){
  hyperGrid0 <- hyperGrid
  outTer <- list()
  calibVal = colnames(calSet)
  
  if("pls" %in% hyperGrid0$Model){
    outTerPLS <- list()

    hyperGrid1 <- expand.grid(unique(hyperGrid0$Trait),unique(hyperGrid0$Treatment), unique(hyperGrid0$Derivative), unique(hyperGrid0$Model)[unique(hyperGrid0$Model)=="pls"], pcs)
    colnames(hyperGrid1) <- c("Trait", "Treatment", "Derivative", "Model", "PC")
    hyperGrid <- arrange(hyperGrid1, Trait)
    v_mod1 <- nrow(hyperGrid)
    
    for (fold in 1:K){ 
      
      resVal <- list()
      for (v in 1:v_mod1){
        #print(v)
        #See what goes where
        #start fresh
        dataSet.fold <- dataSet[!is.na(calSet[,hyperGrid[v,]$Trait]),]
        calSet.fold <- as.data.frame(calSet[!is.na(calSet[,hyperGrid[v,]$Trait]),])
        colnames(calSet.fold) <- hyperGrid[v,]$Trait
        focus <- foldsAll[[hyperGrid[v,]$Trait]]
        trainSet <- which(focus != fold)
        testSet <- setdiff(1:length(focus), trainSet)
        
        valSet.fold<- dataSet.fold[testSet,]
        valQual.fold <- calSet.fold[testSet,,drop = FALSE]
        
        dataSet.fold <- dataSet.fold[trainSet,]
        calSet.fold <- calSet.fold[trainSet,,drop=FALSE]
        
        traitCalib <- calSet.fold[,hyperGrid[v,]$Trait]
        traitRem <- valQual.fold[,hyperGrid[v,]$Trait]
        #Remove outliers
        dataSet.foldRem <- dataSet.fold[!(abs((traitCalib - mean(traitCalib))/sd(traitCalib)) > 3),] #farther than 3 SD from the mean
        traitCalib <- traitCalib[!(abs((traitCalib - mean(traitCalib))/sd(traitCalib)) > 3)]
        
        MSC(dataSet.foldRem)
        EMSCmy(dataSet.foldRem)
        
        dataSet.combo <- applySP(hyperGrid[v,2], dataSet.foldRem)
        #spectraPlot(dataSet.combo, trait)
        #Apply math treatment
        dataSet.combo <- applyMT(hyperGrid[v,3], dataSet.combo)
        
        param <- calibVal[hyperGrid[v,]$Trait]
        
        df_train <- data.frame(Y=traitCalib, X=I(dataSet.combo))
        alpha <- hyperGrid[v,5]
        model <- applyModel(hyperGrid[v,4], df_train, mode=2, pcs = alpha)
        fitModel <- model$model
        
        #Apply spectral pretreatments
        #1.Extract refspectrum for MSC/EMSC, in case needed
        refSpectra <- attributes(fitModel)["msc.ref"]
        refSpectra <- as.matrix(refSpectra$msc.ref)
        refEMSC <- attributes(fitModel)["emsc.ref"]
        refEMSC <- refEMSC$emsc.ref
        
        #2.Apply spectral transformations
        valSet <- applySP (hyperGrid[v,2], valSet.fold, mod = 2)
        valSet <- applyMT (hyperGrid[v,3], valSet)
        
        #run model
        if(!hyperGrid[v,4]=="plsLOCO"){
          df_test <- data.frame(X=I(valSet))
          obs <- traitRem
          
          preds <- predict(fitModel, newdata = df_test)
        }else{ #if we go the plsLOCO route
          calSet <- applySP (hyperGrid[v,2], dataSet.foldRem, mod = 1)
          calSet <- applyMT (hyperGrid[v,3], calSet)
          
          preds <- plsLOCO(calSet, traitCalib, valSet, traitRem)
          obs <- traitRem
        }
        
        if(length(dim(preds))==3){ #if it is an array from PLS package
          rpd <- RPDmy(preds[,,dim(preds)[3]], obs) #to accomodate for results from PLS
          rmse <- RMSEPmy(preds[,,dim(preds)[3]], obs) 
          preds <- preds[,,dim(preds)[3]]
          
        }else{
          rmse <- RMSEPmy(preds, obs) 
          rpd <- RPDmy(preds, obs) 
        }
        rer <- RER(obs,preds)
        mape <- MAPE(obs,preds)
        r2 <- rsq(obs, preds)
        descriptive <- c(r2, rpd, rmse, rer, mape)
        names(descriptive) <- valStats
        
        resVal[[v]] <- list(stats = descriptive, obs = obs, preds = preds)
      }
      
      fitStats.V <- hyperGrid
      for (fit in 1:length(names(descriptive))){
        #print(names(descriptive)[fit])
        df <- as.data.frame(sapply(resVal, function(x) rbind(x$stats[[valStats[fit]]])))
        colnames(df) <- names(descriptive)[fit]
        fitStats.V <- cbind(fitStats.V, df)
      }
      
      print(paste0("Finishing fold ", fold))
      #outTer[[fold]] <- list(cv = fitStats, val = fitStats.V, all = results, all.V = resVal)
      outTerPLS[[fold]] <- list(val = fitStats.V)
      saveRDS(outTerPLS, "./data.outTerPLS.RDS")
    }
    outTer[["PLS"]] <- outTerPLS
    
  }
  
  if("elasnet" %in% hyperGrid0$Model){
    alphaL=seq(0, 1, by = 0.1) #alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.
    #lambdaL=seq(0, 0.1, by = 0.1)
    lambdaL=seq(0, 1, by = 0.1)
    hyperGrid2 <- expand.grid(unique(hyperGrid0$Trait),unique(hyperGrid0$Treatment), unique(hyperGrid0$Derivative), unique(hyperGrid0$Model)[unique(hyperGrid0$Model)=="elasnet"], alphaL, lambdaL)
    colnames(hyperGrid2) <- c("Trait", "Treatment", "Derivative", "Model", "Alpha", "Lambda")
    hyperGrid <- arrange(hyperGrid2, Trait)
    
    outTerEN <- list()
    v_mod2 <- nrow(hyperGrid)
    
    for (fold in 1:K){ 
      
      resVal <- list()
      for (v in 1:v_mod2){
        #print(v)
        #See what goes where
        #start fresh
        dataSet.fold <- dataSet[!is.na(calSet[,hyperGrid[v,]$Trait]),]
        calSet.fold <- as.data.frame(calSet[!is.na(calSet[,hyperGrid[v,]$Trait]),])
        colnames(calSet.fold) <- hyperGrid[v,]$Trait
        focus <- foldsAll[[hyperGrid[v,]$Trait]]
        trainSet <- which(focus != fold)
        testSet <- setdiff(1:length(focus), trainSet)
        
        valSet.fold<- dataSet.fold[testSet,]
        valQual.fold <- calSet.fold[testSet,,drop = FALSE]
        
        dataSet.fold <- dataSet.fold[trainSet,]
        calSet.fold <- calSet.fold[trainSet,,drop=FALSE]
        
        traitCalib <- calSet.fold[,hyperGrid[v,]$Trait]
        traitRem <- valQual.fold[,hyperGrid[v,]$Trait]
        #Remove outliers
        dataSet.foldRem <- dataSet.fold[!(abs((traitCalib - mean(traitCalib))/sd(traitCalib)) > 3),] #farther than 3 SD from the mean
        traitCalib <- traitCalib[!(abs((traitCalib - mean(traitCalib))/sd(traitCalib)) > 3)]
        
        MSC(dataSet.foldRem)
        EMSCmy(dataSet.foldRem)
        
        dataSet.combo <- applySP(hyperGrid[v,2], dataSet.foldRem)
        #spectraPlot(dataSet.combo, trait)
        #Apply math treatment
        dataSet.combo <- applyMT(hyperGrid[v,3], dataSet.combo)
        
        param <- calibVal[hyperGrid[v,]$Trait]
        
        df_train <- data.frame(Y=traitCalib, X=I(dataSet.combo))
        alpha <- hyperGrid[v,5]
        lambda <- hyperGrid[v,6]
        model <- applyModel(hyperGrid[v,4], df_train, mode=2, alpha=alpha, lambda=lambda)
        fitModel <- model$model
        
        #Apply spectral pretreatments
        #1.Extract refspectrum for MSC/EMSC, in case needed
        refSpectra <- attributes(fitModel)["msc.ref"]
        refSpectra <- as.matrix(refSpectra$msc.ref)
        refEMSC <- attributes(fitModel)["emsc.ref"]
        refEMSC <- refEMSC$emsc.ref
        
        #2.Apply spectral transformations
        valSet <- applySP (hyperGrid[v,2], valSet.fold, mod = 2)
        valSet <- applyMT (hyperGrid[v,3], valSet)
        
        #run model
        if(!hyperGrid[v,4]=="plsLOCO"){
          df_test <- data.frame(X=I(valSet))
          obs <- traitRem
          
          preds <- predict(fitModel, newdata = df_test)
        }else{ #if we go the plsLOCO route
          calSet <- applySP (hyperGrid[v,2], dataSet.foldRem, mod = 1)
          calSet <- applyMT (hyperGrid[v,3], calSet)
          
          preds <- plsLOCO(calSet, traitCalib, valSet, traitRem)
          obs <- traitRem
        }
        
        if(length(dim(preds))==3){ #if it is an array from PLS package
          rpd <- RPDmy(preds[,,dim(preds)[3]], obs) #to accomodate for results from PLS
          rmse <- RMSEPmy(preds[,,dim(preds)[3]], obs) 
          preds <- preds[,,dim(preds)[3]]
          
        }else{
          rmse <- RMSEPmy(preds, obs) 
          rpd <- RPDmy(preds, obs) 
        }
        rer <- RER(obs,preds)
        mape <- MAPE(obs,preds)
        r2 <- rsq(obs, preds)
        descriptive <- c(r2, rpd, rmse, rer, mape)
        names(descriptive) <- valStats
        
        resVal[[v]] <- list(stats = descriptive, obs = obs, preds = preds)
      }
      
      fitStats.V <- hyperGrid
      for (fit in 1:length(names(descriptive))){
        #print(names(descriptive)[fit])
        df <- as.data.frame(sapply(resVal, function(x) rbind(x$stats[[valStats[fit]]])))
        colnames(df) <- names(descriptive)[fit]
        fitStats.V <- cbind(fitStats.V, df)
      }
      
      print(paste0("Finishing fold ", fold))
      #outTer[[fold]] <- list(cv = fitStats, val = fitStats.V, all = results, all.V = resVal)
      outTerEN[[fold]] <- list(val = fitStats.V)
      saveRDS(outTerEN, "./data.outTerEN.RDS")
    }
    outTer[["elasnet"]] <- outTerEN
    
  }
  

  
  if("lwplsr" %in% hyperGrid0$Model){
    outTerLWPLSR <- list()

    k <- c(10,20,30,40,50) # for tuning
    nlv <- 1:10 #for tuning
    hyperGrid3 <- expand.grid(unique(hyperGrid0$Trait),unique(hyperGrid0$Treatment), unique(hyperGrid0$Derivative), unique(hyperGrid0$Model)[unique(hyperGrid0$Model)=="lwplsr"], k ,nlv)
    colnames(hyperGrid3) <- c("Trait", "Treatment", "Derivative", "Model", "K", "NLV")
    hyperGrid <- arrange(hyperGrid3, Trait)
    v_mod3 <- nrow(hyperGrid)
    
    for (fold in 1:K){ 

      resVal <- list()
      for (v in 1:v_mod3){
        #print(v)
        #See what goes where
        #start fresh
        dataSet.fold <- dataSet[!is.na(calSet[,hyperGrid[v,]$Trait]),]
        calSet.fold <- as.data.frame(calSet[!is.na(calSet[,hyperGrid[v,]$Trait]),])
        colnames(calSet.fold) <- hyperGrid[v,]$Trait
        focus <- foldsAll[[hyperGrid[v,]$Trait]]
        trainSet <- which(focus != fold)
        testSet <- setdiff(1:length(focus), trainSet)
        
        valSet.fold<- dataSet.fold[testSet,]
        valQual.fold <- calSet.fold[testSet,,drop = FALSE]
        
        dataSet.fold <- dataSet.fold[trainSet,]
        calSet.fold <- calSet.fold[trainSet,,drop=FALSE]
        
        traitCalib <- calSet.fold[,hyperGrid[v,]$Trait]
        traitRem <- valQual.fold[,hyperGrid[v,]$Trait]
        #Remove outliers
        dataSet.foldRem <- dataSet.fold[!(abs((traitCalib - mean(traitCalib))/sd(traitCalib)) > 3),] #farther than 3 SD from the mean
        traitCalib <- traitCalib[!(abs((traitCalib - mean(traitCalib))/sd(traitCalib)) > 3)]
        
        MSC(dataSet.foldRem)
        EMSCmy(dataSet.foldRem)
        
        dataSet.combo <- applySP(hyperGrid[v,2], dataSet.foldRem)
        #spectraPlot(dataSet.combo, trait)
        #Apply math treatment
        dataSet.combo <- applyMT(hyperGrid[v,3], dataSet.combo)
        
        param <- calibVal[hyperGrid[v,]$Trait]
        
        df_train <- data.frame(Y=traitCalib, X=I(dataSet.combo))
        k <- hyperGrid[v,5]
        nlv <- hyperGrid[v,6]
        model <- applyModel(hyperGrid[v,4], df_train, mode=2, k = k, nlv = nlv)
        fitModel <- model$model
        
        #Apply spectral pretreatments
        #1.Extract refspectrum for MSC/EMSC, in case needed
        refSpectra <- attributes(fitModel)["msc.ref"]
        refSpectra <- as.matrix(refSpectra$msc.ref)
        refEMSC <- attributes(fitModel)["emsc.ref"]
        refEMSC <- refEMSC$emsc.ref
        
        #2.Apply spectral transformations
        valSet <- applySP (hyperGrid[v,2], valSet.fold, mod = 2)
        valSet <- applyMT (hyperGrid[v,3], valSet)
        
        #run model

          df_test <- data.frame(X=I(valSet))
          obs <- traitRem
          
          preds <- predict(fitModel, df_test)
          preds <- preds$pred

          rmse <- RMSEPmy(preds, obs) 
          rpd <- RPDmy(preds, obs) 

        rer <- RER(obs,preds)
        mape <- MAPE(obs,preds)
        r2 <- rsq(obs, preds)
        descriptive <- c(r2, rpd, rmse, rer, mape)
        names(descriptive) <- valStats
        
        resVal[[v]] <- list(stats = descriptive, obs = obs, preds = preds)
      }
      
      fitStats.V <- hyperGrid
      for (fit in 1:length(names(descriptive))){
        #print(names(descriptive)[fit])
        df <- as.data.frame(sapply(resVal, function(x) rbind(x$stats[[valStats[fit]]])))
        colnames(df) <- names(descriptive)[fit]
        fitStats.V <- cbind(fitStats.V, df)
      }
      
      print(paste0("Finishing fold ", fold))
      #outTer[[fold]] <- list(cv = fitStats, val = fitStats.V, all = results, all.V = resVal)
      outTerLWPLSR[[fold]] <- list(val = fitStats.V)
      saveRDS(outTerLWPLSR, "./data.outTerLWPLSR.RDS")
    }
    outTer[["LWPLSR"]] <- outTerLWPLSR
    
  }
  
  

  return(outTer)
}

hyperTune <- function(dataSet, qualSet, valGrid){
  valGrid$PC <- NA
  valGrid$RMSE <- NA
  v_mod <- nrow(valGrid)
  for(models in 1:v_mod){
    calibSet.fold <- dataSet
    calibQual.fold <- qualSet
    
    #set up ref spectra for MSC and EMSC per fold, in case it is needed
    MSC(calibSet.fold)
    EMSCmy(calibSet.fold)
    
    param <- calibVal
    trait <- calibQual.fold[[param]]
    
    #remove outliers
    calibSet.foldRem <- calibSet.fold[!(abs((trait - mean(trait))/sd(trait)) > 3),] #farther than 3 SD from the mean
    trait <- trait[!(abs((trait - mean(trait))/sd(trait)) > 3)]
    
    calibSet.combo <- applySP(valGrid[models,2], calibSet.foldRem)
    #spectraPlot(calibSet.combo, trait)
    #Apply math treatment
    calibSet.combo <- applyMT(valGrid[models,3], calibSet.combo)
    
    df_train <- data.frame(Y=trait, X=I(calibSet.combo))
    
    res <- applyModel(valGrid[models,4], df_train)
    stats <- res$stats[3] #RMSE
    valGrid[models,]$PC <- res$model$ncomp
    valGrid[models,]$RMSE <- stats
    valGrid <- valGrid %>%
      arrange(RMSE)

  }
  return(valGrid)
}

findBest <- function(hyperGrid){
  valGrid <- hyperGrid
  max <- valGrid[1,]$RMSE+valGrid[1,]$RMSE*0.05
  valGrid2 <- valGrid %>%
    filter(RMSE < max) %>%
    arrange(PC)
  return(valGrid2)
}



splitFolds <- function(n,dataframe, quals){
  foldsAll <- list()
  calibVal <- colnames(quals)

  if (any(is.na(quals[]))){
    print("There are missing values in the quality dataset.")
    
    for (p in calibVal){
      #print(p)
      N <- nrow(dataframe[!is.na(quals[,p]),])
      print(paste0(p, ": ", N, " complete rows"))
      folds <- rep( 1:n, ceiling(N/n) )
      folds <- sample(folds) 
      folds <- folds[1:N] 
      foldsAll[[p]] <- folds
    }
  }else{
    N <- nrow(dataframe)
    #print(N)
    folds <- rep(1:n, ceiling(N/n) )
    folds <- sample(folds) 
    folds <- folds[1:N] 
    for (p in calibVal){
      foldsAll[[p]] <- folds
    }
  }
  return(foldsAll)
}

fig_label <- function(text, region="plot", pos="topleft", cex=NULL, ...) {
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  text(x1+100, y1, text, cex=2, font = 2, ...)
  return(invisible(c(x,y)))
}

resPlot <- function(results,bestModels, fold){
  par(mfrow = c(2,2))
  #plot results
  for (it in 1:length(calibVal)){
    bestModel <- bestModels[bestModels$Trait ==calibVal[it],]
    bestModelID <- bestModel[1,]$ID
    grid <- testGrid[testGrid$Treatment==bestModel$Treatment[1] 
                     & testGrid$Derivative==bestModel$Derivative[1] & testGrid$Model==bestModel$Model[1],]
    bestModelID.cv <- as.integer(rownames(grid))
    grid <- valGrid[valGrid$Trait ==calibVal[it] & valGrid$Treatment==bestModel$Treatment[1] 
                    & valGrid$Derivative==bestModel$Derivative[1] & valGrid$Model==bestModel$Model[1],]
    bestModelID.cv2 <- as.integer(rownames(grid))
    #bestModelID.cv <- which(rownames(valGrid[valGrid$Trait==calibVal[it],])==bestModelID)
    x <- results[[fold]]$all[[it]]
    #preds <- results[[fold]]$all[[calibVal[it]]]$allMods[[bestModelID]]$model$pred$pred
    
    preds <- attributes(x$allMods[[bestModelID.cv]]$model)["preds"]$preds
    obs <- attributes(x$allMods[[bestModelID.cv]]$model)["obs"]$obs
    
    r2 <- rsq(obs, preds) #from cv
    #rpd <- mean(bestModel[bestModel$Stats=="RPD",]$Value) #mean rpd across folds
    rpd <- RPDmy(results[[fold]]$all.V[[bestModelID.cv2]]$obs, results[[fold]]$all.V[[bestModelID.cv2]]$preds) #from validation
    rmsec <- RMSEPmy(preds, obs)
    #rmsep <- mean(bestModel[bestModel$Stats=="RMSE",]$Value) #mean rmsec across folds
    rmsep <- results[[fold]]$val$RMSE[[bestModelID.cv2]]
    xlim = c(min(obs)*0.95, max(obs) *1.05)
    ylim = c(min(preds)*0.95, max(preds)*1.05)
    plot(obs, preds, pch = 19,
         main=paste0(calibVal[it], " (", paste0(as.vector(bestModel$Treatment[1]), "+", as.vector(bestModel$Derivative[1]), "+", as.vector(bestModel$Model[1]), ")")),
         axes=F, col = "white", xlab = "Observations", ylab = "Predictions", xlim = xlim, ylim=ylim)
    
    box(bty="l")
    xmin <- par("usr")[1]
    xmax <- par("usr")[2]
    ymin <- par("usr")[3]
    ymax <- par("usr")[4]
    abline(coef = c(0,1), col = "darkgrey")#, lty = "dashed")
    axis(2)
    axis(1)
    points(obs, preds,   col = alpha("658296", 0.6), pch = 19) # blue "0284db")
    points(results[[fold]]$all.V[[bestModelID.cv2]]$obs, results[[fold]]$all.V[[bestModelID.cv2]]$preds, col = alpha("red", 1), pch = 19)
    rp = vector('expression',4)
    
    #rpdC <- RPD.2(predict(fitModel, df_train), trait)
    rp[1] = substitute(expression(R^2 == MYVALUE),
                       list(MYVALUE = format(r2,dig=2)))[2]
    
    rp[2] = substitute(expression(RMSEC == MYVALUE),
                       list(MYVALUE = format(rmsec,dig=3)))[2] 
    rp[3] = substitute(expression(RMSEP == MYVALUE),
                       list(MYVALUE = format(rmsep,dig=3)))[2]
    rp[4] = substitute(expression(RPD == MYVALUE),
                       list(MYVALUE = format(rpd,dig=3)))[2]
    # legend(x=xmin, y=ymax, legend =  rp, cex = 0.9, bty = "n",
    #        y.intersp=0.2) #interline in the legend)
    legend(x=xmin, y=ymax, legend =  rp, cex = 0.9, bty = "n", y.intersp=0.8) #interline in the legend)
    # legend("topleft", legend =  rp, cex = 0.9, bty = "n", y.intersp=0.2) #interline in the legend)
    
    legend(x= range(obs)[2]-(range(obs)[2]-range(obs)[1])/4, y=mean(range(preds))-(mean(range(preds))-min(range(preds)))/3,
           legend = c("Training set", "Test set"),
           col = c(alpha("658296", 0.7), alpha("red", 1)),
           pch = 19,
           cex = 0.9, bty = "n",
           y.intersp=0.8)
  }
  
  par(mfrow = c(1,1))
}
