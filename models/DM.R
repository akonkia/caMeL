list.of.packages <- c("glmnet", "caret", "signal", "mdatools", "chemometrics","EMSC","pls")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library('hyperSpec')
library(caret)
library(pls)
source('<folderpath>/funcsFQual5.R')

getMD <- function(pcamodel, newSample){
  spc.x <- pcamodel$x
  perc_var <- cumsum(pcamodel$sdev^2/sum(pcamodel$sdev^2))
  pecs <- grep(0.99, perc_var)[1]
  
  #Project the new sample onto the PCA space
  newpoint <- scale(newSample, pcamodel$center, pcamodel$scale) %*% pcamodel$rotation
  newMah <- rbind(spc.x,newpoint) 
  new.mah <- chemometrics::Moutlier(newMah[,1:pecs],quantile =0.975, plot=FALSE) 
  newMD <- new.mah$md[length(new.mah$md)]
  
  oldMD <- new.mah$md[-length(new.mah$md)]
  #Find local friends
  closestDist <- min(abs(newMD-oldMD)) 
  res <- c(newMD, closestDist)
  res <- round(res, 2)
  #return(list(MD = newMD, LN = closestDist))
  return(res)
}

# # #Loop over models
# # allModels <- base::apply(models,1, function(x) {
# #   #library(tidyverse)
# #   readRDS(paste0("/Users/akonkia/Desktop/caMeL-backup-Teagasc/caMeL/models/", unname(x[1])))
# # })

modelPred <- readRDS(paste0("/Users/akonkia/Desktop/caMeL-backup-Teagasc/caMeLmodelFree/", "models/DM.RDS"))
res <- modelPred
name <- res$name
fitModel <<- res$model
#extract EMSC, MSC, pca, combo from fitModel
refSpectra <- attributes(fitModel)["msc.ref"]
refSpectra <<- as.matrix(refSpectra$msc.ref)
refEMSC <- attributes(fitModel)["emsc.ref"]
refEMSC <<- refEMSC$emsc.ref
combo <- res$combo
combo <<- combo
spc.pca <- res$pca
spc.pca <<-spc.pca
#predict on new spectra
df_test <- data.frame(X=I(as.matrix(applyMT(combo$Derivative, applySP(combo$Treatment, newSpec.spc, mod = 2)))))
df_test <<- df_test
if(combo$Model =="pls"){
  preds <- predict(fitModel, newdata = df_test)
  preds <- preds[,,dim(preds)[3]]
}else  if(combo$Model =="elasnet"){
  preds <- predict(fitModel, newdata = df_test)

}else if(combo$Model == "lwplsr"){
  preds <- predict(fitModel,  df_test)
  preds <- preds$pred
}
predMD <<-  getMD(spc.pca,newSpec.spc)
predicted <<- round(preds,2)
