#create settings file
settings <- data.frame(matrix(ncol=2,nrow=1, dimnames=list(NULL, c("ModelPrefix", "Script"))))

settings[1,] <- c("DM", "/models/DM.R")

#Add more models in new rows. for example:
#settings <- data.frame(matrix(ncol=2,nrow=6, dimnames=list(NULL, c("ModelPrefix", "Script"))))
#settings[1,] <- c("DM", "/models/DM.R")
#settings[2,] <- c("YourModelPrefix", "pathToYourModelScript")


write.table(settings, "./models/modelSettings.txt")

