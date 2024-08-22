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

#create settings file
settings <- data.frame(matrix(ncol=2,nrow=1, dimnames=list(NULL, c("ModelPrefix", "Script"))))

settings[1,] <- c("DM", "/models/DM.R")

#Add more models in new rows. for example:
#settings <- data.frame(matrix(ncol=2,nrow=6, dimnames=list(NULL, c("ModelPrefix", "Script"))))
#settings[1,] <- c("DM", "/models/DM.R")
#settings[2,] <- c("YourModelPrefix", "pathToYourModelScript")


write.table(settings, "./models/modelSettings.txt")

