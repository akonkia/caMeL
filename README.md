caMeL - chemometrics applied machine learning interface

Notes on the caMeL app:

If you want to use the caMeL app, there are the following steps  [(-) steps are optional]:
- download the full caMeL repository and folder with R 4.2.0 to your drive of choice (url: https://cran.r-project.org/bin/windows/base/old/4.2.0/).
- remove file Options.RDS, if present, as it may contain old settings (it is in the /models folder).
- the code in app.R, run.R and camel.bat needs to be changed to reflect the <filepath>. Go for an absolute path.
- move the caMeL shortcut file to your Desktop.
(-) if the icon does not present a camel, change the shortcut icon to camelapha.ico by going into "Properties" of the shortcut and changing icon path.
(-) if the terminal stays maximized after the app launches, you can change this behaviour by going into shortcut>"Properties" (set to Minimise window).

To start the app, double-click the shortcut and wait.
When the app is pointed to a folder to scan for files, it will create:
- a /converted folder,
- prediction.csv file,
- distanceMs file,
- fileList.txt file.

How the app works:
- you use the shortcut icon to launch a .bat file (works on Windows, create a .command file on macOS).
- a Shiny R app will be launched in your browser.
- point the caMeL app to a folder of your choice and start scanning for .csv of .spc files.
- change the app settings if you want to merge files.

If something crashes for some reason in the middle and the app breaks, it might be wise to remove the /converted/ folder and the three files to start anew. This applies if you want to re-do the analysis for any folder of choice.

If you choose the file merging option, the app will create a merged.csv file. You cannot re-run the merging for as long as the merged.csv will be present in the folder. Delete the merged.csv or move it to be able to run again.
If something crashes, the run.Rout file should have the printed messages to troubleshoot.
The /lookhere folder contains .spc files that were used for testing and demonstration purposes.
The distanceMs file contains Mahalanobis distance measure (MD) for each sample as compared to the calibration subset used to create the predictive model. LN is the distance to the closest local neighbour. If the sample does not have close neighbours, it will be marked with a blue dot and halo. It can be picked up to extend the calibration. If the sample has a very distant MD score and no local neighbours, it will be flagged with red color and a halo.


This app works best with R 4.2.0. Newer R releases destroyed several dependencies in packages used by the app. If you have R, but the app still does not work, install Rstudio and try running the app.R script in RStudio to troubloeshoot.

This repository contains one model, predicting dry matter content (in %) - DM. This model was built to work with .spc output from Corona Zeiss Extreme NIRS spectrometer. You can test the app on supplied .spc files in the "look here" folder and make modifications as you need.

Want to create you own models and use them with caMeL?
- caMeL calls external scripts to run the models. It was done to allow for more freedom in creating models (you can change this behaviour with a little programming experience to launch models globally - this will speed up the process).
- a model is called by a function predictSpectra with two parameters: spectraPred (an object holding the spectra of a singular sample, read from the .spc or .csv file), modelscript (a pointer to path where the script to be sourced is placed). These two objects are created by the app for your. You need to provide the logic to return the predictions. Optionally, return the calculations for predMD as well (mahalanobis based calculations to check if your new sample is similar to the samples used in the model training).
- use createModelsFile.R to create a .txt file that will point the app to the models you want to use. It is necessary to follow the schema: you need the model prefix and the script that will be called by the app to perform calculation.
- in your model script, contain all information necessary to return predictions (a "predicted" object). See DM.R for an example of a callable model compatible with caMeL.

The caMeL app was created as part of a secondment funded from the European Union's Horizon 2020 research and innovation program under the Marie Sklodowska-Curie grant agreement No. 841882.
