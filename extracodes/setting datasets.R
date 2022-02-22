
#1. First example -----------------------------------------------------------

x<-read.table("D:/Cloud/Dropbox/Artigos_ModelosSimétricos/1-Article_LandOceanTemperature/Codes/Land-Ocean Temperature Index (C)_new_1880_2021.txt",h=T)
temperature<-ts(x$No_Smoothing,start=c(1880))
temperature
save(temperature, file="data/temperature.RData")

#write data.R description.

library(roxygen2); # Read in the roxygen2 R package
roxygenise();
