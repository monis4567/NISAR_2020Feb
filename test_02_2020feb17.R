#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#################################################################################
#
# Code developed during NISAR project (2019-2020) by Steen W. Knudsen at NIVA-DK
# R code for analysing filtered and extracted eDNA samples collected by the
# Danish environmental agency (Miløstyrelsen) (MST) over 2017-2018.
# 
# The water samples collected by MST are to be analysed for eDNA levels from 20
# non indeginuous marine species in Danish seas.
#
# The 20 non indeginuous marine species targeted for analysis are the same 20 
# species that the MONIS2, MONIS3 and MONIS4 project, carried out at NIVA-DK
#
# This code is prepared for analysis of the qPCR data obtained in the MONIS5 
# project by laboratorial work carried out in 2018-2019 by SWK
#
# The overall purpose of the present code is to analyse the eDNA levels inferred
# in the qPCR setups performed over 2018-2019.


# This code contains the following sections:
# 1 - Preparation of standard curves with qPCR eDNA levels plotted on to curves. 
#   These plots are prepared as copy numbers versus cycle treshold (Ct) values
# 2 - Plot of sampling locations on maps, for each species,
#   with indication of eDNA intensity for each location monitored
# 3 - Make tables of evaluations of eDNA levels with categories inferred from
#   limit of detection (LOD) and limit of quantification (LOQ)
# 3 - Make tables of evaluations of eDNA replicates analysed in qPCR
#   set ups. With categories assigned from ealier inferred 
#   limit of detection (LOD) and limit of quantification (LOQ).
#
#


#################################################################################

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Use ipdw package to interpolate between marine sampling locations
# interpolate between sampling locations using coastlines as barriers 
#-  as the fish swims, not as the crow flies!
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#code is prepared in 2019-Aug by Steen W. Knudsen
# and is able to run in :
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# > R.Version()
# $platform
# [1] "x86_64-apple-darwin13.4.0"
# $arch
# [1] "x86_64"
# $os
# [1] "darwin13.4.0"
# $system
# [1] "x86_64, darwin13.4.0"
# $status
# [1] ""
# $major
# [1] "3"
# $minor
# [1] "3.3"
# $year
# [1] "2017"
# $month
# [1] "03"
# $day
# [1] "06"
# $`svn rev`
# [1] "72310"
# $language
# [1] "R"
# $version.string
# [1] "R version 3.3.3 (2017-03-06)"
# $nickname
# [1] "Another Canoe"
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#remove everything in the working environment, without a warning!!
rm(list=ls())

#libr.path <- "/home/sknu003/uoa00029_runs/Rplot_tryout" 
#.libPaths( c( libr.path, .libPaths()) ) 

#libr.path <- "/home/sknu003/uoa00029_runs/Rplot_tryout"
#libr.path <- "/scale_wlg_persistent/filesets/home/sknu003/R/x86_64-pc-linux-gnu-library/3.5"

#libr.path <- "/scale_wlg_persistent/filesets/home/sknu003/R/x86_64-pc-linux-gnu-library/3.6"
#.libPaths( c( .libPaths(), libr.path) )

#.libPaths()

#.libPaths( c( libr.path , .libPaths() ) )
#.libPaths()
#.libPaths(libr.path)
#.libPaths()
#chooseCRANmirror(graphics=FALSE)
#chooseCRANmirror(4)
#'chooseCRANmirror(graphics=FALSE, ind=4)'


#________________________________________________________________________________
# Install packages for making standard curve plots and exporting html tables
#________________________________________________________________________________
#see this
#website
#on how to only install required packages
#https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  scales, 
  fields, 
  gplots,
  plyr,
  ReporteRs)



## install the package 'scales', which will allow you to make points on your plot more transparent
#install.packages("scales")
if(!require(scales)){
  install.packages("scales")
  library(scales)
}
library(scales)

#install.packages("fields")
if(!require(fields)){
  install.packages("fields")
  library(fields)
}
library(fields)

## install the package 'gplots', to be able to translate colors to hex - function: col2hex
#install.packages("gplots")
if(!require(gplots)){
  install.packages("gplots")
  library(gplots)
}
library(gplots)

## install the package 'glad', to be able to color using the function 'myPalette'
#install.packages("glad")
#library(glad)

require(graphics)

## install the package 'marmap', which will allow you to plot bathymetric maps
#install.packages("marmap")
#library(marmap)

#get the package that enables the function 'subplot'
#install.packages("TeachingDemos")
#library(TeachingDemos)

#get package to make maps
#install.packages("rworldmap")
#require (rworldmap)

#install.packages("rworldxtra")
#require(rworldxtra)

#get package to read excel files
#install.packages("readxl")
#library(readxl)

#get package to do count number of observations that have the same value at earlier records:
# see this website: https://stackoverflow.com/questions/11957205/how-can-i-derive-a-variable-in-r-showing-the-number-of-observations-that-have-th
#install.packages("plyr")
if(!require(plyr)){
  install.packages("plyr")
  library(plyr)
}
library(plyr)

#get package to make maps - see this website: http://www.molecularecologist.com/2012/09/making-maps-with-r/
#install.packages("mapdata")
#library(mapdata)

#get package to make maps - see this website: http://www.molecularecologist.com/2012/09/making-maps-with-r/
#install.packages("maps")
#library(maps)
# #get package for shapefiles see this website: http://www.molecularecologist.com/2012/09/making-maps-with-r/
# install.packages(maptools)
# library(maptools)  #for shapefiles

# #get package for adding pies on the map
#install.packages("mapplots")
#library(mapplots)

#get the packages required for exporting to a table to word
#install.packages("ReporteRs")
if(!require(ReporteRs)){
  install.packages("ReporteRs")
  library(ReporteRs)
}

devtools::install_github("davidgohel/ReporteRs")
devtools::install_github("davidgohel/officer")

if(!require(officer)){
  install.packages("officer")
  library(officer)
}
library(ReporteRs)
library(officer)


#install.packages("tableHTML")
#https://cran.r-project.org/web/packages/tableHTML/vignettes/tableHTML.html
if(!require(tableHTML)){
  install.packages("tableHTML")
  library(tableHTML)
}
require(tableHTML)





#####################################################################################

# #get package for adding pies and bars on the map
if(!require(mapplots)){
  install.packages("mapplots")
  library(mapplots)
}
#install.packages("mapplots")
library(mapplots)
#get package to make maps - see this website: http://www.molecularecologist.com/2012/09/making-maps-with-r/
if(!require(mapdata)){
  install.packages("mapdata")
  library(mapdata)
}
#install.packages("mapdata")
library(mapdata)
#get package to make maps - see this website: http://www.molecularecologist.com/2012/09/making-maps-with-r/
#install.packages("maps")
if(!require(maps)){
  install.packages("maps")
  library(maps)
}
library(maps)
#get the package that enables the function 'subplot'
#install.packages("TeachingDemos")
if(!require(TeachingDemos)){
  install.packages("TeachingDemos")
  library(TeachingDemos)
}
library(TeachingDemos)
#get package to make maps
#install.packages("rworldmap")
if(!require(rworldmap)){
  install.packages("rworldmap")
  library(rworldmap)
}
require (rworldmap)
#get another package to make maps
#install.packages("rworldxtra")
if(!require(rworldxtra)){
  install.packages("rworldxtra")
  library(rworldxtra)
}
require(rworldxtra)
## install the package 'scales', 
#which will allow you to make points on your plot more transparent
#install.packages("scales")
if(!require(scales)){
  install.packages("scales")
  library(scales)
}
library(scales)
#install.packages("fields")
if(!require(fields)){
  install.packages("fields")
  library(fields)
}
library(fields)
## install the package 'marmap', which will allow you to plot bathymetric maps
#install.packages("marmap")
if(!require(marmap)){
  install.packages("marmap")
  library(marmap)
}
library(marmap)



#________________________________________________________________________________
# Install packages for making maps with eDNA levels mapped and interpolation
# between points
#________________________________________________________________________________

#________________________________________________________________________________
# - use the two spatial dataframes as in this example    https://jsta.github.io/ipdw/articles/ipdw2.html

# also check out this website: https://globalfishingwatch.org/data-blog/working-with-our-downloadable-public-data-in-r/
# and this website: https://www.molecularecologist.com/2015/07/marmap/
#________________________________________________________________________________
# get the rgeos package
if(!require(rgeos)){
  install.packages("rgeos", repos='http://cran.us.r-project.org')
}  
library(rgeos)

# get the ipdw package
if(!require(ipdw)){
  install.packages("ipdw", repos='http://cran.us.r-project.org')
}  
library(ipdw)

# get the scales package
if (!requireNamespace("scales", quietly=TRUE))
  install.packages("scales", repos='http://cran.us.r-project.org')
library(scales)


# get the sf package
if (!requireNamespace("sf", quietly=TRUE))
  install.packages("sf", repos='http://cran.us.r-project.org')
library(sf)

# get the rnaturalearth package
if (!requireNamespace("rnaturalearth", quietly=TRUE))
  install.packages("rnaturalearth", repos='http://cran.us.r-project.org')
library(rnaturalearth)

#Read in the rgdal library
if (!requireNamespace("rgdal", quietly=TRUE))
  install.packages("rgdal", repos='http://cran.us.r-project.org')
library(rgdal)

# get the sp package
if (!requireNamespace("sp", quietly=TRUE))
  install.packages("sp", repos='http://cran.us.r-project.org')
library(sp)

#https://www.rdocumentation.org/packages/biogeo/versions/1.0/topics/dms2dd
# get biogeo package to be able to use 'dms2dd' function
if(!require(biogeo)){
  install.packages("biogeo", repos='http://cran.us.r-project.org')
}
library(biogeo)
#https://www.rdocumentation.org/packages/rgdal/versions/1.3-6/topics/readOGR

#:::::::::::::::::example data set below to try out - start :::::::::::::
# #some values for a dataframe 
# X_UTM <- c(494687.0, 575538.9, 609334.6, 605383.9, 642761.1, 467847.1, 541688.4, 594828.5, 598755.3, 599125.8, 746990.6, 613889.7, 612131.3)
# Y_UTM <- c(6393366, 6452774, 6409207, 6361790, 6336866, 6374410, 6381022, 6354321, 6332203, 6294628, 6085877, 6053287, 6177412)
# zval <- c(1.5, 2.3, 5.6, 8.1, 7.3, 4.6, 3.7, 6.1, 1.8, 2.3, 2.7, 2.8, 7.6)
# #assemble to a dataframe
# loc_DK <- data.frame(X_UTM,Y_UTM,zval)
# #convert from UTM to decimal lat long # https://stackoverflow.com/questions/45230760/convert-utm-to-decimal-degree-in-r
# utmcoor<-SpatialPoints(cbind(loc_DK$X_UTM,loc_DK$Y_UTM), proj4string=CRS("+proj=utm +zone=32"))
# longlatcoor<-spTransform(utmcoor,CRS("+proj=longlat"))
# #extract from spatial points datafram and append back
# loc_DK$declon2 <- coordinates(longlatcoor)[,1]
# loc_DK$declat2 <- coordinates(longlatcoor)[,2]
#:::::::::::::::::example data set below to try out - end :::::::::::::

###########################################################################
# Set working directory and read in csv files

# set working directory to data folder
# setwd("pathToDirHere")
wd ="/home/hal9000/NISAR_analysis"
#wd = "E:/NISAR_analysis/"
#wd <- "/Users/steenknudsen/Documents/Documents/NIVA_Ansaettelse_2019sep/NISAR_analysis"
setwd(wd)
getwd()
#read in the qPCR data with the eDNA reads 
MONIS5eDNA01_df <- 
  read.csv("outfile02_merged_mxpro_txtrepfls_MONIS5.csv", sep = ";",
           stringsAsFactors = FALSE)
MONIS5eDNA02_df <- MONIS5eDNA01_df
#read in the species specific assays used
MONIS3.ls.assays01_df <- 
  read.csv("MONIS3_lst_of_spcfc_eDNA_assays.csv", sep = ",",
           stringsAsFactors = FALSE)
#read in the collected and filtered and extracted water samples
MST_smpls01_df <- 
  read.csv("MST_samples_2017_2018.csv", sep = ",",
           stringsAsFactors = FALSE)
#read in the collected and filtered and extracted water samples
tx_hierc_df <- 
  read.csv("NISAR_spcs_tax_hierarchy01.csv", sep = ",",
           stringsAsFactors = FALSE)


###########################################################################
#keep only selected columns
colnames(MONIS3.ls.assays01_df)
keeps <- c("Assay_ID",
           "Common_name_Danish",
           "Lat_Species",
           "Genus",
           "species")
#keep only selected columns
MONIS3.ls.assays02_df <- MONIS3.ls.assays01_df[keeps]
#get number of columns
n.col.MO3 <- length(MONIS3.ls.assays02_df)
#get only unique columns
MONIS3.ls.assays03_df <- unique(MONIS3.ls.assays02_df[,1:n.col.MO3])
#replace the Danish letters
com.nm.Dan01 <- gsub("\xbf","oe",MONIS3.ls.assays03_df$Common_name_Danish)
com.nm.Dan02 <- gsub("\x8c","aa",com.nm.Dan01)
com.nm.Dan03 <- gsub(" ","_",com.nm.Dan02)
#append back to df
MONIS3.ls.assays03_df$Common_name_Danish <- com.nm.Dan03
# replace "Magallana" with "Crassostrea" # the previous genus name used to be "Crassostrea"
MONIS3.ls.assays03_df$Genus <- gsub("Magallana","Crassostrea",MONIS3.ls.assays03_df$Genus)
MONIS3.ls.assays03_df$Lat_Species <- gsub("Magallana","Crassostrea",MONIS3.ls.assays03_df$Lat_Species)
# replace "serruculata" with "verruculosa" ## the previous species name used to be "serruculata"
MONIS3.ls.assays03_df$species <- gsub("serruculata","verruculosa",MONIS3.ls.assays03_df$species)
MONIS3.ls.assays03_df$Lat_Species <- gsub("serruculata","verruculosa",MONIS3.ls.assays03_df$Lat_Species)

# split text - see: https://stevencarlislewalker.wordpress.com/2013/02/13/remove-or-replace-everything-before-or-after-a-specified-character-in-r-strings/
# and concatenate text - see: https://stackoverflow.com/questions/7201341/how-can-2-strings-be-concatenated 
# to get 6 letter abbr of latin speciesnames
ls.abbr.spcnm <-  paste(
  substr(sub('\\_.*', '', MONIS3.ls.assays03_df$Genus), 1, 3),
  substr(sub('.*\\_', '', MONIS3.ls.assays03_df$species), 1, 3),
  sep="."
)
#remove point with gsub
ls.abbr.spcnm <- gsub("\\.","",ls.abbr.spcnm)
#add back on to latin name dataframe
MONIS3.ls.assays03_df$abbr.nm <- ls.abbr.spcnm
#remove blanks #NOTE!! This will remove all NTC's with "No Ct"
MONIS5eDNA02_df<-na.omit(MONIS5eDNA02_df)
#remove "No Ct"
MONIS5eDNA02_df<-MONIS5eDNA02_df[!grepl("NoCt", MONIS5eDNA02_df$Quantitycopies),]
#change x into numeric variable
MONIS5eDNA02_df$CtdRn=as.numeric(as.character(MONIS5eDNA02_df$CtdRn))
MONIS5eDNA02_df$Quantitycopies=as.numeric(as.character(MONIS5eDNA02_df$Quantitycopies))
#head(MONIS5eDNA02_df,5)
#match between dataframes to add latin species names and DK common names
MONIS5eDNA02_df$gen_specnm <- MONIS3.ls.assays03_df$Lat_Species[match(MONIS5eDNA02_df$speciesabbr, MONIS3.ls.assays03_df$abbr.nm)]
MONIS5eDNA02_df$Genus <- MONIS3.ls.assays03_df$Genus[match(MONIS5eDNA02_df$speciesabbr, MONIS3.ls.assays03_df$abbr.nm)]
MONIS5eDNA02_df$species <- MONIS3.ls.assays03_df$species[match(MONIS5eDNA02_df$speciesabbr, MONIS3.ls.assays03_df$abbr.nm)]
MONIS5eDNA02_df$dk_comnm <- MONIS3.ls.assays03_df$Common_name_Danish[match(MONIS5eDNA02_df$speciesabbr, MONIS3.ls.assays03_df$abbr.nm)]
#check if both lists have the same species - and identify the missing species
sp.A <- unique(MONIS5eDNA02_df$gen_specnm)
sp.B <- MONIS3.ls.assays03_df$Lat_Species
unique(sp.B[! sp.B %in% sp.A])
#match between dataframes
MONIS5eDNA02_df$Lok_omr01 <- MST_smpls01_df$Lok_omr01[match(MONIS5eDNA02_df$smpltp, MST_smpls01_df$U_Pr_Nr)]
MONIS5eDNA02_df$Dato_inds <- MST_smpls01_df$Dato_inds[match(MONIS5eDNA02_df$smpltp, MST_smpls01_df$U_Pr_Nr)]
MONIS5eDNA02_df$Vwf_mL <- MST_smpls01_df$Vwf_mL[match(MONIS5eDNA02_df$smpltp, MST_smpls01_df$U_Pr_Nr)]
MONIS5eDNA02_df$lok_pos_lat <- MST_smpls01_df$lok_pos_lat[match(MONIS5eDNA02_df$smpltp, MST_smpls01_df$U_Pr_Nr)]
MONIS5eDNA02_df$lok_pos_lon <- MST_smpls01_df$lok_pos_lon[match(MONIS5eDNA02_df$smpltp, MST_smpls01_df$U_Pr_Nr)]


#unique(MONIS5eDNA02_df$smpltp)
#paste a new column based on variables separated by point
MONIS5eDNA02_df$Lok_omr01.Welltype <- paste(MONIS5eDNA02_df$Lok_omr01, MONIS5eDNA02_df$WellType,  sep=".")
#get the unique smpl names for locations and WellTypes
#this will allow you to assign a fixed colour to the sampling locations
unHaWT <- unique(MONIS5eDNA02_df$Lok_omr01.Welltype)
# make a transparent color
transp_col <- rgb(0, 0, 0, 0)
#transp_col <- as.character("#FFFFFF")
HaWTnoNA <- addNA(unHaWT)
col.01<-as.numeric(as.factor(unHaWT))
#make a small dataframe w harbours and standards and numbers assigned, 
#use the col2hex in gplot pacakge to convert the 'red' color name to hex-color
col.02 <- col2hex(palette(rainbow(length(col.01))))
lok.cols <- cbind(unHaWT,col.01, col.02)
length(unHaWT)
length(col.01)
length(col.02)
#replace the colour for the standard dilution sample type with the transparent colour
col.03<-replace(col.02, unHaWT=="NA.Standard", transp_col)
col.04 <- cbind(lok.cols,col.03)
#this dataframe now has unique standardized colors for each sampling locality
colfor.lok <- as.data.frame(col.04)
#match to main data frame and add as new color
MONIS5eDNA02_df$col.06 <- colfor.lok$col.03[match(MONIS5eDNA02_df$Lok_omr01.Welltype, colfor.lok$unHaWT)]
#insert the transparent color for all matches with "NA.Standard"
MONIS5eDNA02_df$col.06[MONIS5eDNA02_df$Lok_omr01.Welltype=="NA.Standard"] <- transp_col


######################################################################################
#   Appendix A - start
######################################################################################
####################################################################################
#
# prepare std dilution curve plots for each for species
#
####################################################################################
#first get unique species names 
#get the unique species names
latspecnm <- unique(MONIS5eDNA02_df$gen_specnm)
#match the assay number to the data frame with species
AIfps <- MONIS3.ls.assays03_df$Assay_ID[match(latspecnm, MONIS3.ls.assays03_df$Lat_Species)]
#pad with zeros to two characters
#see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
AIfps <-stringr::str_pad(AIfps, 2, pad = "0")
#make a new data frame with assay Id No and species
nlspnm <- data.frame(AIfps,latspecnm)
#reorder by the column 'AssayIDNo'
nlspnm<- nlspnm[order(nlspnm$AIfps),]
#make a list of numbers for the unique species
no.latspc <- seq(1:length(latspecnm))
#add a new column with no to use for appendix numbering
nlspnm <- cbind(nlspnm, no.latspc) 
#use the new order of latin species names for producing plots
latspecnm <- unique(nlspnm$latspecnm)




######################################################################################
#   make standard curve plots for each species for each season 
######################################################################################
#Use a copy of the data frame to iterate over
amp <- MONIS5eDNA02_df
########################################################
# for loop start here
########################################################
#latspecnm  <- "Mnemiopsis leidyi"
#spec.lat  <- "Mnemiopsis leidyi"
# loop over all species names in the unique list of species, and make plots. 
#Notice that the curly bracket ends after the pdf file is closed
for (spec.lat in latspecnm){
  #print(spec.lat)
  #}
  #get the Danish commom name
  #first split the string by the dot
  #https://stackoverflow.com/questions/33683862/first-entry-from-string-split
  #and escape the dot w two backslashes
  latnm <- sapply(strsplit(spec.lat,"\\."), `[`, 1)
  #get Danish common name
  sbs.dknm <- MONIS3.ls.assays03_df$Common_name_Danish[match(latnm, MONIS3.ls.assays03_df$Lat_Species)]
  #get AssIDNo
  sbs.AssIDNo <- MONIS3.ls.assays03_df$Assay_ID[match(latnm, MONIS3.ls.assays03_df$Lat_Species)]
  #get the number for the appendix plot number
  no.spc.app.plot <- nlspnm$no.latspc[match(spec.lat, nlspnm$latspecnm)]
  #pad with zeros to two characters
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  no.spc.app.plot <-stringr::str_pad(no.spc.app.plot, 2, pad = "0")
  #get the latin species nam without underscore
  spec.lat.no_undersc <- paste(sub('_', ' ', spec.lat))
  # Exporting PFD files via postscript()           
  pdf(c(paste("App_A",no.spc.app.plot,"_plot_qpcr_MONIS5_AssID",sbs.AssIDNo,"_",spec.lat,"_std_dilution_series.pdf",  sep = ""))
      ,width=(1.6*8.2677),height=(2*1.6*2*2.9232))
  #op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
  op <- par(mfrow=c(1,1), # set number of panes inside the plot - i.e. c(2,2) would make four panes for plots
            oma=c(1,1,0,0), # set outer margin (the margin around the combined plot area) - higher numbers increase the number of lines
            mar=c(5,5,5,5) # set the margin around each individual plot 
  )
  #subset based on variable values, subset by species name 
  sbs.amp <- amp[ which(amp$gen_specnm==spec.lat), ]
  #identify LOD
  lod.id.df<-sbs.amp[(sbs.amp$WellType=='Standard'),]
  lod.val<-min(lod.id.df$Quantitycopies)
  #identify LOQ
  #limit the dataframe to only well type that equals standard
  zc<-sbs.amp[(sbs.amp$WellType=='Standard'),]
  #count the occurences of dilution steps - i.e. the number of succesful replicates
  #see this webpage: https://www.miskatonic.org/2012/09/24/counting-and-aggregating-r/
  #zd<-count(zc, "WellName")
  zd<-count(zc, "Quantitycopies")
  #turn this into a dataframe
  ze<-as.data.frame(zd)
  #match the dilution step to the number of occurences -i.e. match between the two dataframes
  no.occ <- ze$freq[match(zc$Quantitycopies,ze$Quantitycopies)]
  #add this column with counted occurences to the limited dataframe
  zg <- cbind.data.frame(zc,no.occ)
  #exlude all observations where less than 3 replicates amplified
  zh<-zg[(zg$no.occ>=3),]
  #get the lowest dilution step that succesfully ampllified on all 3 repliactes
  loq.val=min(zh$Quantitycopies)
  #Conditionally Remove Dataframe Rows with R
  #https://stackoverflow.com/questions/8005154/conditionally-remove-dataframe-rows-with-r
  sbs.pamp<-sbs.amp[!(sbs.amp$WellType=='Standard' & sbs.amp$Quantitycopies<=5),]
  
  #__________________# plot1   - triangles________________________________________
  # Exporting EPS files via postscript()
  # postscript(c(paste("plot_qpcr_MONIS3_",sbs.AssIDNo,"_",spec.lat,"_std_dilution_series.eps", sep = "")),
  #             width=(1.6*8.2677),height=(2*1.6*2.9232),
  #             #family = "Arial", 
  #             paper = "special", onefile = FALSE,
  #             horizontal = FALSE)
  ##  Create a data frame with eDNA
  y.sbs.amp <- sbs.amp$CtdRn
  x.sbs.amp <- sbs.amp$Quantitycopies
  d.sbs.famp <- data.frame( x.sbs.amp = x.sbs.amp, y.sbs.amp = y.sbs.amp )
  #get( getOption( "device" ) )()
  plot(
    y.sbs.amp ~ x.sbs.amp,
    data = d.sbs.famp,
    type = "n",
    log  = "x",
    las=1, # arrange all labels horizontal
    xaxt='n', #surpress tick labels on x-axis
    yaxt='n', #surpress tick labels on y-axis
    #main=c(paste("qPCR standard curve - for ",sbs.AssIDNo,"\n-",spec.lat,seas,"(",sbs.dknm,")"),  sep = ""), 
    
    #add a title with bquote
    main=c(bquote('qPCR standard curve for'~italic(.(spec.lat.no_undersc))
                  ~'('~.(sbs.dknm)~'), '
                  ~'AssayNo'~.(sbs.AssIDNo) #~', '
                  #~.(eng.seas)
    )),
    #offset = 2,
    #sub="sub-title",
    xlab="target-eDNA in extract. (copy/qPCR-reaction)",
    ylab="Ct",
    #xlim = c( 0.1, 1000000000 ),
    #ylim = c( 10, 50 )
    xlim = c( 0.234, 0.428*1000000000 ),
    ylim = c( 9.55, 48.446 )
  )
  #add labels to the points
  pos_vector <- rep(3, length(sbs.amp$Lok_omr01))
  #pos_vector[sbs.amp$Harbour %in% c("Roedby", "Aalborgportland", "KalundborgStatiolHavn")] <- 4
  #pos_vector[sbs.amp$Harbour %in% c("AalborgHavn")] <- 2
  text(x.sbs.amp, y.sbs.amp, labels=sbs.amp$Lok_omr01, cex= 0.8, pos=pos_vector, las=3)
  ##  Put grid lines on the plot, using a light blue color ("lightsteelblue2").
  # add horizontal lines in grid
  abline(
    h   = c( seq( 8, 48, 2 )),
    lty = 1, lwd =0.6,
    col = colors()[ 225 ]
  )
  
  # add vertical lines in grid
  abline(
    v   = c( 
      seq( 0.1, 1, 0.1 ),
      seq( 1e+0, 1e+1, 1e+0 ),
      seq( 1e+1, 1e+2, 1e+1 ),
      seq( 1e+2, 1e+3, 1e+2 ),
      seq( 1e+3, 1e+4, 1e+3 ),
      seq( 1e+4, 1e+5, 1e+4 ), 
      seq( 1e+5, 1e+6, 1e+5 ),
      seq( 1e+6, 1e+7, 1e+6 ),
      seq( 1e+7, 1e+8, 1e+7 ),
      seq( 1e+8, 1e+9, 1e+8 )),
    lty = 1, lwd =0.6,
    col = colors()[ 225 ]
  )
  # add line for LOQ
  abline(v=loq.val, lty=2, lwd=1, col="black")
  text(loq.val*0.7,15,"LOQ",col="black",srt=90,pos=1, font=1)
  
  # add line for LOD 
  abline(v=lod.val, lty=1, lwd=1, col="red")
  text(lod.val*0.7,22,"LOD",col="red",srt=90,pos=1, font=1)
  
  # add line for Ct-cut-off
  abline(h=seq(41,100,1000), lty=1, lwd=3, col="darkgray")
  text(10,40.6,"cut-off",col="darkgray",srt=0,pos=3, font=2, cex=1.2)
  
  # make a transparent color
  #transp_col <- rgb(0, 0, 0, 0)
  #make numbers for the sample type
  #convert NAs to a number 
  # https://stackoverflow.com/questions/27195956/convert-na-into-a-factor-level
  #sbs.amp.stndnm <- addNA(sbs.amp$Harbour)
  #col.01<-as.numeric(as.factor(sbs.amp.stndnm))
  #make a small dataframe w harbours and standards and numbers assigned, 
  #check that the standard is matched up with the transparent color - currently no 17 or 16 ?
  #harbourcols <- cbind(sbs.amp.stndnm,col.01,sbs.amp$Harbour)
  #replace the colour for the standard dilution sample type with the transparent colour
  #col.02<-replace(col.01, col.01==16, transp_col)
  #col.04 <- colforharb$col.02[match(sbs.amp$Harbour.Welltype, colforharb$unHaWT)]
  
  ##  Draw the points over the grid lines.
  points( y.sbs.amp ~ x.sbs.amp, data = d.sbs.famp, 
          pch=c(24), lwd=1, cex=1.8,
          bg=as.character(sbs.amp$col.06)
  )
  #edit labels on the x-axis
  ticks <- seq(-1, 9, by=1)
  labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
  axis(1, at=c(0.1, 1, 10, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7, 1e+8, 1e+9), pos=8, labels=labels)
  #edit labels on the y-axis
  axis(side=2, at=seq(8, 50, by = 2), las=1, pos=0.1)
  
  #estimate a model for each STD subset incl below LOQ
  sbs.amp$x <- sbs.amp$Quantitycopies
  sbs.amp$y<- sbs.amp$CtdRn
  logEst.amp_STD <- lm(y~log(x),sbs.amp)
  
  #add log regresion lines to the plot
  with(as.list(coef(logEst.amp_STD)),
       curve(`(Intercept)`+`log(x)`*log(x),add=TRUE,
             lty=1))
  
  #estimate a model for each STD subset for dilution steps above LOQ
  ab.loq.sbs.amp<-zh # get the previously limited dataframe from identifying LOQ
  ab.loq.sbs.amp$x <- ab.loq.sbs.amp$Quantitycopies
  ab.loq.sbs.amp$y<- ab.loq.sbs.amp$CtdRn
  logEst.abloqamp_STD <- lm(y~log(x),ab.loq.sbs.amp) #make a linear model
  
  #add log regresion lines to the plot
  with(as.list(coef(logEst.abloqamp_STD)),
       curve(`(Intercept)`+`log(x)`*log(x),add=TRUE,
             lty=1, col="red"))
  
  #add 95% confidence intervals around each fitted line
  #inspired from this webpage
  #https://stat.ethz.ch/pipermail/r-help/2007-November/146285.html
  #for the first line - with below LOQ
  newx<-seq(lod.val,1e+6,1000)
  prdlogEst.amp_STD<-predict(logEst.amp_STD,newdata=data.frame(x=newx),interval = c("confidence"), 
                             level = 0.95, type="response")
  prd2logEst.amp_STD<- prdlogEst.amp_STD
  #polygon(c(rev(newx), newx), c(rev(prd2[ ,3]), prd2[ ,2]), col = 'grey80', border = NA)
  lines(newx,prd2logEst.amp_STD[,2],col="black",lty=2)
  lines(newx,prd2logEst.amp_STD[,3],col="black",lty=2)
  #add 95% conf. intervals for the second line - only above LOQ
  newx<-seq(loq.val,1e+6,100)
  prdlogEst.abloqamp_STD<-predict(logEst.abloqamp_STD,newdata=data.frame(x=newx),interval = c("confidence"), 
                                  level = 0.95 , type="response")
  prd2logEst.abloqamp_STD<- prdlogEst.abloqamp_STD
  #polygon(c(rev(newx), newx), c(rev(prd2[ ,3]), prd2[ ,2]), col = 'grey80', border = NA)
  lines(newx,prd2logEst.abloqamp_STD[,2],col="red",lty=2)
  lines(newx,prd2logEst.abloqamp_STD[,3],col="red",lty=2)
  # add a legend for colors on points
  legend(1e+7*0.5,49,
         unique(sbs.amp$Lok_omr01.Welltype),
         pch=c(24),
         bg="white",
         #NOTE!! the hex color numbers must be read as characters to translate into hex colors
         pt.bg = as.character(unique(sbs.amp$col.06)),
         y.intersp= 0.7, cex=0.9)
  # add a second legend for types of regression lines
  legend(1000,49,
         c("incl below LOQ","excl below LOQ"),
         #pch=c(24), #uncomment to get triangles on the line in the legend
         cex=0.8,
         bg="white",
         lty=c(1), col=c("black","red"),
         y.intersp= 0.7)
  #title(main=c(paste("qPCR standard curve - for ",spec.lat,"\n-(",sbs.dknm,")"),  sep = ""), 
  #        col.main="red",
  #    sub="My Sub-title", col.sub="blue",
  #    xlab="My X label", ylab="My Y label",
  #    col.lab="green", cex.lab=0.75)
  
  # add title for the pdf-page
  mtext(c(paste("Appendix A",no.spc.app.plot,"."),  sep = ""), outer=TRUE, 
        #use at , adj and padj to adjust the positioning
        at=par("usr")[1]+0.15*diff(par("usr")[1:2]),
        adj=3.4,
        padj=2,
        #use side to place it in te top
        side=3, cex=1.6, line=-1.15)
  
  #apply the par settings for the plot as defined above.
  par(op)
  # end pdf file to save as
  dev.off()  
  ########################################################
  # for loop on species end here
  ########################################################
  
}
######################################################################################
######################################################################################
#   Appendix A - end
######################################################################################












#################################################################################

#################################################################################
#################################################################################
# #get package for adding pies and bars on the map
if(!require(mapplots)){
  install.packages("mapplots")
  library(mapplots)
}
#install.packages("mapplots")
library(mapplots)
#get package to make maps - see this website: http://www.molecularecologist.com/2012/09/making-maps-with-r/
if(!require(mapdata)){
  install.packages("mapdata")
  library(mapdata)
}
#install.packages("mapdata")
library(mapdata)
#get package to make maps - see this website: http://www.molecularecologist.com/2012/09/making-maps-with-r/
#install.packages("maps")
if(!require(maps)){
  install.packages("maps")
  library(maps)
}
library(maps)
#get the package that enables the function 'subplot'
#install.packages("TeachingDemos")
if(!require(TeachingDemos)){
  install.packages("TeachingDemos")
  library(TeachingDemos)
}
library(TeachingDemos)
#get package to make maps
#install.packages("rworldmap")
if(!require(rworldmap)){
  install.packages("rworldmap")
  library(rworldmap)
}
require (rworldmap)
#get another package to make maps
#install.packages("rworldxtra")
if(!require(rworldxtra)){
  install.packages("rworldxtra")
  library(rworldxtra)
}
require(rworldxtra)
## install the package 'scales', 
#which will allow you to make points on your plot more transparent
#install.packages("scales")
if(!require(scales)){
  install.packages("scales")
  library(scales)
}
library(scales)
#install.packages("fields")
if(!require(fields)){
  install.packages("fields")
  library(fields)
}
library(fields)
## install the package 'marmap', which will allow you to plot bathymetric maps
#install.packages("marmap")
if(!require(marmap)){
  install.packages("marmap")
  library(marmap)
}
library(marmap)



#################################################################################
#  prepare data frames for  concentration in copies per Liter of seawater
#################################################################################
#copy the data frame
MONIS5eDNA03_df <- MONIS5eDNA02_df
#set NA blanks to zero
MONIS5eDNA03_df$CtdRn[is.na(MONIS5eDNA03_df$CtdRn)] <- 0
MONIS5eDNA03_df$Quantitycopies[is.na(MONIS5eDNA03_df$Quantitycopies)] <- 0
#make sure numbers are numbers
MONIS5eDNA03_df$Quantitycopies <- as.numeric(as.character(MONIS5eDNA03_df$Quantitycopies))
MONIS5eDNA03_df$Vwf_mL <- as.numeric(as.character(MONIS5eDNA03_df$Vwf_mL))
#add column with copies per Liter of filtered water
#Ae = (Cqpcr /Fe) /Vwf. 
#’Ae’ number of  eDNA-copies per volumen filtered water, 
#’Cqpcr’ number of copies detected in the qPCR-well, #smpls02.1$meanQuantitycopies 
#’Fe’ the ratio of the eluted extrated filtrate used in a qPCR-well #5/350
# For the qPCR 5 uL of extracted eDNA from the filter sample was used
# For the extraction of eDNA from the filters the extracted eDNA was eluated in 350 uL
#’Vwf’ is volumen of seawater filtered.
rt <- 5/350
# get the template volume used for qPCR setups
MONIS5eDNA03_df$templvol2 <- as.numeric(substr(MONIS5eDNA03_df$templvol, 1, 1))
# relate the volume of template used in qPCR to recalculate the ratio of
#extracted eluate used
rt <- (MONIS5eDNA03_df$templvol2/350)
#per mL
MONIS5eDNA03_df$copies_per_mLwater <- (MONIS5eDNA03_df$Quantitycopies/(rt))/MONIS5eDNA03_df$Vwf_mL
#per Liter
MONIS5eDNA03_df$copies_per_Lwater <- MONIS5eDNA03_df$copies_per_mLwater*1000
#replace nas with zeros
MONIS5eDNA03_df$copies_per_Lwater[is.na(MONIS5eDNA03_df$copies_per_Lwater)]<-0
#add one to be able to do logarithmic scales
MONIS5eDNA03_df$copies_per_Lwater_plone<- MONIS5eDNA03_df$copies_per_Lwater+1
#take log10 to all copies
MONIS5eDNA03_df$log.10_copies_L <- log10(MONIS5eDNA03_df$copies_per_Lwater_plone)
# following this example: https://stackoverflow.com/questions/41336606/accessing-element-of-a-split-string-in-r
#split the collection date, and get the third element, to get the year
MONIS5eDNA03_df$year_inds <- MONIS5eDNA03_df$Dato_inds %>%
  strsplit( "/" ) %>%
  sapply( "[", 3 )
#check what years the samples have been collected
unique(MONIS5eDNA03_df$year_inds)
#paste latin species names and year for collection together
MONIS5eDNA03_df$gen_specnm.year <- paste(MONIS5eDNA03_df$gen_specnm,MONIS5eDNA03_df$year_inds, sep=".")
#get the LOD (Limit of detection) for each set of qPCR run
#use the function aggregate to get the minimum value for a group
lodtable1 <- aggregate(MONIS5eDNA03_df[, "Quantitycopies"], list(MONIS5eDNA03_df$gen_specnm.year, MONIS5eDNA03_df$WellType), min)
#subset this table by group
lodtable2 <- lodtable1[ which(lodtable1$Group.2=="Standard"), ]
#rename the column names
colnames(lodtable2) <- c("spc.year","WellT","LOD")
#identify LOQ for each qPCR run
#limit the dataframe to only well type that equals standard
oc<-MONIS5eDNA03_df[(MONIS5eDNA03_df$WellType=='Standard'),]
#add a new column that merges two columns
oc$Quan.spc.years <- paste(oc$Quantitycopies, oc$gen_specnm.year,  sep=".")
#count the occurences of dilution steps - i.e. the number of succesful replicates
#see this webpage: https://www.miskatonic.org/2012/09/24/counting-and-aggregating-r/
#and this webpage: https://stackoverflow.com/questions/9809166/count-number-of-rows-within-each-group
od<-count(oc, c("Quantitycopies","gen_specnm.year"))
#turn this into a dataframe
oe<-as.data.frame(od)
#add a new column that merges two columns
oe$Quan.spc.years <- paste(as.character(oe$Quantitycopies), oe$gen_specnm.year,  sep=".")
#match the dilution step to the number of occurences -i.e. match between the two dataframes
no.occ <- oe$freq[match(oc$Quan.spc.years,oe$Quan.spc.years)]
#add this column with counted occurences to the limited dataframe
og <- cbind.data.frame(oc,no.occ)
#exlude all observations where less than 3*3 replicates amplified
oh<-og[(og$no.occ>=3),]
#get the lowest dilution step that succesfully ampllified on all 3*3 repliactes
# there were 3 qPCR runs per species, and each run included 
#3 replicates of each standard dilution step
#use aggregate to get the minimum for each
loqtable1 <- aggregate(oh[, "Quantitycopies"], list(oh$gen_specnm.year), min)
#change the column names
colnames(loqtable1) <- c("spc.year","LOQ")
#copy the LOD table and add the corresponding LOQ values
loq.lod.table <- lodtable2
#match back to data frame
loq.lod.table$LOQ <- loqtable1$LOQ[match(lodtable2$spc.year,loqtable1$spc.year)]

#make a copy of the original data frame
# The previous copy of the original data frame had all wells with zero
# qPCR reads removed -  you will need to include the zero detections
MONIS5eDNA04_df <- MONIS5eDNA01_df
#Replace all NoCts with zero
MONIS5eDNA04_df$Quantitycopies[MONIS5eDNA04_df$Quantitycopies=="NoCt"] <- 0
#delete rows that have "NoCtforFAMStandards" in the column 'smpls03$Quantitycopies'
#this error is caused by the wrong qpcr plate setup for  'Pseudochattonella_farcimen'
MONIS5eDNA04_df<-MONIS5eDNA04_df[!(MONIS5eDNA04_df$Quantitycopies=="NoCtforFAMStandards"),]
#change variable into numeric variable
MONIS5eDNA04_df$CtdRn=as.numeric(as.character(MONIS5eDNA04_df$CtdRn))
# replace NAs with zeros
MONIS5eDNA04_df$CtdRn[is.na(MONIS5eDNA04_df$CtdRn)] <- 0
MONIS5eDNA04_df$Quantitycopies=as.numeric(as.character(MONIS5eDNA04_df$Quantitycopies))
#match between dataframes to add latin species names and DK common names
MONIS5eDNA04_df$gen_specnm <- MONIS3.ls.assays03_df$Lat_Species[match(MONIS5eDNA04_df$speciesabbr, MONIS3.ls.assays03_df$abbr.nm)]
MONIS5eDNA04_df$Genus <- MONIS3.ls.assays03_df$Genus[match(MONIS5eDNA04_df$speciesabbr, MONIS3.ls.assays03_df$abbr.nm)]
MONIS5eDNA04_df$species <- MONIS3.ls.assays03_df$species[match(MONIS5eDNA04_df$speciesabbr, MONIS3.ls.assays03_df$abbr.nm)]
MONIS5eDNA04_df$dk_comnm <- MONIS3.ls.assays03_df$Common_name_Danish[match(MONIS5eDNA04_df$speciesabbr, MONIS3.ls.assays03_df$abbr.nm)]
#match between data frames
MONIS5eDNA04_df$Lok_omr01 <-    MST_smpls01_df$Lok_omr01[match(MONIS5eDNA04_df$smpltp, MST_smpls01_df$U_Pr_Nr)]
MONIS5eDNA04_df$Dato_inds <-    MST_smpls01_df$Dato_inds[match(MONIS5eDNA04_df$smpltp, MST_smpls01_df$U_Pr_Nr)]
MONIS5eDNA04_df$Vwf_mL <-       MST_smpls01_df$Vwf_mL[match(MONIS5eDNA04_df$smpltp, MST_smpls01_df$U_Pr_Nr)]
MONIS5eDNA04_df$lok_pos_lat <-  MST_smpls01_df$lok_pos_lat[match(MONIS5eDNA04_df$smpltp, MST_smpls01_df$U_Pr_Nr)]
MONIS5eDNA04_df$lok_pos_lon <-  MST_smpls01_df$lok_pos_lon[match(MONIS5eDNA04_df$smpltp, MST_smpls01_df$U_Pr_Nr)]
# make variables numeric
MONIS5eDNA04_df$Quantitycopies <- as.numeric(as.character(MONIS5eDNA04_df$Quantitycopies))
MONIS5eDNA04_df$Vwf_mL <- as.numeric(as.character(MONIS5eDNA04_df$Vwf_mL))
# following this example: https://stackoverflow.com/questions/41336606/accessing-element-of-a-split-string-in-r
#split the collection date, and get the third element, to get the year
MONIS5eDNA04_df$year_inds <- MONIS5eDNA04_df$Dato_inds %>%
  strsplit( "/" ) %>%
  sapply( "[", 3 )
#paste a new column based on variables separated by point
MONIS5eDNA04_df$WellType.year_inds <- paste(MONIS5eDNA04_df$WellType, MONIS5eDNA04_df$year_inds,  sep=".")
#paste a new column based on variables separated by point
MONIS5eDNA04_df$gen_specnm.year_inds <- paste(MONIS5eDNA04_df$gen_specnm, MONIS5eDNA04_df$year_inds,  sep=".")
#paste a new column based on variables separated by point
MONIS5eDNA04_df$Lok_omr01.Welltype <- paste(MONIS5eDNA04_df$Lok_omr01, MONIS5eDNA04_df$WellType,  sep=".")
#get the unique smpl names for Harbours and WellTypes
unHaWT <- unique(MONIS5eDNA04_df$Lok_omr01.Welltype)
#transp_col <- as.character("#FFFFFF")
HaWTnoNA <- addNA(unHaWT)
col.01<-as.numeric(as.factor(HaWTnoNA))
#make a small dataframe w harbours and standards and numbers assigned, 
#use the col2hex in gplot pacakge to convert the 'red' color name to hex-color
col.02 <- col2hex(palette(rainbow(length(col.01)+1)))

lokomrcols <- cbind(unHaWT,col.01, col.02)
length(unHaWT)
length(col.01) 
length(col.02)
#replace the colour for the standard dilution sample type with the transparent colour
col.03<-replace(col.02, unHaWT=="NA.Standard", transp_col)
col.04 <- cbind(lokomrcols,col.03)
colforlokomr <- as.data.frame(col.04)
#match to main data frame and add as new color
MONIS5eDNA04_df$col.06 <- colforlokomr$col.03[match(MONIS5eDNA04_df$Lok_omr01.Welltype, colforlokomr$unHaWT)]
#insert the transparent color for all matches with "NA.Standard"
MONIS5eDNA04_df$col.06[MONIS5eDNA04_df$Lok_omr01.Welltype=="NA.Standard"] <- transp_col
# following this example: https://stackoverflow.com/questions/41336606/accessing-element-of-a-split-string-in-r
#split the collection date, and get the third element, to get the year
loq.lod.table$gen_spcnm <- loq.lod.table$spc.year %>%
  strsplit( "\\." ) %>%
  sapply( "[", 1 )
#match between the smpls03 dataframe and the LOD and LOQ table
MONIS5eDNA04_df$LOD <- loq.lod.table$LOD[match(MONIS5eDNA04_df$gen_specnm, loq.lod.table$gen_spcnm)]
MONIS5eDNA04_df$LOQ <- loq.lod.table$LOQ[match(MONIS5eDNA04_df$gen_specnm, loq.lod.table$gen_spcnm)]

#add column with copies per Liter of filtered water
#Ae = (Cqpcr /Fe) /Vwf. 
#’Ae’ number of  eDNA-copies per volumen filtered water, 
#’Cqpcr’ number of copies detected in the qPCR-well, #smpls02.1$meanQuantitycopies 
#’Fe’ the ratio of the eluted extrated filtrate used in a qPCR-well #5/350
# For the qPCR 5 uL of extracted eDNA from the filter sample was used
# For the extraction of eDNA from the filters the extracted eDNA was eluated in 350 uL
#’Vwf’ is volumen of seawater filtered.
rt <- 5/350
# get the template volume used for qPCR setups
MONIS5eDNA04_df$templvol2 <- as.numeric(substr(MONIS5eDNA04_df$templvol, 1, 1))
# relate the volume of template used in qPCR to recalculate the ratio of
#extracted eluate used
rt <- (MONIS5eDNA04_df$templvol2/350)
#per mL
MONIS5eDNA04_df$copies_per_mLwater <- (MONIS5eDNA04_df$Quantitycopies/(rt))/MONIS5eDNA04_df$Vwf_mL
#per Liter
MONIS5eDNA04_df$copies_per_Lwater <- MONIS5eDNA04_df$copies_per_mLwater*1000
#replace nas with zeros
MONIS5eDNA04_df$copies_per_Lwater[is.na(MONIS5eDNA04_df$copies_per_Lwater)]<-0
#add one to be able to do logarithmic scales
MONIS5eDNA04_df$copies_per_Lwater_plone<- MONIS5eDNA04_df$copies_per_Lwater+1
#take log10 to all copies
MONIS5eDNA04_df$log.10_copies_L <- log10(MONIS5eDNA04_df$copies_per_Lwater_plone)
#make a subset 
MONIS5eDNA05_df <- MONIS5eDNA04_df[ which(MONIS5eDNA04_df$WellType=="Unknown"), ]
# exclude the rows that have NAs for years
MONIS5eDNA06_df <- MONIS5eDNA05_df[which(!is.na(MONIS5eDNA05_df$year_inds)), ] 
#head(MONIS5eDNA06_df,5)
######################################################################################################
# Get mean for each set of 3 technical qPCR replicates per species per season per harbour
######################################################################################################
#get the mean quantity for each species per season per port
MONIS5eDNA07_df <- aggregate(MONIS5eDNA06_df[, "Quantitycopies"], list(MONIS5eDNA06_df$gen_specnm.year_inds, MONIS5eDNA06_df$smpltp), mean)
#change the column names
colnames(MONIS5eDNA07_df) <- c("spc.year","MSTsmpl","meanQuantitycopies")
# following this example: https://stackoverflow.com/questions/41336606/accessing-element-of-a-split-string-in-r
#split the collection date, and get the third element, to get the year
MONIS5eDNA07_df$gen_spcnm <- MONIS5eDNA07_df$spc.year %>%
  strsplit( "\\." ) %>%
  sapply( "[", 1 )
#match with LOD and LOQ
MONIS5eDNA07_df$LOD <- loq.lod.table$LOD[match(MONIS5eDNA07_df$gen_spcnm, loq.lod.table$gen_spcnm)]
MONIS5eDNA07_df$LOQ <- loq.lod.table$LOQ[match(MONIS5eDNA07_df$gen_spcnm, loq.lod.table$gen_spcnm)]
#add an empty column with just NAs
MONIS5eDNA07_df[,"eDNA_eval_mean"] <- NA
#replace in the empty column, the order is important, as you otherwise will end up with the last evaluations
MONIS5eDNA07_df$eDNA_eval_mean[MONIS5eDNA07_df$meanQuantitycopies>MONIS5eDNA07_df$LOQ] <- "aboveLOQ"
MONIS5eDNA07_df$eDNA_eval_mean[MONIS5eDNA07_df$meanQuantitycopies<MONIS5eDNA07_df$LOQ] <- "AbLOD_BeLOQ"
MONIS5eDNA07_df$eDNA_eval_mean[MONIS5eDNA07_df$meanQuantitycopies<MONIS5eDNA07_df$LOD & !MONIS5eDNA07_df$meanQuantitycopies==0] <- "belowLOD"
MONIS5eDNA07_df$eDNA_eval_mean[MONIS5eDNA07_df$meanQuantitycopies==0] <- "NoCt"
#Match WellType back to dataframe
MONIS5eDNA07_df$WellType <- MONIS5eDNA06_df$WellType[match(MONIS5eDNA07_df$MSTsmpl, MONIS5eDNA06_df$smpltp)]
#add an empty column with just NAs to fil with color codings
MONIS5eDNA07_df[,"eDNA_col_eval_mean"] <- NA
#replace in the empty column, the order is important, as you otherwise will end up with the last evaluations
MONIS5eDNA07_df$eDNA_col_eval_mean[MONIS5eDNA07_df$meanQuantitycopies>MONIS5eDNA07_df$LOQ] <- "black" # "green"
MONIS5eDNA07_df$eDNA_col_eval_mean[MONIS5eDNA07_df$meanQuantitycopies<MONIS5eDNA07_df$LOQ] <- "red" #"orange" 
MONIS5eDNA07_df$eDNA_col_eval_mean[MONIS5eDNA07_df$meanQuantitycopies<MONIS5eDNA07_df$LOD & !MONIS5eDNA07_df$meanQuantitycopies==0] <-"yellow" #"yellow"
MONIS5eDNA07_df$eDNA_col_eval_mean[MONIS5eDNA07_df$meanQuantitycopies==0] <- "white" #,"red"
#head(MONIS5eDNA07_df,6)
#reorder the dataframe, to make it easier to look at
MONIS5eDNA08_df <- MONIS5eDNA07_df[ order(MONIS5eDNA07_df$spc.year, MONIS5eDNA07_df$MSTsmpl), ]
# following this example: https://stackoverflow.com/questions/41336606/accessing-element-of-a-split-string-in-r
#split the collection date, and get the third element, to get the year
MONIS5eDNA08_df$year.inds <- MONIS5eDNA08_df$spc.year %>%
  strsplit( "\\." ) %>%
  sapply( "[", 2 )
#head(MONIS5eDNA08_df,4)
#match between data frames
MONIS5eDNA08_df$Lok_omr01 <-    MST_smpls01_df$Lok_omr01[match(MONIS5eDNA08_df$MSTsmpl, MST_smpls01_df$U_Pr_Nr)]
MONIS5eDNA08_df$Dato_inds <-    MST_smpls01_df$Dato_inds[match(MONIS5eDNA08_df$MSTsmpl, MST_smpls01_df$U_Pr_Nr)]
MONIS5eDNA08_df$Vwf_mL <-       MST_smpls01_df$Vwf_mL[match(MONIS5eDNA08_df$MSTsmpl, MST_smpls01_df$U_Pr_Nr)]
MONIS5eDNA08_df$lok_pos_lat <-  MST_smpls01_df$lok_pos_lat[match(MONIS5eDNA08_df$MSTsmpl, MST_smpls01_df$U_Pr_Nr)]
MONIS5eDNA08_df$lok_pos_lon <-  MST_smpls01_df$lok_pos_lon[match(MONIS5eDNA08_df$MSTsmpl, MST_smpls01_df$U_Pr_Nr)]
#get unique rows for more than one variable
tmp.vol <- unique(MONIS5eDNA06_df[c("gen_specnm", "templvol2")])
colnames(MONIS5eDNA06_df)
#get unique rows for more than one variable
smpltp.Vwf_mL <- unique(MONIS5eDNA06_df[c("smpltp", "Vwf_mL")])
# add template volume column to data frame
MONIS5eDNA08_df$templvol2 <-  tmp.vol$templvol2[match(MONIS5eDNA08_df$gen_spcnm, tmp.vol$gen_specnm)]
MONIS5eDNA08_df$Vwf_mL <-  smpltp.Vwf_mL$Vwf_mL[match(MONIS5eDNA08_df$MSTsmpl, smpltp.Vwf_mL$smpltp)]
#Ae = (Cqpcr /Fe) /Vwf. 
#’Ae’ number of  eDNA-copies per volumen filtered water, 
#’Cqpcr’ number of copies detected in the qPCR-well, #smpls20$meanQuantitycopies 
#’Fe’ the ratio of the eluted extrated filtrate used in a qPCR-well #5/350
#’Vwf’ is volumen of seawater filtered. #smpls20$volfilt_mL
# get the template volume used for qPCR setups
# relate the volume of template used in qPCR to recalculate the ratio of
#extracted eluate used
rt <- (MONIS5eDNA08_df$templvol2/350)
#per mL
MONIS5eDNA08_df$copies_per_mLwater <- (MONIS5eDNA08_df$meanQuantitycopies/(5/350))/MONIS5eDNA08_df$Vwf_mL
#per Liter
MONIS5eDNA08_df$copies_per_Lwater <- MONIS5eDNA08_df$copies_per_mLwater*1000
#replace nas with zeros
MONIS5eDNA08_df$copies_per_Lwater[is.na(MONIS5eDNA08_df$copies_per_Lwater)]<-0
#add one to be able to do logarithmic scales
MONIS5eDNA08_df$copies_per_Lwater_plone<- MONIS5eDNA08_df$copies_per_Lwater+1
#take log10 to all copies
MONIS5eDNA08_df$log.10_copies_L <- log10(MONIS5eDNA08_df$copies_per_Lwater_plone)

#add an empty column with just NAs to fil with evaluations
MONIS5eDNA08_df[,"no_for_log.10_eDNAlvls"] <- NA

#replace in the empty column, the order is important, 
#as you otherwise will end up with the last evaluations
MONIS5eDNA08_df$no_for_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L==0 ] <- 1 #if No Ct

MONIS5eDNA08_df$no_for_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L<=
                                         log10((MONIS5eDNA08_df$LOD/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) 
                                       & MONIS5eDNA08_df$log.10_copies_L>0] <- 2 #if below LOD but above zero

MONIS5eDNA08_df$no_for_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                         log10((MONIS5eDNA08_df$LOD/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) ] <- 3 #if above LOD

MONIS5eDNA08_df$no_for_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                         log10((MONIS5eDNA08_df$LOQ/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) ] <- 4 #if above LOQ

MONIS5eDNA08_df$no_for_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                         log10((MONIS5eDNA08_df$LOQ/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) &
                                         MONIS5eDNA08_df$copies_per_Lwater>=1] <- 5 #if above LOQ, and within 1-10 copies per L

MONIS5eDNA08_df$no_for_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                         log10((MONIS5eDNA08_df$LOQ/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) &
                                         MONIS5eDNA08_df$copies_per_Lwater>=10] <- 6 #if above LOQ, and within 10-100 copies per L

MONIS5eDNA08_df$no_for_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                         log10((MONIS5eDNA08_df$LOQ/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) &
                                         MONIS5eDNA08_df$copies_per_Lwater>=100] <- 7 #if above LOQ, and within 100-1000 copies per L

MONIS5eDNA08_df$no_for_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                         log10((MONIS5eDNA08_df$LOQ/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) &
                                         MONIS5eDNA08_df$copies_per_Lwater>=1000] <- 8 #if above LOQ, and within 1e3-1e4 copies per L

MONIS5eDNA08_df$no_for_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                         log10((MONIS5eDNA08_df$LOQ/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) &
                                         MONIS5eDNA08_df$copies_per_Lwater>=1E4] <- 9 #if above LOQ, and within 1e4-1e5 copies per L

MONIS5eDNA08_df$no_for_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                         log10((MONIS5eDNA08_df$LOQ/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) &
                                         MONIS5eDNA08_df$copies_per_Lwater>=1E5] <- 10 #if above LOQ, and within 1e5-1e6 copies per L

MONIS5eDNA08_df$no_for_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                         log10((MONIS5eDNA08_df$LOQ/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) &
                                         MONIS5eDNA08_df$copies_per_Lwater>=1E6] <- 11 #if above LOQ, and within 1e6-1e7 copies per L

MONIS5eDNA08_df$no_for_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                         log10((MONIS5eDNA08_df$LOQ/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) &
                                         MONIS5eDNA08_df$copies_per_Lwater>=1E7] <- 12 #if above LOQ, and within 1e7-1e8 copies per L

unique(MONIS5eDNA08_df$no_for_log.10_eDNAlvls)

#add an empty column with just NAs to fil with evaluations
MONIS5eDNA08_df[,"col_log.10_eDNAlvls"] <- NA

#replace in the empty column, the order is important, 
#as you otherwise will end up with the last evaluations
MONIS5eDNA08_df$col_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L==0 ] <- 0 #if No Ct

MONIS5eDNA08_df$col_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L<=
                                      log10((MONIS5eDNA08_df$LOD/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) 
                                    & MONIS5eDNA08_df$log.10_copies_L>0] <- 0 #if below LOD but above zero

MONIS5eDNA08_df$col_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                      log10((MONIS5eDNA08_df$LOD/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) ] <- 0 #if above LOD

MONIS5eDNA08_df$col_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                      log10((MONIS5eDNA08_df$LOQ/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) ] <- 1 #if above LOQ

MONIS5eDNA08_df$col_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                      log10((MONIS5eDNA08_df$LOQ/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) &
                                      MONIS5eDNA08_df$copies_per_Lwater>=1] <- 2 #if above LOQ, and within 1-10 copies per

MONIS5eDNA08_df$col_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                      log10((MONIS5eDNA08_df$LOQ/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) &
                                      MONIS5eDNA08_df$copies_per_Lwater>=10] <- 3 #if above LOQ, and within 10-100 copies per L

MONIS5eDNA08_df$col_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                      log10((MONIS5eDNA08_df$LOQ/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) &
                                      MONIS5eDNA08_df$copies_per_Lwater>=100] <- 4 #if above LOQ, and within 100-1000 copies per L

MONIS5eDNA08_df$col_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                      log10((MONIS5eDNA08_df$LOQ/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) &
                                      MONIS5eDNA08_df$copies_per_Lwater>=1000] <- 5 #if above LOQ, and within 1e3-1e4 copies per L

MONIS5eDNA08_df$col_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                      log10((MONIS5eDNA08_df$LOQ/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) &
                                      MONIS5eDNA08_df$copies_per_Lwater>=1E4] <- 6 #if above LOQ, and within 1e4-1e5 copies per L

MONIS5eDNA08_df$col_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                      log10((MONIS5eDNA08_df$LOQ/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) &
                                      MONIS5eDNA08_df$copies_per_Lwater>=1E5] <- 7 #if above LOQ, and within 1e5-1e6 copies per L

MONIS5eDNA08_df$col_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                      log10((MONIS5eDNA08_df$LOQ/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) &
                                      MONIS5eDNA08_df$copies_per_Lwater>=1E6] <- 8 #if above LOQ, and within 1e6-1e7 copies per L

MONIS5eDNA08_df$col_log.10_eDNAlvls[MONIS5eDNA08_df$log.10_copies_L>=
                                      log10((MONIS5eDNA08_df$LOQ/(rt))/MONIS5eDNA08_df$Vwf_mL*1000) &
                                      MONIS5eDNA08_df$copies_per_Lwater>=1E7] <- 9 #if above LOQ, and within 1e7-1e8 copies per L
#check unique values
maxcolval <- max(unique(MONIS5eDNA08_df$col_log.10_eDNAlvls))






######################################################################################################
# prepare a data frame with eDNA evalution for higest detection lvl for each of 
# 3 techincal qPCR replicates per species per season per harbour
######################################################################################################

#copy the data frame
MONIS5eDNA09_df <- MONIS5eDNA06_df
#make anew column with evaluations in
MONIS5eDNA09_df[,"eDNA_eval_mean"] <- NA
#replace in the empty column, the order is important, as you otherwise will end up with the last evaluations
MONIS5eDNA09_df$eDNA_eval_mean[MONIS5eDNA09_df$Quantitycopies>MONIS5eDNA09_df$LOQ] <- "aboveLOQ"
MONIS5eDNA09_df$eDNA_eval_mean[MONIS5eDNA09_df$Quantitycopies<MONIS5eDNA09_df$LOQ] <- "AbLOD_BeLOQ"
MONIS5eDNA09_df$eDNA_eval_mean[MONIS5eDNA09_df$Quantitycopies<MONIS5eDNA09_df$LOD & !MONIS5eDNA09_df$Quantitycopies==0] <- "belowLOD"
MONIS5eDNA09_df$eDNA_eval_mean[MONIS5eDNA09_df$Quantitycopies==0] <- "NoCt"


colnames(MONIS5eDNA09_df)
head(MONIS5eDNA09_df,4)
#count the number of evaluations per species per harbour per season for each set of 3 replicates
pc<-MONIS5eDNA09_df #copy the data frame
#add a new column with unified values from columns
pc$gen_specnm.year_inds.eDNA_eval_mean.smpltp <- paste(pc$gen_specnm.year_inds, pc$eDNA_eval_mean, pc$smpltp, sep=".")
#count the occurences of dilution steps - i.e. the number of succesful replicates
#see this webpage: https://www.miskatonic.org/2012/09/24/counting-and-aggregating-r/
#and this webpage: https://stackoverflow.com/questions/9809166/count-number-of-rows-within-each-group
pd<-count(pc, c("gen_specnm.year_inds.eDNA_eval_mean.smpltp"))
#turn this into a dataframe
pe<-as.data.frame(pd)
#match the dilution step to the number of occurences -i.e. match between the two dataframes
no.occ <- pe$freq[match(pc$gen_specnm.year_inds.eDNA_eval_mean.smpltp,pe$gen_specnm.year_inds.eDNA_eval_mean.smpltp)]
#add this column with counted occurences to the limited dataframe
pg <- cbind.data.frame(pc,no.occ)
#add another column
pg$spc.year.smpltp <- paste(pg$gen_specnm.year_inds, pg$smpltp, sep=".")

unique(pg$no.occ)
#make subsets
pg_aboveLOQ <- pg[ which(pg$eDNA_eval_mean=="aboveLOQ"), ]
pg_AbLOD_BeLOQ <- pg[ which(pg$eDNA_eval_mean=="AbLOD_BeLOQ"), ]
pg_belowLOD <- pg[ which(pg$eDNA_eval_mean=="belowLOD"), ]
pg_NoCt <- pg[ which(pg$eDNA_eval_mean=="NoCt"), ]

#count for each of the subsets
ph_aboveLOQ<-count(pg_aboveLOQ, c("spc.year.smpltp"))
ph_AbLOD_BeLOQ<-count(pg_AbLOD_BeLOQ, c("spc.year.smpltp"))
ph_belowLOD<-count(pg_belowLOD, c("spc.year.smpltp"))
ph_NoCt<-count(pg_NoCt, c("spc.year.smpltp"))
#turn into data frames
pi_aboveLOQ<-as.data.frame(ph_aboveLOQ)
pi_AbLOD_BeLOQ<-as.data.frame(ph_AbLOD_BeLOQ)
pi_belowLOD<-as.data.frame(ph_belowLOD)
pi_NoCt<-as.data.frame(ph_NoCt)
#paste columns together
MONIS5eDNA09_df$spc.year.smpltp <- paste(MONIS5eDNA09_df$gen_specnm.year_inds,MONIS5eDNA09_df$smpltp,sep=".")
#check th sample types
unique(MONIS5eDNA09_df$smpltp)
#count number of species per year sample site 
c.sp.ye.splt <-as.data.frame(table(MONIS5eDNA09_df$spc.year.smpltp))
#check the species names match what they should
unique(MONIS5eDNA09_df$gen_specnm)
#check the number of replicates is not higher than what was used in the qPCR setups
unique(c.sp.ye.splt$Freq)
#match counts back with data frame
MONIS5eDNA09_df$freq_NoCt <- pi_NoCt$freq[match(MONIS5eDNA09_df$spc.year.smpltp, pi_NoCt$spc.year.smpltp)]
MONIS5eDNA09_df$freq_belowLOD <- pi_belowLOD$freq[match(MONIS5eDNA09_df$spc.year.smpltp, pi_belowLOD$spc.year.smpltp)]
MONIS5eDNA09_df$freq_AbLOD_BeLOQ <- pi_AbLOD_BeLOQ$freq[match(MONIS5eDNA09_df$spc.year.smpltp, pi_AbLOD_BeLOQ$spc.year.smpltp)]
MONIS5eDNA09_df$freq_aboveLOQ <- pi_aboveLOQ$freq[match(MONIS5eDNA09_df$spc.year.smpltp, pi_aboveLOQ$spc.year.smpltp)]
#replace NAs with zeros
MONIS5eDNA09_df$freq_NoCt[is.na(MONIS5eDNA09_df$freq_NoCt)] <- 0
MONIS5eDNA09_df$freq_belowLOD[is.na(MONIS5eDNA09_df$freq_belowLOD)] <- 0
MONIS5eDNA09_df$freq_AbLOD_BeLOQ[is.na(MONIS5eDNA09_df$freq_AbLOD_BeLOQ)] <- 0
MONIS5eDNA09_df$freq_aboveLOQ[is.na(MONIS5eDNA09_df$freq_aboveLOQ)] <- 0
#make a new column that shows the replicates quality levels
MONIS5eDNA09_df$freq_repl_eval <- as.character(paste("'",MONIS5eDNA09_df$freq_NoCt,"/", MONIS5eDNA09_df$freq_belowLOD,"/", MONIS5eDNA09_df$freq_AbLOD_BeLOQ,"/", MONIS5eDNA09_df$freq_aboveLOQ,"'", sep=""))
#add an empty column with just NAs to fil with color codings
MONIS5eDNA09_df[,"eDNA_eval_p_repl_col"] <- NA
#replace in the empty column, the order is important, as you otherwise will end up with the last evaluations
MONIS5eDNA09_df$eDNA_eval_p_repl_col[MONIS5eDNA09_df$freq_NoCt>=0] <- 1 #"white" #"NoCt" #0  
MONIS5eDNA09_df$eDNA_eval_p_repl_col[MONIS5eDNA09_df$freq_belowLOD>=1] <- 2 #"yellow" #"belowLOD" #1  
MONIS5eDNA09_df$eDNA_eval_p_repl_col[MONIS5eDNA09_df$freq_AbLOD_BeLOQ>=1] <- 3 # "azure3" #"AbLOD_BeLOQ" #2 
MONIS5eDNA09_df$eDNA_eval_p_repl_col[MONIS5eDNA09_df$freq_aboveLOQ>=1] <- 4 #"azure4" #"one aboveLOQ" #3 
MONIS5eDNA09_df$eDNA_eval_p_repl_col[MONIS5eDNA09_df$freq_aboveLOQ>=3] <- 5 #"black" #"all 3 above LOQ" #4
#add an empty column with just NAs to fil with evaluations
MONIS5eDNA09_df[,"eDNA_eval_p_repl_descr"] <- NA
#replace in the empty column, the order is important, as you otherwise will end up with the last evaluations
MONIS5eDNA09_df$eDNA_eval_p_repl_descr[MONIS5eDNA09_df$freq_NoCt>=0] <- "NoCt" #0  
MONIS5eDNA09_df$eDNA_eval_p_repl_descr[MONIS5eDNA09_df$freq_belowLOD>=1] <- "belowLOD" #1  
MONIS5eDNA09_df$eDNA_eval_p_repl_descr[MONIS5eDNA09_df$freq_AbLOD_BeLOQ>=1] <- "AbLOD_BeLOQ" #2 
MONIS5eDNA09_df$eDNA_eval_p_repl_descr[MONIS5eDNA09_df$freq_aboveLOQ>=1] <- "1aboveLOQ" #3 
MONIS5eDNA09_df$eDNA_eval_p_repl_descr[MONIS5eDNA09_df$freq_aboveLOQ>=3] <- "3aboveLOQ" #4 
#add an empty column with just NAs to fil with color codings
MONIS5eDNA09_df[,"eDNA_eval_t_repl_col"] <- NA
#replace in the empty column, the order is important, as you otherwise will end up with the last evaluations
MONIS5eDNA09_df$eDNA_eval_t_repl_col[MONIS5eDNA09_df$eDNA_eval_p_repl_col==1] <- "white" #"white" #"NoCt" #0  
MONIS5eDNA09_df$eDNA_eval_t_repl_col[MONIS5eDNA09_df$eDNA_eval_p_repl_col==2] <- "yellow" #"yellow" #"belowLOD" #1  
MONIS5eDNA09_df$eDNA_eval_t_repl_col[MONIS5eDNA09_df$eDNA_eval_p_repl_col==3] <- "orange" # "azure3" #"AbLOD_BeLOQ" #2 
MONIS5eDNA09_df$eDNA_eval_t_repl_col[MONIS5eDNA09_df$eDNA_eval_p_repl_col==4] <- "red" #"azure4" #"one aboveLOQ" #3 
MONIS5eDNA09_df$eDNA_eval_t_repl_col[MONIS5eDNA09_df$eDNA_eval_p_repl_col==5] <- "black" #"black" #"all 3 above LOQ" #4
#check the sample sites
unique(MONIS5eDNA09_df$smpltp)
#check the categories assigned
unique(MONIS5eDNA09_df$freq_repl_eval)
#make a column with AssayIDcode numbers
MONIS5eDNA09_df$AssayIDcode <- gsub("AssID","",MONIS5eDNA09_df$AssayID)
# get the number of the month . See this question on stackoverflow: https://stackoverflow.com/questions/22603847/how-to-extract-month-from-date-in-r
date02 <- as.factor(MONIS5eDNA09_df$Dato_inds)
MONIS5eDNA09_df$month_inds <- lubridate::month(as.POSIXlt(date02, format="%d/%m/%Y"))
# get the abbreviation for the month the sample was collected. See this website: https://stackoverflow.com/questions/22058393/convert-a-numeric-month-to-a-month-abbreviation/22063912
month_M09 <- as.numeric(MONIS5eDNA09_df$month_inds)
MONIS5eDNA09_df$month_inds2 <- month.abb[month_M09]
#The idea is to use these abbreviations later on in the maps, and plot the months sampled underneath each point sampled
#check the distribution of months sampled
hist(month_M09)
# The distribution of samples can be used later on for determining which months 
# to subset the data frame by, for making seasonal plots
#based on the range of months add a column with season categories
# start with an empty column, only with NAs
MONIS5eDNA09_df[,"season_cat"] <- NA
# based on evaluations of the months, add up with seasonal category
MONIS5eDNA09_df$season_cat[MONIS5eDNA09_df$month_inds<=6 ] <- "season_1" #season_1 is months numbered 1 to 6
MONIS5eDNA09_df$season_cat[MONIS5eDNA09_df$month_inds>=7 ] <- "season_2" #season_2 is months numbered 7 and above
# Now "season_1" equals the samples from the spring
# and "season_2" equals the samples from the fall


#count the number of season to loop over
no.of.seasons <- length(unique(MONIS5eDNA09_df$season_cat))
# make a sequence of numbers to use in a data frame
no_for_season <- seq(1:no.of.seasons)
#get the names of the seasons -  to use in the loop below
categories.of.seasons <- sort(unique(MONIS5eDNA09_df$season_cat))
# make names for the seasons
names.of.seasons <- c("spring","fall")
# bind to a data frame
seaons_nms_df <- as.data.frame(cbind(no_for_season,categories.of.seasons,names.of.seasons))
# make one of the columns numeric
seaons_nms_df$no_for_season <- as.numeric(seaons_nms_df$no_for_season)














####################################################################################
# Start Appendix B
####################################################################################


########################################################################################
#
# prepare maps with eDNA categories mapped on harbours
# In a previous version of this code I included a section in the middle, 
# that allows for plotting bars on the harbours
# to reflect the intensity of the eDNA levels. This middle section has been deleted
# but should be possible to append from some of my previous R code on eDNA levels
#
########################################################################################

#check if the object with the bathymetric map does not exists
if (!exists("dk.sea"))
{
  ## draw a map of the world, but limit the x-axis to between 8 to 14 ° E lon 
  ## and the y-axis to between 54 and 58 N° lat , 
  ## colour the land and sea in hexadecimal colors
  ## import bathymetric data for the region, 
  ## with lon1 and lon2 specifying the western end eastern 
  ## boundaries of the bathymetric plot 
  ## set the resolution to be 2, a higher resulotion takes longer to download
  # All maps are to be based on the same bathymetric map , so this can be reused
  dk.sea <- getNOAA.bathy(lon1 = 6, lon2 = 18,
                          lat1 = 54, lat2 = 58, resolution = 2)
  
  #close the if test above - i.e. if the 'dk.sea' object does not exist, then get it
}

# All maps are to be based on the same bathymetric map , so this can be reused
########################################################################################
#first get the species names, and Assay No's and assign numbers to each, to use 
# for numbering appendices

#get the unique species names
latspecnm <- unique(MONIS5eDNA09_df$gen_specnm)
#match the assay number to the data frame with species
AIfps <- MONIS3.ls.assays03_df$Assay_ID[match(latspecnm, MONIS3.ls.assays03_df$Lat_Species)]
#make a new data frame with assay Id No and species
nlspnm <- data.frame(AIfps,latspecnm)
#reorder by the column 'AssayIDNo'
nlspnm<- nlspnm[order(nlspnm$AIfps),]
#make a list of numbers for the unique species
no.latspc <- seq(1:length(latspecnm))
#add a new column with no to use for appendix numbering
nlspnm <- cbind(nlspnm, no.latspc) 
#use the new order of latin species names for producing plots
latspecnm <- unique(nlspnm$latspecnm)
#remove NA from the factors
#see this webpage: https://www.r-bloggers.com/r-drop-factor-levels-in-a-dataset/
latspecnm <- latspecnm[latspecnm!="NA"]
# loop over all species names in the unique list of species, and make plots. 
#Notice that the curly bracket ends after the pdf file is closed

#use a single species to check out if the loop works
#latspecnm <- "Mya arenaria"
#spec.lat <- "Mya arenaria"
########################################################################################
#loop over all species in data frame
for (spec.lat in latspecnm){
  print(spec.lat)
  #}
  #get the Danish commom name
  sbs.dk.nm <- MONIS3.ls.assays03_df$Common_name_Danish[match(spec.lat, MONIS3.ls.assays03_df$Lat_Species)]
  #replace the first underscore with a space in the dk common name
  sbs.dk.nm <- sub("_"," ",sbs.dk.nm)
  #get the AssayIDNo
  sbs.AssIDNo <- MONIS3.ls.assays03_df$Assay_ID[match(spec.lat, MONIS3.ls.assays03_df$Lat_Species)]
  #get the number for the appendix plot number
  no.spc.app.plot <- nlspnm$no.latspc[match(spec.lat, nlspnm$latspecnm)]
  #pad with zeros to two characters
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  AIfps <-stringr::str_pad(no.spc.app.plot, 2, pad = "0")
  #subset based on variable values, subset by species name
  sbs.MONIS5eDNA09_df <- MONIS5eDNA09_df[ which(MONIS5eDNA09_df$gen_specnm==spec.lat), ]
  #count using the plyr-package - see: https://www.miskatonic.org/2012/09/24/counting-and-aggregating-r/
  sbs.tot_smpl <- count(sbs.MONIS5eDNA09_df, c("Lok_omr01", "lok_pos_lon", "lok_pos_lat","year_inds"))
  #get the unique years
  year.dv <- unique(sbs.MONIS5eDNA09_df$year_inds)
  # get the upper and lower year in the range of the years with eDNA monitoring
  min.year <- min(year.dv)
  max.year <- max(year.dv)
  #make a line with years to use for the file name
  year_rng <- paste(min.year,"_to_",max.year,sep="")
  #count the number of years to use for the panes
  nppy <- length(year.dv)
  
  # use single years to test out loop below
  #year.dv <- "2017"
  #year.ind <- "2017"
  
  
  #make a pdf file nad give it a file name based on variables
  pdf(c(paste("App_B",AIfps,"_edna_lvls_MONIS5.AssayNo",sbs.AssIDNo,"_",spec.lat,"_",year_rng,"_02.pdf",  sep = ""))
      ,width=(2*1.6*8.2677),height=(nppy*1.6*2*2.9232))
  #set plotting margins
  op <-
    par(
      mfrow = c(nppy, 2), # set number of panes inside the plot - i.e. c(2,2) would make four panes for plots
      oma = c(1, 1, 0, 0), # set outer margin (the margin around the combined plot area) - higher numbers increase the number of lines
      mar = c(4, 4, 4, 2) # set the margin around each individual plot. with sides ordered as: " c(bottom,left,top,right)"
    )
  
  ##library for colours 
  library(RColorBrewer)
  #https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/
  #display.brewer.pal(9, "GnBu")
  blues <- rev(brewer.pal(9, "GnBu"))
  #dev.off()
  
  #loop over the years sampled
  for (year.ind in year.dv){
    print(year.ind)
    #} 
    #year.ind <- "2018"
    #subset based on variable values
    # subset among the years
    sbs.pery.MONIS5eDNA09_df <- sbs.MONIS5eDNA09_df[ which(sbs.MONIS5eDNA09_df$year_inds==year.ind ), ]
    
    ################################################################################
    # comment section below to out  to remove bathymetric isobars on maps.
    ################################################################################
    
    #to try out the loop assign only one category
    #categories.of.seasons <- "season_1"
    #loop over the seasons
    for (season in categories.of.seasons){
      print(season)
      #}
      #use match to match the season with a data frame and get the name for the season
      spcfc_seaon_name <- seaons_nms_df$names.of.seasons[match(season, seaons_nms_df$categories.of.seasons)]
      spcfc_seaon_name <- as.character(spcfc_seaon_name)
      
      # subset among the months to get from 1 to 6
      sbs.pery.MONIS5eDNA09_season_df <- sbs.pery.MONIS5eDNA09_df[ which(sbs.pery.MONIS5eDNA09_df$season_cat== season), ]
      #check if the data frame is empty . See this question : https://stackoverflow.com/questions/35366187/how-to-write-if-else-statements-if-dataframe-is-empty
      # check for multiple situations - first check if the data frame is empty - https://www.datamentor.io/r-programming/if-else-statement/
      if (dim(sbs.pery.MONIS5eDNA09_season_df)[1] == 0) {
        print(paste("data frame for",spcfc_seaon_name,year.ind,"is empty", sep=" "))
        #}
        #if subsetted data fram for spring is empty - no samples
        #then create an empty data frame with zeroes and no color for points
        #c_M09_season_no <- length(sbs.pery.MONIS5eDNA09_season_df)
        #sbs.pery.MONIS5eDNA09_season2_df <- rbind(sbs.pery.MONIS5eDNA09_season_df,(rep(0,c_M09_season_no)))
        #colnames(sbs.pery.MONIS5eDNA09_season2_df) <- colnames(sbs.pery.MONIS5eDNA09_season_df)
        #sbs.pery.MONIS5eDNA09_season_df <- sbs.pery.MONIS5eDNA09_season2_df
        #sbs.pery.MONIS5eDNA09_season_df$eDNA_eval_t_repl_col <- "#000000FF"
        #
        plot.new()
      }
      # check for multiple situations - second check if the data frame has more than one dimension - https://www.datamentor.io/r-programming/if-else-statement/
      # if the dataframe does have more than one dimension, then try and plot it
      else if (dim(sbs.pery.MONIS5eDNA09_season_df)[1] >= 1) 
        #start curly bracket for 'else if' test testing whether the dimensions on the data frame is more than one, if it is then make the plot on the map
      {
        #____________________________________________________________________________________
        # start -make the plot for the season
        #____________________________________________________________________________________
        #plot a color gradient for depth of the sea
        plot(dk.sea,  add=F, lwd = c(0, 0), col=transp_col, image = TRUE, bpal = c(alpha(blues, 0.4)),
             xlim = c(6, 18), ylim = c(54, 58),
             asp=1.4, #change from 1 to 1.6 for stretchin the map along the latitude
             cex=2.8, #vfont = c("sans serif", "plain"),
             las=1
        )
        #plot bathymetric lines
        ################################################################################
        # start  : comment section below to out  to remove bathymetric isobars on maps.
        ################################################################################
        plot(dk.sea, add=TRUE, lwd = c(0.8, 1.6), lty = c(1, 1),
             xlim = c(6, 18), ylim = c(54, 58),
             deep = c(-500, 0), shallow = c(-50, 0), step = c(50, 0),
             cex=2.8, vfont = c("sans serif", "plain"),
             asp=1.4, #change from 1 to 1.6 for stretchin the map along the latitude
             col = c("#00009B", "black"), drawlabels = c(TRUE, TRUE),
             las=1) # set labels horizontal to axis : see : https://stackoverflow.com/questions/1828742/rotating-axis-labels-in-r
        ################################################################################
        # end  : comment section below to out  to remove bathymetric isobars on maps.
        ################################################################################
        #plot land on map
        map('worldHires', add=TRUE, fill=TRUE, 
            xlim = c(6, 18), ylim = c(54, 58),
            col="grey",
            asp=1.4, #change from 1 to 1.6 for stretchin the map along the latitude
            bg=transp_col,
            las=1)
        
        #deduct a bit from the latitude, to lower the positioning of the label
        redlatpt1 <- (max(sbs.pery.MONIS5eDNA09_season_df$lok_pos_lat)-min(sbs.pery.MONIS5eDNA09_season_df$lok_pos_lat))/35.6
        #add text to the same lon lat position for harbour
        text(sbs.pery.MONIS5eDNA09_season_df$lok_pos_lon, 
             sbs.pery.MONIS5eDNA09_season_df$lok_pos_lat-redlatpt1, 
             sbs.pery.MONIS5eDNA09_season_df$month_inds2, 
             pos=1,
             cex= 0.9
        )
        #add a point to each position - with colour for the evaluation of the eDNA level
        #for eDNA evaluation
        points (sbs.pery.MONIS5eDNA09_season_df$lok_pos_lon, 
                sbs.pery.MONIS5eDNA09_season_df$lok_pos_lat, 
                pch = 22, 
                bg=c(alpha(c(as.character(sbs.pery.MONIS5eDNA09_season_df$eDNA_eval_t_repl_col)))),
                col="black", #set color of point
                lwd=2,
                #pos = 1,
                cex= 2.4)
        #add a title
        title(main = c(paste("eDNA detected in the",spcfc_seaon_name,"in",year.ind
                             , sep=" "))
              ,line=-2)
        #get the latin species name without underscore
        spec.lat.no_undersc <- paste(sub('_', ' ', spec.lat))
        # add legend for eDNA evaluation # see for placement: http://www.sthda.com/english/wiki/add-legends-to-plots-in-r-software-the-easiest-way
        legend("topright", "(x,y)", 
               bg="white",
               c("No Ct","below LOD", "above LOD and below LOQ" ,"1 above LOQ", "3 above LOQ"),
               ncol=1,
               pch = c(22,22, 22, 22, 22), #set type of point
               col= c("black", "black", "black", "black", "black"), #set color of point
               #pt.bg=c(alpha(c("white", "yellow", "black"), 0.6)), #set background color of point
               pt.bg=c(c("white", "yellow", "orange","red", "black")), #set background color of point
               pt.lwd=c(1.0),      
               title = "eDNA evalution ",
               cex=1.1,
               inset = 0.02)
        
        #end curly bracket for 'else if' test testing whether the dimensions on the data frame is more than one, if it is then make the plot on the map
      }
      
      #____________________________________________________________________________________
      # end -make the plot for season
      #____________________________________________________________________________________
      
      
      
      #end loop over seasons
    }
    #end loop over years
  }
  #apply the par settings for the plot as defined above.
  par(op)
  #add a title for the pdf plot
  #mtext(c(bquote('Appendix B'
  #         ~.(AIfps)~', '
  #         ~italic(.(spec.lat))~', '
  #         ~'('~.(sbs.dknm)~'), '
  #)), outer=TRUE,
  mtext(c(paste("Appendix B",
                AIfps," ",
                spec.lat,
                " (",sbs.dk.nm,")"),  sep = ""), outer=TRUE, 
        #use at , adj and padj to adjust the positioning
        #at=par("usr")[1]+0.15*diff(par("usr")[1:2]),
        adj=0.02,#1,
        padj=0,#2,
        #use side to place it in the top
        side=3, cex=1.2, line=-2.15)
  # end pdf file to save as
  dev.off()
  #
} #end loop over species
## 
# here the loop over all species names ends
####################################################################################
# End Appendix B
####################################################################################















####################################################################################
# Start Appendix C
####################################################################################

##########################################################################################
#
# Make a table that looks somewhat similar to Table 5 presented by :
#  Li, J, Hatton‐Ellis, TW, Lawson Handley, L‐J, et al. Ground‐truthing of a fish‐based environmental DNA metabarcoding method for assessing the quality of lakes. J Appl Ecol. 2019; 56: 1232– 1244. https://doi.org/10.1111/1365-2664.13352 
#
# Make a table for each year sampled
##########################################################################################
#copy the data frame
MONIS5eDNA10_df <- MONIS5eDNA09_df
#define the columns to keep 
keeps <- c("smpltp",
           "gen_specnm.year_inds",
           "eDNA_eval_p_repl_descr",
           "month_inds2",
           "season_cat")
#keep only selected columns
MONIS5eDNA11_df <- MONIS5eDNA10_df[keeps]
#count number of rows
nrow(MONIS5eDNA11_df)
#keep unique rows only
MONIS5eDNA12_df <- MONIS5eDNA11_df %>% dplyr::distinct(smpltp, gen_specnm.year_inds, eDNA_eval_p_repl_descr, .keep_all = TRUE)
# Sort by vector name [smpltp] then [gen_specnm.year_inds] # https://chartio.com/resources/tutorials/how-to-sort-a-data-frame-by-multiple-columns-in-r/
MONIS5eDNA11_df <- MONIS5eDNA11_df[
  with(MONIS5eDNA11_df, order(smpltp, gen_specnm.year_inds)),
  ]


# Sort by vector name [smpltp] then [gen_specnm.year_inds] # https://chartio.com/resources/tutorials/how-to-sort-a-data-frame-by-multiple-columns-in-r/
MONIS5eDNA12_df <- MONIS5eDNA12_df[
  with(MONIS5eDNA12_df, order(smpltp, gen_specnm.year_inds)),
  ]
# make a new column that fuses MST sample number together with month
MONIS5eDNA12_df$smpltp.month <- paste(MONIS5eDNA12_df$smpltp,".",MONIS5eDNA12_df$month_inds2,sep="")
#split column by delimiter, and turn in to data frame # https://www.rdocumentation.org/packages/splitstackshape/versions/1.4.8/topics/cSplit
MONIS5eDNA13_df <- as.data.frame(splitstackshape::cSplit(MONIS5eDNA12_df,"gen_specnm.year_inds", sep = "."))
#Rename specific column # see :  https://stackoverflow.com/questions/7531868/how-to-rename-a-single-column-in-a-data-frame
names(MONIS5eDNA13_df)[names(MONIS5eDNA13_df) == 'gen_specnm.year_inds_2'] <- 'yrs_smpl'
names(MONIS5eDNA13_df)[names(MONIS5eDNA13_df) == 'gen_specnm.year_inds_1'] <- 'gen_specnm'

#use the year listed in a vector previously
yrs <- year.dv
# use just a single year to start with for testing the loop
#yrs <- "2017"
#yr_smpl <- "2017"
#loop over years sampled -  to produce individual tables per year sampled
for (yr_smpl in yrs){
  print(yr_smpl)
  #}
  #subset based on variable values - only retain rows where the column that match the criteria 
  sbs.MO13y_df <- MONIS5eDNA13_df[ which(MONIS5eDNA13_df$yrs_smpl==yr_smpl), ]
  #to try out the loop assign only one category
  #categories.of.seasons <- "season_1"
  #loop over the seasons
  for (season in categories.of.seasons){
    print(season)
    #}
    #subset based on variable values - only retain rows where the column that match the criteria 
    sbs.MO14ym_df <- sbs.MO13y_df[ which(sbs.MO13y_df$season_cat==season), ]
    
    head(sbs.MO14ym_df,4)
    #define the columns to keep 
    keeps <- c("gen_specnm",
               "eDNA_eval_p_repl_descr",
               "smpltp.month")
    #keep only selected columns
    sbs.MO15ym_df <- sbs.MO14ym_df[keeps]
    
    #reshape the data frame to have smpls for columns
    sbs.MO16ym_df <- reshape(sbs.MO15ym_df, idvar = "gen_specnm", timevar = "smpltp.month", direction = "wide")
    #Replace characters in column names gsub : # https://stackoverflow.com/questions/39670918/replace-characters-in-column-names-gsub
    names(sbs.MO16ym_df) <- gsub(x = names(sbs.MO16ym_df), pattern = "eDNA_eval_p_repl_descr\\.", replacement = "")  
    #count the number of columns
    nc.MO16 <- ncol(sbs.MO16ym_df)
    
    
    #use match to match the season with a data frame and get the name for the season
    spcfc_seaon_name <- seaons_nms_df$names.of.seasons[match(season, seaons_nms_df$categories.of.seasons)]
    spcfc_seaon_name <- as.character(spcfc_seaon_name)
    #use match to match the season with a data frame and get the category number for the season
    spcfc_seaon_no <- seaons_nms_df$no_for_season[match(season, seaons_nms_df$categories.of.seasons)]
    spcfc_seaon_no <- as.numeric(spcfc_seaon_no)
    
    if (dim(sbs.MO16ym_df)[1] == 0) {
      print(paste("data frame for",spcfc_seaon_name,yr_smpl,"is empty", sep=" "))
      sbs.MO16ym_df <- as.data.frame(rbind(c("MST_smpl01","MST_smpl02"),c("no_data","sampled")))
    }
    
    #replace values in entire data frame
    sbs.MO16ym_df[sbs.MO16ym_df=="NoCt"]<-"NoCq"
    sbs.MO16ym_df[sbs.MO16ym_df=="belowLOD"]<-"bLOD" 
    sbs.MO16ym_df[sbs.MO16ym_df=="AbLOD_BeLOQ"]<-"aLODbLOQ"
    sbs.MO16ym_df[sbs.MO16ym_df=="1aboveLOQ"]<-"1aLOQ"
    sbs.MO16ym_df[sbs.MO16ym_df=="3aboveLOQ"]<-"3aLOQ"
    #get the unique years sampled
    yrs <- yr_smpl
    # Remove columns from dataframe where ALL values are NA # https://stackoverflow.com/questions/2643939/remove-columns-from-dataframe-where-all-values-are-na
    sbs.MO16ym_df <- sbs.MO16ym_df[,colSums(is.na(sbs.MO16ym_df))<nrow(sbs.MO16ym_df)]
    # delete the first column from the data frame
    #sbs.MO_df[,1] <- NULL
    
    sbs.MO_df <- sbs.MO16ym_df
    #https://cran.r-project.org/web/packages/tableHTML/vignettes/tableHTML.html
    if(!require(tableHTML)){
      install.packages("tableHTML")
      library(tableHTML)
    }
    require(tableHTML)
    #try the tableHTML with no border
    tableHTML <- sbs.MO_df %>% 
      tableHTML(border = 0) 
    #count the number of columns in the dataframe
    l.s.MO <- length(sbs.MO_df)
    #get unique cell values in dataframe : see : http://r.789695.n4.nabble.com/Retrieve-distinct-values-within-a-whole-data-frame-td1460205.html
    #apart from the first column
    unique(unlist(sbs.MO_df[2:l.s.MO]))
    #make lists of the words in the cells to color using the 'add_css_conditional_column' function
    eDNA.lvl01 <- c("NoCq") #white
    eDNA.lvl02 <- c("bLOD") #yellow
    eDNA.lvl03 <- c("aLODbLOQ") #orange
    eDNA.lvl04 <- c("1aLOQ") #red
    eDNA.lvl05 <- c("3aLOQ") #black
    #place the data frame in a tableHTML object
    tableHTML <- sbs.MO_df %>% 
      tableHTML()
    # for eDNA.lvl01 <- c("NoCq") #white
    words <- eDNA.lvl01 #<- c("NoCq") #white
    col.f.cell <- "white"
    for (word in words) {
      tableHTML <- tableHTML %>% 
        add_css_conditional_column(columns = 2:l.s.MO, #make it work on column 2 to the last column
                                   conditional = "contains",
                                   value = word,
                                   css = list(c("background-color"),
                                              c(col.f.cell)))
    }
    # for eDNA.lvl02 <- c("beLOD") #yellow
    words <- eDNA.lvl02
    col.f.cell <- "yellow"
    for (word in words) {
      tableHTML <- tableHTML %>% 
        add_css_conditional_column(columns = 2:l.s.MO, #make it work on column 2 to the last column
                                   conditional = "contains",
                                   value = word,
                                   css = list(c("background-color"),
                                              c(col.f.cell)))
    }
    # for eDNA.lvl03 <- c("abLODbeLOQ") #orange
    words <- eDNA.lvl03
    col.f.cell <- "orange"
    for (word in words) {
      tableHTML <- tableHTML %>% 
        add_css_conditional_column(columns = 2:l.s.MO, #make it work on column 2 to the last column
                                   conditional = "contains",
                                   value = word,
                                   css = list(c("background-color"),
                                              c(col.f.cell)))
    }
    # for eDNA.lvl04 <- c("1abLOQ") #red
    words <- eDNA.lvl04
    col.f.cell <- "red"
    for (word in words) {
      tableHTML <- tableHTML %>% 
        add_css_conditional_column(columns = 2:l.s.MO, #make it work on column 2 to the last column
                                   conditional = "contains",
                                   value = word,
                                   css = list(c("background-color"),
                                              c(col.f.cell)))
    }
    # for eDNA.lvl05 <- c("3abLOQ") #black
    words <- eDNA.lvl05
    col.f.cell <- "black" # use this color for the cell
    col.f.font <- "white" #use this color for the font
    for (word in words) {
      tableHTML <- tableHTML %>% 
        add_css_conditional_column(columns = 2:l.s.MO, #make it work on column 2 to the last column
                                   conditional = "contains",
                                   value = word,
                                   css = list(c("background-color","color"),
                                              c(col.f.cell,col.f.font)))
    }
    
    t.HTML17 <- tableHTML
    
    
    #pad with zeros to two characters
    #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
    no.e4 <-stringr::str_pad(spcfc_seaon_no, 2, pad = "0")
    #head(t.HTML17)
    #and to export in a file
    write_tableHTML(t.HTML17, file = paste("App_C_",yr_smpl,"_",no.e4,"_",spcfc_seaon_name,"_table_eDNA_evalu.html", sep=""))
    
    #end loop over seasons
  }
  #end loop over years
}

####################################################################################
# End Appendix C
####################################################################################










####################################################################################
# Start Appendix D
####################################################################################
#
# Make html tables with color coding for evaluation of eDNA levels per category
# Make a table for each year sampled
#
####################################################################################
#copy the data frame
MONIS5eDNA18_df <- MONIS5eDNA09_df
#head(MONIS5eDNA18_df,9)
#define the columns to keep 
keeps <- c("smpltp",
           "gen_specnm.year_inds",
           "freq_repl_eval",
           "month_inds2",
           "season_cat")
#keep only selected columns
MONIS5eDNA19_df <- MONIS5eDNA18_df[keeps]
#count number of rows
nrow(MONIS5eDNA19_df)
#keep unique rows only
MONIS5eDNA20_df <- MONIS5eDNA19_df %>% dplyr::distinct(smpltp, gen_specnm.year_inds, freq_repl_eval, .keep_all = TRUE)
# Sort by vector name [smpltp] then [gen_specnm.year_inds] # https://chartio.com/resources/tutorials/how-to-sort-a-data-frame-by-multiple-columns-in-r/
MONIS5eDNA19_df <- MONIS5eDNA19_df[
  with(MONIS5eDNA19_df, order(smpltp, gen_specnm.year_inds)),
  ]
# Sort by vector name [smpltp] then [gen_specnm.year_inds] # https://chartio.com/resources/tutorials/how-to-sort-a-data-frame-by-multiple-columns-in-r/
MONIS5eDNA20_df <- MONIS5eDNA20_df[
  with(MONIS5eDNA20_df, order(smpltp, gen_specnm.year_inds)),
  ]
# make a new column that fuses MST sample number together with month
MONIS5eDNA20_df$smpltp.month <- paste(MONIS5eDNA20_df$smpltp,".",MONIS5eDNA20_df$month_inds2,sep="")
#split column by delimiter, and turn in to data frame # https://www.rdocumentation.org/packages/splitstackshape/versions/1.4.8/topics/cSplit
MONIS5eDNA21_df <- as.data.frame(splitstackshape::cSplit(MONIS5eDNA20_df,"gen_specnm.year_inds", sep = "."))
#Rename specific column # see :  https://stackoverflow.com/questions/7531868/how-to-rename-a-single-column-in-a-data-frame
names(MONIS5eDNA21_df)[names(MONIS5eDNA21_df) == 'gen_specnm.year_inds_2'] <- 'yrs_smpl'
names(MONIS5eDNA21_df)[names(MONIS5eDNA21_df) == 'gen_specnm.year_inds_1'] <- 'gen_specnm'
#use the year listed in a vector previously
yrs <- year.dv
# use just a single year to start with for testing the loop
#yrs <- "2017"
#yr_smpl <- "2017"
#loop over years sampled -  to produce individual tables per year sampled
for (yr_smpl in yrs){
  print(yr_smpl)
  #}
  #subset based on variable values - only retain rows where the column that match the criteria 
  sbs.MO21y_df <- MONIS5eDNA21_df[ which(MONIS5eDNA21_df$yrs_smpl==yr_smpl), ]
  #to try out the loop assign only one category
  #categories.of.seasons <- "season_1"
  #loop over the seasons
  for (season in categories.of.seasons){
    print(season)
    #}
    #subset based on variable values - only retain rows where the column that match the criteria 
    sbs.MO22ym_df <- sbs.MO21y_df[ which(sbs.MO21y_df$season_cat==season), ]
    
    #head(sbs.MO21ym_df,4)
    #define the columns to keep 
    keeps <- c("gen_specnm",
               "freq_repl_eval",
               "smpltp.month")
    #keep only selected columns
    sbs.MO23ym_df <- sbs.MO22ym_df[keeps]
    
    #reshape the data frame to have smpls for columns
    sbs.MO24ym_df <- reshape(sbs.MO23ym_df, idvar = "gen_specnm", timevar = "smpltp.month", direction = "wide")
    #Replace characters in column names gsub : # https://stackoverflow.com/questions/39670918/replace-characters-in-column-names-gsub
    names(sbs.MO24ym_df) <- gsub(x = names(sbs.MO24ym_df), pattern = "freq_repl_eval\\.", replacement = "")  
    #count the number of columns
    nc.MO24 <- ncol(sbs.MO24ym_df)
    #use match to match the season with a data frame and get the name for the season
    spcfc_seaon_name <- seaons_nms_df$names.of.seasons[match(season, seaons_nms_df$categories.of.seasons)]
    spcfc_seaon_name <- as.character(spcfc_seaon_name)
    #use match to match the season with a data frame and get the category number for the season
    spcfc_seaon_no <- seaons_nms_df$no_for_season[match(season, seaons_nms_df$categories.of.seasons)]
    spcfc_seaon_no <- as.numeric(spcfc_seaon_no)
    
    if (dim(sbs.MO24ym_df)[1] == 0) {
      print(paste("data frame for",spcfc_seaon_name,yr_smpl,"is empty", sep=" "))
      sbs.MO24ym_df <- as.data.frame(rbind(c("MST_smpl01","MST_smpl02"),c("no_data","sampled")))
    }
    #head(sbs.MO24ym_df,9)
    #get the number of columns
    nc.MO24 <- ncol(sbs.MO24ym_df)
    # get unique values in the data frame from a selected range of columns # see : https://stackoverflow.com/questions/40003028/extracting-unique-values-from-data-frame-using-r
    unqval.MO24 <- unique(unlist((sbs.MO24ym_df)[,2:nc.MO24]))
    #remove NA from vector - this will exclde the non-evaluated categories
    unqval.MO24 <- unqval.MO24[!is.na(unqval.MO24)]
    
    # With the grep function in R the different elements in 'freq_repl_eval' 
    # can be categorized for 
    # identification later on in the preparation of the html table that is to show
    # the colored categories for eDNA levels detected.
    # The coding for the elements in 'freq_repl_eval' are:
    # ' no of replicates with no Ct/no of replicates below LOD /no of replicates above LOD but below LOQ / no of replicates above LOQ ' 
    # the color coding for these elements in 'freq_repl_eval' are
    # ' white/yellow /orange / red or black '
    # red if a minimum if 1 replicate is above LOQ (disregarding if any lower levels are detected)
    # black if all replicates are above LOQ (this will automatically equal all lower categories being zero)
    # To try out the grep function, I have here below tried grepping for different elements in the list
    #grep for elements that begin with '0
    grep("^'0", unqval.MO24, value=T)
    #grep for elements that end with 0'
    grep("0'$", unqval.MO24, value=T)
    # grep for elements starting '0/0/0/ and then 1 or any higher number
    b_cat <- grep("'0/0/0/[1-9]+", unqval.MO24, value=T) # this will equal the black category with all replicates amplyfiying
    # then grep for all elements with not zero in the last category
    br_cat <- grep("'[0-9]+/[0-9]+/[0-9]+/[1-9]+", unqval.MO24, value=T) # this will equal both the red and balck category
    # find the difference between these two vectors - i.e. subtract the black category from the fused red-black category
    r_cat <- setdiff(br_cat,b_cat)
    # grep for all elements with not zero in the first category
    w_cat <- grep("'[1-9]+/0/0/0'", unqval.MO24, value=T) # this will equal all replicates not amplifying - i.e. this will equal the white category
    # find the difference between two vectors
    wyo_cat <- setdiff(unqval.MO24,br_cat)
    # find the difference between two vectors
    yo_cat <- setdiff(wyo_cat,w_cat)
    # grep for all elements with not zero in the last and second last category
    y_cat <- grep("'[0-9]+/[1-9]+/0/0'", unqval.MO24, value=T)
    # grep for all elements with not zero in thelast category
    o_cat <- grep("'[0-9]+/[0-9]+/[1-9]+/0'", unqval.MO24, value=T)
    #these categories are used below to identify the color coding in the html tables
    
    sbs.MO_df <- sbs.MO24ym_df
    #https://cran.r-project.org/web/packages/tableHTML/vignettes/tableHTML.html
    if(!require(tableHTML)){
      install.packages("tableHTML")
      library(tableHTML)
    }
    require(tableHTML)
    #try the tableHTML with no border
    tableHTML <- sbs.MO_df %>% 
      tableHTML(border = 0) 
    #count the number of columns in the dataframe
    l.s.MO <- length(sbs.MO_df)
    #get unique cell values in dataframe : see : http://r.789695.n4.nabble.com/Retrieve-distinct-values-within-a-whole-data-frame-td1460205.html
    #apart from the first column
    unique(unlist(sbs.MO_df[2:l.s.MO]))
    #make lists of the words in the cells to color using the 'add_css_conditional_column' function
    class(eDNA.lvl01)
    eDNA.lvl01 <- w_cat #c("/0/0/0'") #white
    eDNA.lvl02 <- y_cat #c("'0/") #yellow
    eDNA.lvl03 <- o_cat #c("aLODbLOQ") #orange
    eDNA.lvl04 <- r_cat #c("/1'") #red
    eDNA.lvl05 <- b_cat #c("'0/0/0/") #black
    #place the data frame in a tableHTML object
    tableHTML <- sbs.MO_df %>% 
      tableHTML()
    # for eDNA.lvl01 <- c("NoCq") #white
    words <- eDNA.lvl01 #<- c("NoCq") #white
    col.f.cell <- "white"
    for (word in words) {
      tableHTML <- tableHTML %>% 
        add_css_conditional_column(columns = 2:l.s.MO, #make it work on column 2 to the last column
                                   conditional = "contains",
                                   value = word,
                                   css = list(c("background-color"),
                                              c(col.f.cell)))
    }
    # for eDNA.lvl02 <- c("beLOD") #yellow
    words <- eDNA.lvl02
    col.f.cell <- "yellow"
    for (word in words) {
      tableHTML <- tableHTML %>% 
        add_css_conditional_column(columns = 2:l.s.MO, #make it work on column 2 to the last column
                                   conditional = "contains",
                                   value = word,
                                   css = list(c("background-color"),
                                              c(col.f.cell)))
    }
    # for eDNA.lvl03 <- c("abLODbeLOQ") #orange
    words <- eDNA.lvl03
    col.f.cell <- "orange"
    for (word in words) {
      tableHTML <- tableHTML %>% 
        add_css_conditional_column(columns = 2:l.s.MO, #make it work on column 2 to the last column
                                   conditional = "contains",
                                   value = word,
                                   css = list(c("background-color"),
                                              c(col.f.cell)))
    }
    # for eDNA.lvl04 <- c("1abLOQ") #red
    words <- eDNA.lvl04
    col.f.cell <- "red"
    for (word in words) {
      tableHTML <- tableHTML %>% 
        add_css_conditional_column(columns = 2:l.s.MO, #make it work on column 2 to the last column
                                   conditional = "contains",
                                   value = word,
                                   css = list(c("background-color"),
                                              c(col.f.cell)))
    }
    # for eDNA.lvl05 <- c("3abLOQ") #black
    words <- eDNA.lvl05
    col.f.cell <- "black" # use this color for the cell
    col.f.font <- "white" #use this color for the font
    for (word in words) {
      tableHTML <- tableHTML %>% 
        add_css_conditional_column(columns = 2:l.s.MO, #make it work on column 2 to the last column
                                   conditional = "contains",
                                   value = word,
                                   css = list(c("background-color","color"),
                                              c(col.f.cell,col.f.font)))
    }
    
    t.HTML17 <- tableHTML
    
    #pad with zeros to two characters
    #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
    no.e4 <-stringr::str_pad(spcfc_seaon_no, 2, pad = "0")
    #head(t.HTML17)
    #and to export in a file
    write_tableHTML(t.HTML17, file = paste("App_D_",yr_smpl,"_",no.e4,"_",spcfc_seaon_name,"_table_eDNA_evalu.html", sep=""))
    
    #end loop over seasons
  }
  #end loop over years
}






####################################################################################
# End Appendix D
####################################################################################











#