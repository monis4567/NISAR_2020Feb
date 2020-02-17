# NISAR_2020Feb
## Plot eDNA levels on map 


## Code developed during NISAR project (2019-2020) by Steen W. Knudsen at NIVA-DK
####### R code for analysing filtered and extracted eDNA samples collected by the Danish environmental agency (Miløstyrelsen) (MST) over 2017-2018.

####### The water samples collected by MST are to be analysed for eDNA levels from 20 non indeginuous marine species in Danish seas.

####### The 20 non indeginuous marine species targeted for analysis are the same 20 species that the MONIS2, MONIS3 and MONIS4 project, carried out at NIVA-DK
#
## This code is prepared for analysis of the qPCR data obtained in the MONIS5 project by laboratorial work carried out in 2018-2019 by SWK
#
####### The overall purpose of the present code is to analyse the eDNA levels inferred in the qPCR setups performed over 2018-2019.


## ## This code contains the following sections:
## 1 - Preparation of standard curves with qPCR eDNA levels plotted on to curves. 
#######   These plots are prepared as copy numbers versus cycle treshold (Ct) values
## 2 - Plot of sampling locations on maps, for each species, with indication of eDNA intensity for each location monitored
## 3a - Make tables of evaluations of eDNA levels with categories inferred from limit of detection (LOD) and limit of quantification (LOQ)
## 3b - Make tables of evaluations of eDNA replicates analysed in qPCR set ups. With categories assigned from ealier inferred limit of detection (LOD) and limit of quantification (LOQ).

## A note on input files
####### The required input files can be used as examples for how input files are supposed to be organized and prepared
####### These input files stems from eDNA data collected over 2018-2019.
####### The R code provided will need these input files in order to be able to produce the tables and pdf plots
####### The tables generated with this R code are html files and can be viewed in a web browser

#
