
#This script curates blood vessel and karyotype data from spatial transcriptomics to be used as input for ABM
library(ggplot2)
library(data.table)
library(tidyverse)
library(dplyr)
library("ggpubr")
library(rgl)
library("spatstat.geom") 
library(plotly)
#install.packages("scatterplot3d") # Install
library("scatterplot3d") # load
install.packages("collapse")
library(reshape2)




################################### Blood vessel Data Curation ########################
#This part prepares vessel data 
vessels<-as.matrix(read.csv("/Users/john/Documents/IMO/ploidy_abm/blood_vessels.csv",header=FALSE))

#Scale values to align with ABM grid
maxi<-max(vessels) 
if(maxi>100){
vessels<-vessels/maxi*100 
}
vessels<-round(vessels, digits=0) 
vesses<-as.table(vessels)

#write to txt file
write.table(vessels,"~/Documents/IMO/ploidy_abm/vessels.txt",sep="\t",quote=F,col.names=F, row.names=F)



################################### Karyotype Data Curation ########################
#This part prepares karyotype data 
cellKaryotype<-as.matrix(read.csv("/Users/john/Documents/IMO/ploidy_abm/input/karyotype.csv",header=FALSE))

#Check for unique locations

#Find distinct location, (cells cannot be stacked)
distinctcells<-cellKaryotype[!duplicated(cellKaryotype[,1:2]),]

#Scale coordinates to match grid dimensions
maximum<-max(distinctcells[,1:2])
if(maximum>0){
  distinctcells[,1]<-(distinctcells[,1]/maximum)*100
  distinctcells[,2]<-(distinctcells[,2]/maximum)*100
}
#Round coordinate values to integers
distinctcells<-floor(distinctcells)

#Write to txt 
distinctcells<-as.table(distinctcells)
write.table(distinctcells,"~/Documents/IMO/ploidy_abm/input/karyotype.txt",sep="\t",quote=F,col.names=F, row.names=F)
