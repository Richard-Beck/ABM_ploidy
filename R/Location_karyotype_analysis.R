



#This code will be built upon to analyze Cel population  data from ABM model.
#Load library the following libraries separately before running the main code 


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



#Set working direction to load data from
setwd("/Users/john/Documents/IMO_stuff/ploidy_abm")
myData<-list() #Container for data
file_list<-list.files(path = "/Users/john/Documents/IMO_stuff/ploidy_abm",pattern = "*karyotype_Location.*\\.csv$")
#Re-order list of files to match timepoints
file_list<-file_list[order(as.numeric(sub("([0-9]*).*","\\1",file_list)))]
for(i in 1:length(file_list)){
  myData[[i]]<-as.matrix(read.csv(file_list[i],header=FALSE))
}

totalTimpoints<-length(myData)

#dev.off().  To reset plot


#Extract nonzero rows from large matrix
nonZeroRows<-list()
for(k in 1:length(myData)){
 m=c()
for(i in 1:nrow(myData[[k]])){
  #if both x and y coordinates are zero, then there is no cell there, remove
  if(myData[[k]][i,1]==0 && myData[[k]][i,2]==0){
    #Collect zero indexes
    m<-c(m,i)
  }
}
 #Remove zero indexed rows
nonZeroRows[[k]]<-myData[[k]][-m,]
}

  



#Rename columns (X-coord, y-coord,copy numbers)
columnNames<-c("x","y","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10",
               "C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22")


for(k in 1:length(nonZeroRows)){
  colnames(nonZeroRows[[k]])=columnNames
}




#Extract diploid cells
diploid<-c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
diploidCellsAtTimeK=list()

for(k in 1:length(nonZeroRows)){
  n=c()
  for(i in 1:nrow(nonZeroRows[[k]])){
     if(identical(as.numeric(nonZeroRows[[k]][i,3:24]),diploid)){
    #Collect diploid indexes
     n<-c(n,i)
    
   }
  }
  diploidCellsAtTimeK[[k]]<-nonZeroRows[[k]][n,]
}





#create a matrix  A(xi,yi) =d_i where d is the manhattan distance of karyotype ki from diploid state
manhattanDistance <- function(vect1, vect2){
  dist <- abs(vect1 - vect2)
  dist <- sum(dist)
  return(dist)
}


# Distance between karyotype and diploid : [x-coord, y-coord, dist]
#@TODO modify this function below to use dist(, method="manhattan")
distMatrix<-list()
for(k in 1:length(nonZeroRows)) {
  distMatrix[[k]] <- matrix(, nrow = nrow(nonZeroRows[[k]]), ncol = 3)
  #Extract all diploid cells
  for (i in 1:nrow(nonZeroRows[[k]])) {
    distMatrix[[k]][i, 1:2] = nonZeroRows[[k]][i, 1:2] #copying x,y coordinates
    distMatrix[[k]][i, 3] = manhattanDistance(nonZeroRows[[k]][i, 3:24], diploid)
  }
}



#Represent data as a matrix/grid whose (i,j) component is the distance from diploid
distfromDiploidMatrix<-list()

for(k in 1:length(nonZeroRows)){
  distfromDiploidMatrix[[k]]<-matrix(,nrow=100,ncol=100)
for(i in 1:nrow(distMatrix[[k]])){
  distfromDiploidMatrix[[k]][distMatrix[[k]][i,1],distMatrix[[k]][i,2]]=distMatrix[[k]][i,3]
}
}

#replace N/A in matrix with 0
for(k in 1:length(distfromDiploidMatrix)){
distfromDiploidMatrix[[k]]<- distfromDiploidMatrix[[k]] %>% replace(is.na(.), 0)
}
#an.distfromDiploidMatrix[distfromDiploidMatrix]<-0






#Compute physical distance between cells (coordinates) 
#spatialDistBtnCells is a lower triangular matrix whose ith column comprise the distance between cell i and cells i+1, i+2, 
# and so on repectively

spatialDistBtnCells<-list()
for(k in 1:length(nonZeroRows)){
spatialDistBtnCells[[k]]<-dist(nonZeroRows[[k]][,1:2], method = "euclidean", diag = TRUE, upper = FALSE)
spatialDistBtnCells[[k]]<-as.matrix(spatialDistBtnCells[[k]])
}

#Distance between cell by karyotype vectors 
#Same as distance from diploid cells above, except that this computes all pairwaise distances between karyotypes
karyotypeDistBtnCells<-list()
for(k in 1:length(nonZeroRows)){
  karyotypeDistBtnCells[[k]]<-as.matrix(dist(nonZeroRows[[k]][,3:24], method = "manhattan", diag = TRUE, upper = FALSE))
#karyotypeDistBtnCells[[k]]<-as.matrix(karyotypeDistBtnCells[[k]])
}

#plot karyotype distance at timpoint T 
plot(karyotypeDistBtnCells[[10]][1:17,1], main="Cell 1 vs population", xlab="Cells", ylab="Karyotype Dist.") 
#Timepoint 10: Plot of  karyotype distance btn cell 1  and all other cells

#Timepoint 10: plot of karyotype distances vs spatial distance between cell x and all other cells at timepoint y
plot(karyotypeDistBtnCells[[10]][1:80,1],spatialDistBtnCells[[10]][1:80,1],pch=c(8),col="red",
     xlab="cell 1- to- population:  Karyotype Dist.",ylab="Cell 1 -to-population: Spatial distance",
     main="Cell 1 against Population: Karytype distance vs spacial distance plot")




#populate the various cell types

#distinct Karyotypes
distinctKaryotypesPerTimepoint<-list()

#count the number of distinct karyotypes per timepoint
for(k in 1:length(nonZeroRows)){
  distinctKaryotypesPerTimepoint[[k]]  <-unique(nonZeroRows[[k]][,3:24])
}


# Count the number of each distinc karyotype per Timepoint
frequencyOfDistinctKaryotypesPerTimepoint<-list()
#Intege vector where ith value is the  total number karyotypes present at timepoint t=i
numberOfkaryotypesPerTimepoint<-matrix(0,nrow=1, ncol=length(nonZeroRows))
for(k in 1: length(nonZeroRows)){
  df<-as.data.frame(nonZeroRows[[k]][,3:24])
  #counting how many of each unique karyotype there are at each timepoint
   temp<-collapse::fcount(df)      #temp<-df %>% group_by_all() %>% count
  # Count is natural placed as last column, this swtching count vector to column 11
  frequencyOfDistinctKaryotypesPerTimepoint[[k]]<-temp %>% select(N, everything())
  numberOfkaryotypesPerTimepoint[1,k]<-nrow(frequencyOfDistinctKaryotypesPerTimepoint [[k]])
}





columnNames<-c("Count","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10",
               "C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22")

#colnames( countOfDistinctKaryotypesPerTimepoint)=columnNames


allKaryotypes<-as.matrix(frequencyOfDistinctKaryotypesPerTimepoint[[1]])
for(k in 2:length(frequencyOfDistinctKaryotypesPerTimepoint )){
  allKaryotypes<-rbind(allKaryotypes,as.matrix(frequencyOfDistinctKaryotypesPerTimepoint [[k]]))
}

allKaryotypes<-unique(allKaryotypes[,2:23])
allKaryotypes
#Total number of Karyotypes THAT HAS OCCURED DURING simulation
totalNumberOfKaryotypesInDataset<-nrow(allKaryotypes)






countOfEactKaryotypeOverTime<-matrix(0, nrow=nrow(allKaryotypes),ncol=length(frequencyOfDistinctKaryotypesPerTimepoint ))
for(j in 1:nrow(allKaryotypes)){
     for(k in 1:length(frequencyOfDistinctKaryotypesPerTimepoint )){
      n<-which(apply(frequencyOfDistinctKaryotypesPerTimepoint[[k]][,2:23], 1, function(x) return(all(x ==allKaryotypes[j,] ))))
      if(!identical(n,integer(0))){ # If index is not empty, then perform the subseting below
      countOfEactKaryotypeOverTime[j,k]=as.numeric(frequencyOfDistinctKaryotypesPerTimepoint[[k]][n,1])
      }
  }
}

#Replace NA with 0
countOfEactKaryotypeOverTime<- countOfEactKaryotypeOverTime %>% replace(is.na(.), 0)
columnNamesForTimepoints<-paste0("T",c(1:totalTimpoints))
rowNmaesForKaryotypes=paste0("K",c(1:nrow(allKaryotypes)))
#Time<-seq(from=1,to=totalTimpoints)
colnames(countOfEactKaryotypeOverTime)= columnNamesForTimepoints #rbind(Time[1:totalTimpoints],countOfEactKaryotypeOverTime)
rownames(countOfEactKaryotypeOverTime)=rowNmaesForKaryotypes

plot(countOfEactKaryotypeOverTime[4,],xlab="Time",ylab="Population", 
     main=paste("Evolutionof Karyotype:",lapply(paste0(as.character(allKaryotypes[4,]),collapse=" "),
     function(ii) gsub(" ","",ii))),,col="blue")


heatmap(countOfEactKaryotypeOverTime, main="Karyotye Emergence and Evolution")

#convert countOfEactKaryotypeOverTime to "long " matrix to use in ggplot
#For data frame output use melt
countDataFrame<-melt(countOfEactKaryotypeOverTime)
colnames(countDataFrame)<-c("Karyotype","Time","Value")

#Or unvomment the line below for other forms
#dat<-matrix(countOfEactKaryotypeOverTime, dimnames=list(t(outer(colnames(countOfEactKaryotypeOverTime), 
                                                         #  rownames(countOfEactKaryotypeOverTime), FUN=paste)), NULL))
#countDataFrame<-as.data.frame(dat)



Time=c(1:nrow(t(countOfEactKaryotypeOverTime)))
df<-as.data.frame(cbind(time,t(countOfEactKaryotypeOverTime)))

# plot populations over time
ggplot(df, aes(x=Time)) +
  geom_smooth(mapping=aes(y=K1,colour="K1"))+
  labs(x="Time", y="Size",title="Diploid ppopulation")

ggplot(df, aes(x=Time)) +
  geom_point(mapping=aes(y=K6,colour="K7"))+
  geom_point(mapping=aes(y=K2,colour="K2"))+
  geom_point(mapping=aes(y=K4,colour="K4"))+
  geom_point(mapping=aes(y=K10,colour="K10"))+
  geom_point(mapping=aes(y=K10,colour="K10"))+
  labs(x="Time", y="Size",title="Heterogeneity in population")
  #scale_fill_gradient(low = "purple", high = "red")+
  


#Bar chart /Histogram
ggplot(countDataFrame,aes(x=Time,y=Value,fill=Karyotype))+
  geom_bar(stat='identity')+
  labs(x="Time", ylab="Size", main="Heterogeneity in population over time")



#allPlots<- ggarrange(p1,p2,p3,p4,p5,p6,
#                     labels = c("A", "B", "C"),
#                     ncol = 2, nrow = 3)
# allPlots




