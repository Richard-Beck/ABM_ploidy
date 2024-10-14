setwd("~/projects/008_birthrateLandscape/ABM_ploidy/output/deltas/")

ff <- list.files()

get_dt <- function(fi){
  fi <- head(unlist(strsplit(fi,split=".csv")),1)
  as.numeric(head(unlist(strsplit(fi,split="_")),1))
}

tt <- sapply(ff,get_dt)
ff <- ff[order(tt)]
tt <- tt[order(tt)]

x1 <- read.csv(ff[1],header=F)
x2 <- read.csv(ff[2],header=F)
x2$dv <- x2$V3-x1$V3
x2$g5 <- x2$V3>1
library(ggplot2)
p <- ggplot(x2,aes(x=V1,y=V2,fill=dv))+
  geom_raster()
p
p <- ggplot(x2,aes(x=V1,y=V2,fill=g5))+
  geom_raster()
p

stop()

ff <- ff[tt%%20==0]
tt <- tt[tt%%20==0]


df <- do.call(rbind,pbapply::pblapply(ff,function(fi){
  xi <- read.csv(fi,header=F)
  data.frame(max_delta = max(abs(xi$V3)),tot_delta = sum(xi$V3),ntot=sum(xi$V4)) 
}))

df$my_delta <- c(0,diff(df$ntot))


