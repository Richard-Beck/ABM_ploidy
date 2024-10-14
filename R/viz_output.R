setwd("~/projects/008_birthrateLandscape/ABM_ploidy/output/test/run_2/")
#x <- read.csv("summary.csv")
o <- read.csv("oxygen.csv")
o$day <- round(o$time)

smm_o <- aggregate(list(nsteps=o$day),by=list(day=o$day),length)

p <- ggplot(o[o$time<2,],aes(x=time,y=O2))+
  geom_line()
p
p <- ggplot(o,aes(x=time,y=O2))+
  geom_line()
p
p <- ggplot(smm_o,aes(x=day,y=nsteps))+
  geom_line()+
  scale_y_log10()
p

library(ggplot2)
#p <- ggplot(x,aes(x=Time,y=NCells_R))+
 # geom_line()
#p
