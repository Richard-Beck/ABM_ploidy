Test Homeostasis
================
Richard J Beck
2024-10-21

This example demonstrates setup and run for a sweep where 2 parameters
(normal_oxygen_consumption_rate, growth_michaelis_menten_constant) are
simultaneously varied.

For this sweep the configuration file is set up as demonstrated in the
test_diffusion example. Now, however the argument $sweeptype=x$ can be
set. This enables testing all combinations of the optional arguments, in
this case we do:

normal_oxygen_consumption_rate=100:500:100

growth_michaelis_menten_constant=0.0005:0.003:0.0005

After edit the config file so that the paths correspond to the those on
your machine, the following command can be run:

    Rscript .\scripts\runSweep.R .\examples\test_homeostasis\config.txt

After the sweep has finished running, we can analyze the results.
Overall it looks like homeostasis was not reached in any conditions and
we will have to rerun the sweep with a longer time interval!

``` r
library(ggplot2)
path2sweepdata <- "output/test_homeostasis/data/"
runs <- list.files(path2sweepdata)

pars <- c("normal_oxygen_consumption_rate","growth_michaelis_menten_constant")

df <- do.call(rbind,lapply(runs,function(r){
  o2 <- read.csv(paste0(path2sweepdata,r,"/oxygen.csv"))
  parfile <- readLines(paste0(path2sweepdata,r,"/params.txt"))
  for(parName in pars){
    i <- grepl(parName,parfile)
    par <- parfile[i]
    par <- unlist(strsplit(par,split=" "))
    o2$parval <- as.numeric(par[2])
    colnames(o2)[colnames(o2)=="parval"] <- par[1]
  }
  o2 <- split(o2,f=round(o2$time))
  do.call(rbind,lapply(o2,function(oi){
    tmp <- tail(oi,1)
    tmp$Nsteps <- nrow(oi)
    tmp
  }))
}))


p <- ggplot(df,aes(x=time,y=O2,color=normal_oxygen_consumption_rate,
                   group=normal_oxygen_consumption_rate))+
  facet_grid(cols=vars(growth_michaelis_menten_constant))+
  geom_line()+
  scale_color_viridis_c()
p
```

![](README_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
p <- ggplot(df,aes(x=time,y=Nsteps,color=normal_oxygen_consumption_rate,
                   group=normal_oxygen_consumption_rate))+
  facet_grid(cols=vars(growth_michaelis_menten_constant))+
  geom_line()+
  scale_y_log10("num. diffusion steps")+
  scale_color_viridis_c()
p
```

![](README_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

``` r
p <- ggplot(df,aes(x=normal_oxygen_consumption_rate,
                   y=growth_michaelis_menten_constant,
                   fill=O2))+
  geom_raster()+
  scale_fill_viridis_c()
p
```

![](README_files/figure-gfm/unnamed-chunk-1-3.png)<!-- -->

``` r
path2sweepdata <- "output/test_homeostasis/data/"
runs <- list.files(path2sweepdata)

pars <- c("normal_oxygen_consumption_rate","growth_michaelis_menten_constant")

df <- do.call(rbind,lapply(runs,function(r){
  x <- read.csv(paste0(path2sweepdata,r,"/summary.csv"))
  parfile <- readLines(paste0(path2sweepdata,r,"/params.txt"))
  for(parName in pars){
    i <- grepl(parName,parfile)
    par <- parfile[i]
    par <- unlist(strsplit(par,split=" "))
    x$parval <- as.numeric(par[2])
    colnames(x)[colnames(x)=="parval"] <- par[1]
  }
  tail(x,1)
}))

p <- ggplot(df,aes(x=normal_oxygen_consumption_rate,
                   y=growth_michaelis_menten_constant,
                   fill=nNormal))+
  geom_raster()+
  scale_fill_viridis_c()
p
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
path2sweepdata <- "output/test_homeostasis/data/"
runs <- list.files(path2sweepdata)

pars <- c("normal_oxygen_consumption_rate","growth_michaelis_menten_constant")

df <- do.call(rbind,lapply(runs,function(r){
  x <- read.csv(paste0(path2sweepdata,r,"/summary.csv"))
  parfile <- readLines(paste0(path2sweepdata,r,"/params.txt"))
  for(parName in pars){
    i <- grepl(parName,parfile)
    par <- parfile[i]
    par <- unlist(strsplit(par,split=" "))
    x$parval <- as.numeric(par[2])
    colnames(x)[colnames(x)=="parval"] <- par[1]
  }
  x
}))

p <- ggplot(df,aes(x=time,y=nNormal,color=normal_oxygen_consumption_rate,
                   group=normal_oxygen_consumption_rate))+
  facet_grid(cols=vars(growth_michaelis_menten_constant))+
  geom_line()+
  scale_color_viridis_c()
p
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
