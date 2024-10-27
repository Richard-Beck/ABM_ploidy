
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Error: provide path to config file")
}
path2config <- args[1]
config <- readLines(path2config)

config <- do.call(rbind,lapply(config,function(ci) unlist(strsplit(ci,split="="))))
colnames(config) <- c("argname","argval")

required <- c("path2BaseParFile","path2Output","sweeptype","javaCMD")

if(sum(required%in%config[,"argname"])<length(required)){
  stop(paste("config file requires at least: \n",paste(required,collapse=", ")))
}

path2Output <- config[config[,"argname"]=="path2Output","argval"]
dir.create(path2Output)
path2Pars <- paste0(path2Output,"/params/")
path2Sims <- paste0(path2Output,"/data/")
dir.create(path2Pars,recursive=T)
dir.create(path2Sims,recursive=T)


javaCMD <- config[config[,"argname"]=="javaCMD","argval"]

sweepMode <- config[config[,"argname"]=="sweeptype","argval"]
if(!sweepMode%in%c(".","x")) stop("sweepMode is either . or x")

basePars <- readLines(config[config[,"argname"]=="path2BaseParFile","argval"])
basePars <- do.call(rbind,lapply(basePars,function(bi) unlist(strsplit(bi,split=" "))))
colnames(basePars) <- c("parName","parVal")

pars2Sweep <- config[!config[,"argname"]%in%required,"argval"]
names(pars2Sweep) <- config[!config[,"argname"]%in%required,"argname"]

unpack_parvals <- function(parsStr){
  if(grepl(",",parsStr,fixed=T)){
    return(unlist(strsplit(parsStr,split=",")))
  } else if(grepl(":",parsStr,fixed=T)){
    pars <- unlist(strsplit(parsStr,split=":"))
    pars <- as.numeric(pars)
    pars <- seq(pars[1],pars[2],pars[3])
  } else {
    stop("Detected parameter with single value. Set this directly in BaseParFile instead")
  }
}

pars2Sweep <- lapply(pars2Sweep, unpack_parvals)
parPaths <- NULL
if(sweepMode=="x") {
  pars2Sweep <- expand.grid(pars2Sweep)
  stop("sweep mode x is not yet implemented.")  
}
if(sweepMode=="."){
  if(sum(basePars[,"parName"]%in%names(pars2Sweep))<length(pars2Sweep)){
    stop("default values must be provided in BaseParFile if sweepMode=. ")
  }
  parPaths <- unlist(sapply(names(pars2Sweep),function(ppi){
    tmp <- basePars
    i <- which(tmp[,"parName"]==ppi)
    sapply(pars2Sweep[[ppi]],function(j) {
      tmp[i,"parVal"] <- j
      tmp <- paste(tmp[,1],tmp[,2])
      saveVal <- gsub("[.]","p",j)
      saveName <- paste(ppi,saveVal,sep="_")
      savePath <- paste0(path2Pars,saveName,".txt")
      writeLines(tmp,savePath)
      return(savePath)
    })
  }))
}


cmds <- paste(javaCMD, parPaths, path2Sims)
for(cmd in cmds) system(cmd)

