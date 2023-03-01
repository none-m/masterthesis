#-------------------------------------------------------------------------------------
# Analysis of TCC simulation data for comparison between 3 groups (fixed fold change)
#-------------------------------------------------------------------------------------

# Loading scripts
source("../../param.R")
source("../../package.R")
source("../../pipeline.R")


# Setting parameters
## fixed setting
G <- param$g
q <- param$q
data.name <- param$simdata.name
name <- data.name[1]
k.list <- param$k[[3]]

## manual setting
PDEG <- param$pdeg[[1]]
assign <- param$assign$`3group`[[1]]
P1 <- assign[1]
P2 <- assign[2]
P3 <- assign[3]
n1 <- n2 <- n3 <- param$rep[[1]]


# Loading data
in_f <- paste("../simdata_3group_", name, "_", PDEG, "_", P1, "_fixed_n", n1, ".obj", sep = "")
data.all <- readRDS(in_f)


# DE analysis
## preparation
res.list <- list()
res.list$MBCdeg1 <- res.list$MBCdeg2 <- res.list$MBCdeg3 <- vector("list", length = length(k.list))

## main loop
N_trial <- 100
for (i in 1:N_trial){
  print(i)
  data <- data.all[[name]][[i]]$count
  data.cl <- data.all[[name]][[i]]$group$group
  obj <- as.numeric(data.all[[name]][[i]]$simulation$trueDEG != 0) 
  
  ### edgeR
  res.list$edgeR[[i]] <- my.edger(data, data.cl)
  
  ### DESeq2
  res.list$DESeq2[[i]] <- my.deseq2(data, data.cl)
  
  ### TCC
  res.list$TCC[[i]] <- my.tcc(data, data.cl, q = q)
  
  ### getting norm.factors
  norm.factors <- res.list$TCC[[i]]$norm.factors
  ef.libsizes <- colSums(data)*norm.factors
  size.factors <- ef.libsizes/mean(ef.libsizes)
  nf <- 1000000/colSums(data)
  
  ### MBCdeg
  for(k in 1:length(k.list)){
    res.list$MBCdeg1[[k]][[i]] <- my.mbcdeg(data, data.cl, k = k.list[k], q = q)
    res.list$MBCdeg2[[k]][[i]] <- my.mbcdeg(data, data.cl, normalizer = log(size.factors), k = k.list[k], q = q)
    res.list$MBCdeg3[[k]][[i]] <- my.mbcdeg(data, data.cl, normalizer = log(nf), k = k.list[k], q = q)
  }
  res.list$obj[[i]] <- obj
}
names(res.list$MBCdeg1) <- names(res.list$MBCdeg2) <- names(res.list$MBCdeg3) <- paste("K",k.list,sep="")


# Output
out_f <- paste("../res_3group_", name, "_", PDEG, "_", round(P1, digits = 2), "_fixed_n", n1, ".obj", sep = "")
saveRDS(res.list, out_f)
