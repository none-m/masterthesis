#------------------------------------------------------------------------------------
# DE Analysis of simulation data for comparison between 2 groups (fixed fold change)
#------------------------------------------------------------------------------------

# Loading scripts
source("../../param.R")
source("../../package.R")
source("../../pipeline.R")


# Setting parameters
## fixed setting
G <- param$g
q <- param$q
data.name <- param$simdata.name

## manual setting
PDEG <- param$pdeg[[1]]
assign <- param$assign$`2group`[[1]]
P1 <- assign[1]
P2 <- assign[2]
n1 <- n2 <- param$rep[[1]]
k.list <- param$k[[3]]
name <- data.name[1]
N_trial <- 100


# Loading data
in_f <- paste("../simdata_2group_", PDEG, "_", P1, "_fixed_n", n1, ".obj", sep = "")
data.all <- readRDS(in_f)


# DE analysis
## preparation
res.list <- list()
res.list$MBCdeg1 <- res.list$MBCdeg2 <- res.list$MBCdeg3 <- vector("list", length = length(k.list))

## main loop
for (i in 1:N_trial){
  print(i)
  
  if (name == "tcc"){
    data <- data.all[[name]][[i]]$count
    data.cl <- data.all[[name]][[i]]$group$group
    obj <- as.numeric(data.all[[name]][[i]]$simulation$trueDEG != 0)
  }else if (name == "comp"){
    data <- data.all[[name]][[i]]@count.matrix
    data.cl <- data.all[[name]][[i]]@sample.annotations$condition
    obj <- data.all[[name]][[i]]@variable.annotations$differential.expression 
  }else if (name == "pro"){
    data <- data.all[[name]][[i]]$counts
    data.cl <- data.all[[name]][[i]]$designs
    obj <- rep(0,nrow(data))
    for (d in data.all[[name]][[i]]$DEid){
      obj[d] <- 1
    }
  }
  
  ### edgeR
  res.list$edgeR[[i]] <- my.edger(data, data.cl)
  
  ### DESeq2
  res.list$DESeq2[[i]] <- my.deseq2(data, data.cl)
  
  ### TCC
  res.list$TCC[[i]] <- my.tcc(data, data.cl, q = q)
  
  #### getting norm.factors
  norm.factors <- res.list$TCC[[i]]$norm.factors
  ef.libsizes <- colSums(data)*norm.factors
  size.factors <- ef.libsizes/mean(ef.libsizes)
  nf <- 1000000/colSums(data)
  
  ### MBCdeg
  for(k in 1:length(k.list)){
    res.list$MBCdeg1[[k]][[i]] <- my.mbcdeg(data, data.cl, k = k.list[k],q = q)
    res.list$MBCdeg2[[k]][[i]] <- my.mbcdeg(data, data.cl, normalizer = log(size.factors), k = k.list[k],q = q) 
    res.list$MBCdeg3[[k]][[i]] <- my.mbcdeg(data, data.cl, normalizer = log(nf), k = k.list[k], q = q)
  }
  res.list$obj[[i]] <- obj
}
names(res.list$MBCdeg1) <- names(res.list$MBCdeg2) <- names(res.list$MBCdeg3) <- paste("K", k.list, sep = "")


# Output
out_f <- paste("../res_2group_", name, "_", PDEG, "_", P1, "_fixed_n", n1, ".obj", sep = "")
saveRDS(res.list, out_f)
