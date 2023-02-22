#------------------------------
# DE analysis of permuted data
#------------------------------

# Loading scripts
source("../../param.R")
source("../../package.R")
source("../../pipeline.R")

# Loading data
in_f <- "../realdataset.obj"
data.all <- readRDS(in_f)

# Setting parameters
## fixed setting
data.name <- param$realdata.name
k.list <- param$k[[3]]
q <- param$q

## manual setting
name <- data.name[1]
n1 <- data.all[[name]]$rep[1]
n2 <- data.all[[name]]$rep[2]


# DE analysis
## preparation
res.list <- list()
res.list$MBCdeg1 <- res.list$MBCdeg2 <- res.list$MBCdeg3 <- vector("list", length = length(k.list))
data <- data.all[[name]]$data

## Generating of shuffle data labels
shuffledLabel = matrix(nrow = 100, ncol = ncol(data))
set.seed(2022)

for(i in 1:100){               
  x = rep(NA, ncol(data))                      
  x[sample(1:ncol(data), size = n1)] = "G1"      
  x[is.na(x)] = "G2"　　　　                  
  shuffledLabel[i,] = x                       
}
head(shuffledLabel)


## main loop
for (i in 1:nrow(shuffledLabel)){
  print(i)
  data.cl <- shuffledLabel[i,]
  
  ### edgeR
  res <- my.edger(data, data.cl)
  res$table$estimatedDEG <- rep(0, nrow(res))
  res$table$estimatedDEG[res$FDR <= q] <- 1
  res.list$edgeR[[i]] <- res
  
  ### DESeq2
  res <- my.deseq2(data, data.cl)
  res$table$estimatedDEG <- rep(0, nrow(res))
  res$table$estimatedDEG[res$padj <= q] <- 1
  res.list$DESeq2[[i]] <- res
  
  ### TCC
  res <- my.tcc(data, data.cl, q = q)
  res.list$TCC[[i]] <- res
  
  ### getting normfactors
  norm.factors <- res$norm.factors
  ef.libsizes <- colSums(data)*norm.factors
  size.factors <- ef.libsizes/mean(ef.libsizes)
  nf <- 1000000/colSums(data)
  
  ### MBCdeg
  for(k in 1:length(k.list)){
    res.list$MBCdeg1[[k]][[i]] <- my.mbcdeg(data, data.cl, k = k.list[k], q = q)
    res.list$MBCdeg2[[k]][[i]] <- my.mbcdeg(data, data.cl, normalizer = log(size.factors), k = k.list[k], q = q)
    res.list$MBCdeg3[[k]][[i]] <- my.mbcdeg(data, data.cl, normalizer = log(nf), k = k.list[k], q = q)
  }
}
names(res.list$MBCdeg1) <- names(res.list$MBCdeg2) <- names(res.list$MBCdeg3) <- paste("K", k.list, sep = "")

# Output
out_f <- paste("../res_", name, "_per.obj", sep = "")
saveRDS(res.list, out_f)
