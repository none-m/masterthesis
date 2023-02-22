#---------------------------------------------------------
# DE analysis of original data (differential sample size)
#---------------------------------------------------------

# Loading scripts
source("../../param.R")
source("../../package.R")
source("../../pipeline.R")

# Loading data
in_f <- "../realdataset.obj"
data.all <- readRDS(in_f)
data <- data.all$immune$data


# Setting parameters
## fixed setting
n1 <- data.all$immune$rep[1]
n2 <- data.all$immune$rep[2]
k.list <- param$k[[3]]
data.cl <- c(rep(1, n1), rep(2, n2))
q <- param$q
n.subsample <- param$n.subsample


# down sampling
subdata.list <- list()
for (j in 1:length(n.subsample)){
  d.list <- list()
  for (i in 1:100){
    n <- n.subsample[j]
    d <- downsampling(data, n, n1, n2)
    d.list[[i]] <- d
  }
  subdata.list[[j]] <- d.list
}


# DE analysis
## main loop
for (j in 1:length(n.subsample)){
  res.list <- list()
  res.list[[j]] <- vector("list", length = 3)
  res.list[[j]]$MBCdeg1 <- res.list[[j]]$MBCdeg2 <- res.list[[j]]$MBCdeg3 <- vector("list", length = length(k.list))
  
  for (i in 1:2){
    print(i)
    ### Obtaining sub sample data
    subdata <- subdata.list[[j]][[i]]
    data.cl <- c(rep(1, n.subsample[j]), rep(2, n.subsample[j]))
    
    ### edgeR
    res <- my.edger(subdata, data.cl)
    res$table$estimatedDEG <- rep(0, nrow(res))
    res$table$estimatedDEG[res$FDR <= q] <- 1
    res.list[[j]]$edgeR[[i]] <- res

    ### DESeq2
    res <- my.deseq2(subdata, data.cl)
    res$table$estimatedDEG <- rep(0, nrow(res))
    res$table$estimatedDEG[res$padj <= q] <- 1
    res.list[[j]]$DESeq2[[i]] <- res
    
    ### TCC
    res <- my.tcc(subdata, data.cl, q = q)
    res.list[[j]]$TCC[[i]] <- res
    
    ### getting normfactors
    norm.factors <- res$norm.factors
    ef.libsizes <- colSums(subdata)*norm.factors
    size.factors <- ef.libsizes/mean(ef.libsizes)
    nf <- 1000000/colSums(subdata)
    
    ### MBCdeg
    for(k in 1:length(k.list)){
      res.list[[j]]$MBCdeg1[[k]][[i]] <- my.mbcdeg(subdata, data.cl, k = k.list[k], q = q)
      res.list[[j]]$MBCdeg2[[k]][[i]] <- my.mbcdeg(subdata, data.cl, normalizer = log(size.factors), k = k.list[k], q = q)
      res.list[[j]]$MBCdeg3[[k]][[i]] <- my.mbcdeg(subdata, data.cl, normalizer = log(nf), k = k.list[k], q = q)
    }
  }
  names(res.list[[j]]$MBCdeg1) <- names(res.list[[j]]$MBCdeg2) <- names(res.list[[j]]$MBCdeg3) <- paste("K", k.list, sep = "")
}
names(res.list) <- paste("n", n.subsample, sep = "")

# Output
out_f <- "../res_immune_subsample.obj"
saveRDS(res.list, out_f)
