#----------------------------
# Filtering of detected DEGs
#----------------------------

# Loading scripts
source("../../param.R")
source("../../package.R")

# Setting parameters
## fixed setting
data.name <- param$realdata.name
k.list <- param$k[[3]]

## manual setting
name <- data.name[1]


# Loading data
in_f <- paste("../res_", name, "_ori.obj", sep = "")
res.all <- readRDS(in_f)


# Filtering detected DEGs
## preparation
res.mbcdeg <- list(MBCdeg1 = list(K2 = res.all$MBCdeg1$K2[[1]]$table,
                                  K3 = res.all$MBCdeg1$K3[[1]]$table,
                                  K4 = res.all$MBCdeg1$K4[[1]]$table,
                                  K5 = res.all$MBCdeg1$K5[[1]]$table),
                   MBCdeg2 = list(K2 = res.all$MBCdeg2$K2[[1]]$table,
                                  K3 = res.all$MBCdeg2$K3[[1]]$table,
                                  K4 = res.all$MBCdeg2$K4[[1]]$table,
                                  K5 = res.all$MBCdeg2$K5[[1]]$table),
                   MBCdeg3 = list(K2 = res.all$MBCdeg2$K2[[1]]$table,
                                  K3 = res.all$MBCdeg2$K3[[1]]$table,
                                  K4 = res.all$MBCdeg2$K4[[1]]$table,
                                  K5 = res.all$MBCdeg2$K5[[1]]$table))
                   

## Filtering
### MBCdeg1
for (k in 1:length(k.list)){
  for (g in 1:nrow(res.mbcdeg$MBCdeg1$K2)){
    if(max(res.mbcdeg$MBCdeg1[[k]][g, 1:3]) < 0.5){
      res.mbcdeg$MBCdeg1[[k]][g,]$estimatedDEG <- 2
    }
  }
}

### MBCdeg2
for (k in 1:length(k.list)){
  for (g in 1:nrow(res.mbcdeg$MBCdeg2$K2)){
    if(max(res.mbcdeg$MBCdeg2[[k]][g, 1:3]) < 0.5){
      res.mbcdeg$MBCdeg2[[k]][g,]$estimatedDEG <- 2
    }
  }
}

### MBCdeg3
for (k in 1:length(k.list)){
  for (g in 1:nrow(res.mbcdeg$MBCdeg3$K2)){
    if(max(res.mbcdeg$MBCdeg3[[k]][g, 1:3]) < 0.5){
      res.mbcdeg$MBCdeg3[[k]][g,]$estimatedDEG <- 2
    }
  }
}


# Obtaining DEGs'name and number of DEGs
## Preparation
deg.list <- list()
ndeg.list <- list()

## main.loop
for (k in 1:length(k.list)){
  ### MBCdeg1
  deg.list$MBCdeg1[[k]] <- rownames(res.mbcdeg$MBCdeg1[[k]])[res.mbcdeg$MBCdeg1[[k]]$estimatedDEG == 1]
  ndeg.list$MBCdeg1[[k]] <- rbind(ndeg.list$MBCdeg1[[k]], length(deg.list$MBCdeg1[[k]]))
  
  ### MBCdeg2
  deg.list$MBCdeg2[[k]] <- rownames(res.mbcdeg$MBCdeg2[[k]])[res.mbcdeg$MBCdeg2[[k]]$estimatedDEG == 1]
  ndeg.list$MBCdeg2[[k]] <- rbind(ndeg.list$MBCdeg2[[k]], length(deg.list$MBCdeg2[[k]]))
  
  ### MBCdeg3
  deg.list$MBCdeg3[[k]] <- rownames(res.mbcdeg$MBCdeg3[[k]])[res.mbcdeg$MBCdeg3[[k]]$estimatedDEG == 1]
  ndeg.list$MBCdeg3[[k]] <- rbind(ndeg.list$MBCdeg3[[k]], length(deg.list$MBCdeg3[[k]]))
}
names(deg.list$MBCdeg1) <- names(deg.list$MBCdeg2) <- names(deg.list$MBCdeg3) <- paste("K", k.list, sep = "")
names(ndeg.list$MBCdeg1) <- names(ndeg.list$MBCdeg2) <- names(ndeg.list$MBCdeg3) <- paste("K", k.list, sep = "")

# concatenating the results
matome <- data.frame(MBCdeg1_K2 = ndeg.list$MBCdeg1$K2,
                     MBCdeg1_K3 = ndeg.list$MBCdeg1$K3,
                     MBCdeg1_K4 = ndeg.list$MBCdeg1$K4,
                     MBCdeg1_K5 = ndeg.list$MBCdeg1$K5,
                     MBCdeg2_K2 = ndeg.list$MBCdeg2$K2,
                     MBCdeg2_K3 = ndeg.list$MBCdeg2$K3,
                     MBCdeg2_K4 = ndeg.list$MBCdeg2$K4,
                     MBCdeg2_K5 = ndeg.list$MBCdeg2$K5,
                     MBCdeg3_K2 = ndeg.list$MBCdeg3$K2,
                     MBCdeg3_K3 = ndeg.list$MBCdeg3$K3,
                     MBCdeg3_K4 = ndeg.list$MBCdeg3$K4,
                     MBCdeg3_K5 = ndeg.list$MBCdeg3$K5)
  
# Output(DEGs'name)
out_f <- paste("../res_deg_", name, "_ori_fil.obj", sep = "")
saveRDS(deg.list,out_f)

# Output(number of DEGs)
out_f <- paste("../res_ndeg_", name, "_ori_fil.txt", sep = "")
write.table(matome, out_f, sep = "\t", append = F, quote = F, row.names = F) 
