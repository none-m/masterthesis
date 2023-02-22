#----------------------------------
# Obtaining genes detected as DEGs
#----------------------------------

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


# Obtaining DEGs'name and number of DEGs
## preparation
deg.list <- list()
ndeg.list <- list()

## main.loop
for (i in 1:length(res.all$edgeR)){
  ### edgeR
  deg.list$edgeR[[i]] <- rownames(res.all$edgeR[[i]]$table)[res.all$edgeR[[i]]$table$estimatedDEG == 1]
  ndeg.list$edgeR <- rbind(ndeg.list$edgeR, length(deg.list$edgeR[[i]]))
  
  ### DESeq2
  deg.list$DESeq2[[i]] <- rownames(res.all$DESeq2[[i]]$table)[res.all$DESeq2[[i]]$table$estimatedDEG == 1]
  ndeg.list$DESeq2 <- rbind(ndeg.list$DESeq2, length(deg.list$DESeq2[[i]]))
  
  ### TCC
  deg.list$TCC[[i]] <- res.all$TCC[[i]]$table$gene_id[res.all$TCC[[i]]$table$estimatedDEG == 1]
  ndeg.list$TCC <- rbind(ndeg.list$TCC, length(deg.list$TCC[[i]]))
  
  ### MBCdeg
  for(k in 1){
    ### MBCdeg1
    deg.list$MBCdeg1[[k]][[i]] <- rownames(res.all$MBCdeg1[[k]][[i]]$table)[res.all$MBCdeg1[[k]][[i]]$table$estimatedDEG == 1]
    ndeg.list$MBCdeg1[[k]] <- rbind(ndeg.list$MBCdeg1[[k]], length(deg.list$MBCdeg1[[k]][[i]]))
    
    ### MBCdeg2
    deg.list$MBCdeg2[[k]][[i]] <- rownames(res.all$MBCdeg2[[k]][[i]]$table)[res.all$MBCdeg2[[k]][[i]]$table$estimatedDEG == 1]
    ndeg.list$MBCdeg2[[k]] <- rbind(ndeg.list$MBCdeg2[[k]], length(deg.list$MBCdeg2[[k]][[i]]))
    
    ### MBCdeg3
    deg.list$MBCdeg3[[k]][[i]] <- rownames(res.all$MBCdeg3[[k]][[i]]$table)[res.all$MBCdeg3[[k]][[i]]$table$estimatedDEG == 1]
    ndeg.list$MBCdeg3[[k]] <- rbind(ndeg.list$MBCdeg3[[k]], length(deg.list$MBCdeg3[[k]][[i]]))
  }
}
names(deg.list$MBCdeg1) <- names(deg.list$MBCdeg2) <- names(deg.list$MBCdeg3) <- paste("K", k.list, sep = "")
names(ndeg.list$MBCdeg1) <- names(ndeg.list$MBCdeg2) <- names(ndeg.list$MBCdeg3) <- paste("K", k.list, sep = "")

# concatenating the results
matome <- data.frame(edgeR = ndeg.list$edgeR,
                     DESeq2 = ndeg.list$DESeq2,
                     TCC = ndeg.list$TCC,
                     MBCdeg1_K2 = ndeg.list$MBCdeg1$K2,
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
out_f <- paste("../res_deg_", name, "_ori.obj", sep = "")
saveRDS(deg.list, out_f)

# Output(number of DEGs)
out_f <- paste("../res_ndeg_", name, "_ori.txt", sep = "")
write.table(matome, out_f, sep = "\t", append = F, quote = F, row.names = F)
