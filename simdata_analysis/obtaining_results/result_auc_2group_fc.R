#-------------------------------------------------------------------------
# Obtaining DE analysis results for comparison between 2 groups (ROC-AUC)
#-------------------------------------------------------------------------

# Loading scripts
source("../../param.R")
source("../../package.R")

# Setting parameters
## fixed setting
G <- param$g
data.name <- param$simdata.name

## manual setting
PDEG <- param$pdeg[[1]]
assign <- param$assign$`2group`[[1]]
P1 <- assign[1]
P2 <- assign[2]
n1 <- n2 <- param$rep[[1]]
k.list <- param$k[[1]]
name <- data.name[1]

# Loading data
in_f <- paste("../res_2group_", name, "_", PDEG, "_", P1, "_fc_n", n1, ".obj", sep = "")
res.all <- readRDS(in_f)


# Obtaining AUC value
## preparation
auc.list <- list()
auc.list$MBCdeg1 <- auc.list$MBCdeg2 <- auc.list$MBCdeg3 <- vector("list", length = length(k.list))
matome <- NULL

## main loop
N_trial <- length(res.all$edgeR)
for (i in 1:N_trial){
  obj <- res.all$obj
  
  ### edgeR
  ranking <- res.all$edgeR[[i]]$table$rank
  auc <- AUC(rocdemo.sca(truth = obj, data = -ranking))
  auc.list$edgeR <- rbind(auc.list$edgeR, auc)
  
  ### DESeq2
  ranking <- res.all$DESeq2[[i]]$table$rank
  auc <- AUC(rocdemo.sca(truth = obj, data = -ranking))
  auc.list$DESeq2 <- rbind(auc.list$DESeq2, auc)
  
  ### TCC
  ranking <- res.all$TCC[[i]]$table$rank
  auc <- AUC(rocdemo.sca(truth = obj, data = -ranking))
  auc.list$TCC <- rbind(auc.list$TCC, auc)
  
  
  ### MBCdeg
  for (k in 1:length(k.list)){
    ### MBCdeg1
    ranking <- res.all$MBCdeg1[[k]][[i]]$table$rank
    auc <- AUC(rocdemo.sca(truth = obj, data = -ranking))
    auc.list$MBCdeg1[[k]] <- rbind(auc.list$MBCdeg1[[k]], auc)
    
    ### MBCdeg2
    ranking <- res.all$MBCdeg2[[k]][[i]]$table$rank
    auc <- AUC(rocdemo.sca(truth = obj, data = -ranking))
    auc.list$MBCdeg2[[k]] <- rbind(auc.list$MBCdeg2[[k]], auc)
    
    ### MBCdeg3
    ranking <- res.all$MBCdeg3[[k]][[i]]$table$rank
    auc <- AUC(rocdemo.sca(truth = obj, data = -ranking))
    auc.list$MBCdeg3[[k]] <- rbind(auc.list$MBCdeg3[[k]], auc)
  }
}
names(auc.list$MBCdeg1) <- names(auc.list$MBCdeg2) <- names(auc.list$MBCdeg3) <- paste("K", k.list, sep = "")

# concatenating the results
matome <- data.frame(edgeR = auc.list$edgeR,
                     DESeq2 = auc.list$DESeq2,
                     TCC = auc.list$TCC,
                     MBCdeg1_K2 = auc.list$MBCdeg1$K2,
                     MBCdeg1_K3 = auc.list$MBCdeg1$K3,
                     MBCdeg1_K4 = auc.list$MBCdeg1$K4,
                     MBCdeg1_K5 = auc.list$MBCdeg1$K5,
                     MBCdeg1_K10 = auc.list$MBCdeg1$K10,
                     MBCdeg1_K15 = auc.list$MBCdeg1$K15,
                     MBCdeg1_K20 = auc.list$MBCdeg1$K20,
                     MBCdeg2_K2 = auc.list$MBCdeg2$K2,
                     MBCdeg2_K3 = auc.list$MBCdeg2$K3,
                     MBCdeg2_K4 = auc.list$MBCdeg2$K4,
                     MBCdeg2_K5 = auc.list$MBCdeg2$K5,
                     MBCdeg2_K10 = auc.list$MBCdeg2$K10,
                     MBCdeg2_K15 = auc.list$MBCdeg2$K15,
                     MBCdeg2_K20 = auc.list$MBCdeg2$K20,
                     MBCdeg3_K2 = auc.list$MBCdeg3$K2,
                     MBCdeg3_K3 = auc.list$MBCdeg3$K3,
                     MBCdeg3_K4 = auc.list$MBCdeg3$K4,
                     MBCdeg3_K5 = auc.list$MBCdeg3$K5,
                     MBCdeg3_K10 = auc.list$MBCdeg3$K10,
                     MBCdeg3_K15 = auc.list$MBCdeg3$K15,
                     MBCdeg3_K20 = auc.list$MBCdeg3$K20)

# Output
out_f <- paste("../auc_2group_", name, "_", PDEG, "_", P1, "_fc_n", n1, ".txt", sep = "")
write.table(matome, out_f, sep = "\t", append = F, quote = F, row.names = F)
