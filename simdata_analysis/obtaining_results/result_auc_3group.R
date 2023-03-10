#-------------------------------------------------------------------------
# Obtaining DE analysis results for comparison between 3 groups (ROC-AUC)
#-------------------------------------------------------------------------

# Loading scripts
source("../../param.R")
source("../../package.R")

# Setting parameters
## fixed setting
G <- param$g
n1 <- n2 <- n3 <- param$rep[[1]]
data.name <- param$simdata.name

## manual setting
PDEG <- param$pdeg[[1]]
assign <- param$assign$`3group`[[1]]
P1 <- assign[1]
P2 <- assign[2]
P3 <- assign[3]
k.list <- param$k[[3]]
name <- data.name[1]


# Loading data
in_f <- paste("../res_3group_", name, "_", PDEG, "_", round(P1, digits = 2), "_fixed_n", n1, ".obj", sep = "")
res.all <- readRDS(in_f)


# Obtaining AUC value
## preparation
auc.list <- list()
auc.list$MBCdeg1 <- auc.list$MBCdeg2 <- auc.list$MBCdeg3 <- vector("list", length = length(k.list))
matome <- NULL

## main loop
N_trial <- length(res.all$edgeR)
for (i in 1:N_trial){
  obj <- res.all$obj[[i]]
  
  ### edgeR
  ranking <- res.all$edgeR[[i]]$table$rank
  auc <- AUC(rocdemo.sca(truth=obj, data = -ranking))
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

# concatenate the results
matome <- data.frame(edgeR = auc.list$edgeR,
                     DESeq2 = auc.list$DESeq2,
                     TCC = auc.list$TCC,
                     MBCdeg1_K2 = auc.list$MBCdeg1$K2,
                     MBCdeg1_K3 = auc.list$MBCdeg1$K3,
                     MBCdeg1_K4 = auc.list$MBCdeg1$K4,
                     MBCdeg1_K5 = auc.list$MBCdeg1$K5,
                     MBCdeg2_K2 = auc.list$MBCdeg2$K2,
                     MBCdeg2_K3 = auc.list$MBCdeg2$K3,
                     MBCdeg2_K4 = auc.list$MBCdeg2$K4,
                     MBCdeg2_K5 = auc.list$MBCdeg2$K5,
                     MBCdeg3_K2 = auc.list$MBCdeg3$K2,
                     MBCdeg3_K3 = auc.list$MBCdeg3$K3,
                     MBCdeg3_K4 = auc.list$MBCdeg3$K4,
                     MBCdeg3_K5 = auc.list$MBCdeg3$K5)

# Output
out_f <- paste("../auc_3group_", name, "_", PDEG, "_", round(P1, digits = 2), "_fixed_n", n1, ".txt", sep = "")
write.table(matome, out_f, sep = "\t", append = F, quote = F, row.names = F)
