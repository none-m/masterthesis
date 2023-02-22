#--------------------------------------------------------------------------------
# Obtaining the average of 10 computation times for DE analysis of original data
#--------------------------------------------------------------------------------

# Loading scripts
source("../../param.R")
source("../../package.R")

# Setting parameters
## fixed setting
data.name <- param$realdata.name
k.list <- param$k[[3]]


# Obtaining computation times
df_mean <- NULL
for (name in data.name){
  # Loading data
  in_f <- paste("../res_", name, "_ori.obj", sep = "")
  res.all <- readRDS(in_f)
  
  data <- list(edgeR = res.all$edgeR,
               DESeq2 = res.all$DESeq2,
               TCC = res.all$TCC,
               MBCdeg1_K2 = res.all$MBCdeg1$K2,
               MBCdeg1_K3 = res.all$MBCdeg1$K3,
               MBCdeg1_K4 = res.all$MBCdeg1$K4,
               MBCdeg1_K5 = res.all$MBCdeg1$K5,
               MBCdeg2_K2 = res.all$MBCdeg2$K2,
               MBCdeg2_K3 = res.all$MBCdeg2$K3,
               MBCdeg2_K4 = res.all$MBCdeg2$K4,
               MBCdeg2_K5 = res.all$MBCdeg2$K5,
               MBCdeg3_K2 = res.all$MBCdeg3$K2,
               MBCdeg3_K3 = res.all$MBCdeg3$K3,
               MBCdeg3_K4 = res.all$MBCdeg3$K4,
               MBCdeg3_K5 = res.all$MBCdeg3$K5)
  
  # Calculation time to data frame
  df_ptm <- NULL
  for (j in 1:length(data)){
    list_ptm <- NULL
    for (i in 1:length(data$edgeR)){
      ptm <- data[[j]][[i]]$time
      list_ptm <- rbind(list_ptm, ptm)
    }
    df_ptm <- cbind(df_ptm, list_ptm)
  }
  
  colnames(df_ptm) <- names(data)
  
  # Calculating the average computation time
  mean <- colMeans(df_ptm)
  
  # Summarizing average computation time into data frames
  df_mean <- rbind(df_mean,mean)
}
rownames(df_mean) <- data.name

# Output
out_f <- "ptm_realdata_ori.txt"
write.table(df_mean, out_f, sep = "\t", append = F, quote = F, row.names = F)
