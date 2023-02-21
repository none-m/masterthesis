#----------------------------------------------------------------------------------------------------
# Obtaining the average of 50 computation times for comparison between 2 groups (fixed fold change)
#----------------------------------------------------------------------------------------------------

# Loading scripts
source("../../param.R")
source("../../package.R")

# Setting parameters
## fixed setting
G <- param$g
data.name <- param$simdata.name

## manual setting
PDEG <- param$pdeg[[1]]
n1 <- n2 <- param$rep[[1]]
k.list <- param$k[[3]]
name <- data.name[1]


# Obtaining computation times
df_mean <- NULL
for (assign in param$assign$`2group`){
  P1 <- assign[1]
  P2 <- assign[2]
  
  # Loading data
  in_f <- paste("../res_2group_",name,"_",PDEG,"_",P1,"_fixed_n",n1,".obj",sep="")
  data.all <- readRDS(in_f)
  
  data <- list(edgeR = data.all$edgeR,
                    DESeq2 = data.all$DESeq2,
                    TCC = data.all$TCC,
                    MBCdeg1_K2 = data.all$MBCdeg1$K2,
                    MBCdeg1_K3 = data.all$MBCdeg1$K3,
                    MBCdeg1_K4 = data.all$MBCdeg1$K4,
                    MBCdeg2_K2 = data.all$MBCdeg2$K2,
                    MBCdeg2_K3 = data.all$MBCdeg2$K3,
                    MBCdeg2_K4 = data.all$MBCdeg2$K4,
                    MBCdeg3_K2 = data.all$MBCdeg3$K2,
                    MBCdeg3_K3 = data.all$MBCdeg3$K3,
                    MBCdeg3_K4 = data.all$MBCdeg3$K4)
  
  # Calculation time to data frame
  df_ptm <- NULL
  for (j in 1:length(data)){
    list_ptm <- NULL
    for (i in 1:length(data$edgeR)){
      ptm <- data[[j]][[i]]$time
      list_ptm <- rbind(list_ptm,ptm)
    }
    df_ptm <- cbind(df_ptm,list_ptm)
  }
  
  colnames(df_ptm) <- names(data)
  
  # Calculating the average computation time
  mean <- colMeans(df_ptm)
  
  # Summarizing average computation time into data frames
  df_mean <- rbind(df_mean,mean)
}

rownames(df_mean) <- param$assign$`2group`

# Output
out_f <- paste("../ptm_2group_",name,"_",PDEG,"_fixed_n",n1,".txt",sep="")
write.table(df_mean,out_f,sep="\t",append=F,quote=F,row.names=T)

