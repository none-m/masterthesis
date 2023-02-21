#---------------------------------------------------------------------
# Obtaining DE analysis results for comparison between 2 groups (FDR)
#---------------------------------------------------------------------

# Loading scripts
source("../../param.R")
source("../../package.R")

# Setting parameters
## fixed setting
G <- param$g
nominalFDR <- seq(0,1,0.001)
data.name <- param$simdata.name

## manual setting
PDEG <- param$pdeg[[1]]
nDEG <- G*PDEG
assign <- param$assign$`2group`[[1]]
P1 <- assign[1]
P2 <- assign[2]
n1 <- n2 <- param$rep[[1]]
name <- data.name[1]

# Loading data
in_f <- paste("../res_2group_",name,"_",PDEG,"_",P1,"_fixed_n",n1,".obj",sep="")
res.all <- readRDS(in_f)


# Obtaining FDR
## preparation
trueFDR.list <- list()
matome <- NULL

## main loop
N_trial <- length(res.all$edgeR)
for (i in 1:N_trial){
  ### edgeR
  res <- res.all$edgeR[[i]]$table
  res$FDR[is.na(res$FDR)] <- 1 
  trueFDR <- NULL
  for (t in nominalFDR){
    estimatedDEG <- sum(res$FDR<t) 
    falseDEG <- sum(res$FDR[(nDEG+1):G]<t) 
    trueFDR.i <- falseDEG/estimatedDEG  
    trueFDR <- rbind(trueFDR,trueFDR.i)
  }
  trueFDR[is.na(trueFDR)] <- 0 
  trueFDR.list$edgeR <- cbind(trueFDR.list$edgeR,trueFDR) 

    
  ### DESeq2
  res <- res.all$DESeq2[[i]]$table
  res$padj[is.na(res$padj)] <- 1
  trueFDR <- NULL
  for (t in nominalFDR){
    estimatedDEG <- sum(res$padj<t) 
    falseDEG <- sum(res$padj[(nDEG+1):G]<t) 
    trueFDR.i <- falseDEG/estimatedDEG  
    trueFDR <- rbind(trueFDR,trueFDR.i)
  }
  trueFDR[is.na(trueFDR)] <- 0 
  trueFDR.list$DESeq2 <- cbind(trueFDR.list$DESeq2,trueFDR) 

    
  ### MBCdeg1
  res <- res.all$MBCdeg1$K3[[i]]$table
  res.sort <- res[order(res$pp_nonDEG),] 
  cumsum <- cumsum(res.sort$pp_nonDEG) 
  qval <- cumsum/1:nrow(res) 
  res.sort<-cbind(res.sort,qval)
  
  #### qvalue
  trueFDR <- NULL
  for (t in nominalFDR){
    estimatedDEG<-sum(res.sort$qval<t) 
    falseDEG<-sum(as.numeric(gsub("gene_", "", row.names(res.sort)))[1:estimatedDEG]>nDEG) 
    trueFDR.i <- falseDEG/estimatedDEG  
    trueFDR <- rbind(trueFDR,trueFDR.i)
  }
  trueFDR[is.na(trueFDR)] <- 0 
  trueFDR.list$MBCdeg1$q <- cbind(trueFDR.list$MBCdeg1$q,trueFDR) 
  
  #### posterior probability
  trueFDR <- NULL
  for (t in nominalFDR){
    estimatedDEG<-sum(res.sort$pp_nonDEG<t) 
    if (estimatedDEG==0){ 
      falseDEG <- 0
    }else{
      falseDEG<-sum(as.numeric(gsub("gene_", "", row.names(res.sort)))[1:estimatedDEG]>nDEG)
    }
    trueFDR.i<-falseDEG/estimatedDEG 
    trueFDR <- rbind(trueFDR,trueFDR.i) 
  }
  trueFDR[is.na(trueFDR)] <- 0 
  trueFDR.list$MBCdeg1$pos <- cbind(trueFDR.list$MBCdeg1$pos,trueFDR) 
  
  
  ### MBCdeg2
  res <- res.all$MBCdeg2$K3[[i]]$table
  res.sort <- res[order(res$pp_nonDEG),] 
  cumsum <- cumsum(res.sort$pp_nonDEG) 
  qval <- cumsum/1:nrow(res)
  res.sort<-cbind(res.sort,qval)
  
  #### q value
  trueFDR <- NULL
  for (t in nominalFDR){
    estimatedDEG<-sum(res.sort$qval<t) 
    falseDEG<-sum(as.numeric(gsub("gene_", "", row.names(res.sort)))[1:estimatedDEG]>nDEG) 
    trueFDR.i<-falseDEG/estimatedDEG 
    trueFDR <- rbind(trueFDR,trueFDR.i) 
  }
  trueFDR[is.na(trueFDR)] <- 0 
  trueFDR.list$MBCdeg2$q <- cbind(trueFDR.list$MBCdeg2$q,trueFDR)
  
  
  #### posterior probability
  trueFDR <- NULL
  for (t in nominalFDR){
    estimatedDEG<-sum(res.sort$pp_nonDEG<t)
    if (estimatedDEG==0){ 
      falseDEG <- 0
    }else{
      falseDEG<-sum(as.numeric(gsub("gene_", "", row.names(res.sort)))[1:estimatedDEG]>nDEG)
    }
    trueFDR.i<-falseDEG/estimatedDEG 
    trueFDR <- rbind(trueFDR,trueFDR.i) 
  }
  trueFDR[is.na(trueFDR)] <- 0 
  trueFDR.list$MBCdeg2$pos <- cbind(trueFDR.list$MBCdeg2$pos,trueFDR) 
  
  
  ### MBCdeg3
  res <- res.all$MBCdeg3$K3[[i]]$table
  res.sort <- res[order(res$pp_nonDEG),] 
  cumsum <- cumsum(res.sort$pp_nonDEG) 
  qval <- cumsum/1:nrow(res) 
  res.sort<-cbind(res.sort,qval)
  
  #### q value
  trueFDR <- NULL
  for (t in nominalFDR){
    estimatedDEG<-sum(res.sort$qval<t) 
    falseDEG<-sum(as.numeric(gsub("gene_", "", row.names(res.sort)))[1:estimatedDEG]>nDEG) 
    trueFDR.i<-falseDEG/estimatedDEG 
    trueFDR  <- rbind(trueFDR,trueFDR.i) 
  }
  trueFDR[is.na(trueFDR)] <- 0 
  trueFDR.list$MBCdeg3$q <- cbind(trueFDR.list$MBCdeg3$q,trueFDR) 
  
  
  #### posterior probability
  trueFDR <- NULL
  for (t in nominalFDR){
    estimatedDEG<-sum(res.sort$pp_nonDEG<t) 
    if (estimatedDEG==0){ 
      falseDEG <- 0
    }else{
      falseDEG<-sum(as.numeric(gsub("gene_", "", row.names(res.sort)))[1:estimatedDEG]>nDEG)
    }
    trueFDR.i<-falseDEG/estimatedDEG 
    trueFDR <- rbind(trueFDR,trueFDR.i) 
  }
  trueFDR[is.na(trueFDR)] <- 0
  trueFDR.list$MBCdeg3$pos <- cbind(trueFDR.list$MBCdeg3$pos,trueFDR) 
}
# concatenating the results
matome <- data.frame(edgeR = rowMeans(trueFDR.list$edgeR),
                     DESeq2 = rowMeans(trueFDR.list$DESeq2),
                     MBCdeg1_q = rowMeans(trueFDR.list$MBCdeg1$q),
                     MBCdeg2_q = rowMeans(trueFDR.list$MBCdeg2$q),
                     MBCdeg3_q = rowMeans(trueFDR.list$MBCdeg3$q),
                     MBCdeg1_pos = rowMeans(trueFDR.list$MBCdeg1$pos),
                     MBCdeg2_pos = rowMeans(trueFDR.list$MBCdeg2$pos),
                     MBCdeg3_pos = rowMeans(trueFDR.list$MBCdeg3$pos))

# Output
out_f <- paste("../fdr_2group_",name,"_",PDEG,"_",P1,"_fixed_n",n1,".txt",sep="")
write.table(matome,out_f,sep="\t",append=F,quote=F,row.names=F)
