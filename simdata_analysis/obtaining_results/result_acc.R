#--------------------------------------------------------------------------
# Obtaining DE analysis results for comparison between 2 groups (Accuracy)
#--------------------------------------------------------------------------

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
nDEG <- G*PDEG
nDEG1 <- nDEG*P1
n1 <- n2 <- param$rep[[1]]
name <- data.name[1]

# Loading data
in_f <- paste("../res_2group_",name,"_",PDEG,"_",P1,"_fixed_n",n1,".obj",sep="")
res.all <- readRDS(in_f)


# Obtaining Accuracy
## preparation
accuracy.list <- list()
matome <- NULL

## main loop
N_trial <- length(res.all$edgeR)
for (i in 1:N_trial){
  ### MBCdeg1
  res <- res.all$MBCdeg1$K3[[i]]$table
  centers <- res.all$MBCdeg1$K3[[i]]$centers
  k_cls <- nrow(res.all$MBCdeg1$K3[[i]]$centers)
  k_nonDEG <- res.all$MBCdeg1$K3[[i]]$nonDEGcluster
  
  #### Pattern of each gene
  pattern <- NULL
  pattern[1:k_cls] <- "DEG2"
  pattern[centers[,1] > 0] <- "DEG1"
  pattern[k_nonDEG] <- "nonDEG"
  pattern.list <- res$cluster
  for(i in 1:k_cls){
    pattern.list[pattern.list == i] <- pattern[i]
  }
  res$pattern <- pattern.list
  
  #### Identifying the k value that corresponds to the DEG cluster
  cls_DEG1<-which(pattern=='DEG1') 
  cls_DEG2<-which(pattern=='DEG2')
  cls_nonDEG<-which(pattern=='nonDEG')
  
  #### Checking the result of clustering
  hoge <- res$cluster[1:nDEG1]         # DEG1
  trueDEG1 <- sum(hoge == cls_DEG1)
  hoge <- res$cluster[(nDEG1+1):nDEG]  # DEG2
  trueDEG2 <- sum(hoge==cls_DEG2)
  hoge <- res$cluster[(nDEG+1):G]  # nonDEG
  truenonDEG<-sum(hoge==cls_nonDEG)
  true <- trueDEG1+trueDEG2+truenonDEG
  accuracy <- true/G
  accuracy.list$MBCdeg1 <- rbind(accuracy.list$MBCdeg1,accuracy) 
  
  
  ### MBCdeg2
  res <- res.all$MBCdeg2$K3[[i]]$table
  centers <- res.all$MBCdeg2$K3[[i]]$centers
  k_cls <- nrow(res.all$MBCdeg2$K3[[i]]$centers)
  k_nonDEG <- res.all$MBCdeg2$K3[[i]]$nonDEGcluster
  
  #### Pattern of each gene
  pattern <- NULL
  pattern[1:k_cls] <- "DEG2"
  pattern[centers[,1] > 0] <- "DEG1"
  pattern[k_nonDEG] <- "nonDEG"
  pattern.list <- res$cluster
  for(i in 1:k_cls){
    pattern.list[pattern.list == i] <- pattern[i]
  }
  res$pattern <- pattern.list
  
  #### Identifying the k value that corresponds to the DEG cluster
  cls_DEG1<-which(pattern=='DEG1') 
  cls_DEG2<-which(pattern=='DEG2')
  cls_nonDEG<-which(pattern=='nonDEG')
  
  #### Checking the result of clustering
  hoge <- res$cluster[1:nDEG1]         # DEG1
  trueDEG1 <- sum(hoge == cls_DEG1)
  hoge <- res$cluster[(nDEG1+1):nDEG]  # DEG2
  trueDEG2 <- sum(hoge==cls_DEG2)
  hoge <- res$cluster[(nDEG+1):G]  # nonDEG
  truenonDEG<-sum(hoge==cls_nonDEG)
  true <- trueDEG1+trueDEG2+truenonDEG
  accuracy <- true/G
  accuracy.list$MBCdeg2 <- rbind(accuracy.list$MBCdeg2,accuracy)
  
  
  ### MBCdeg3
  res <- res.all$MBCdeg3$K3[[i]]$table
  centers <- res.all$MBCdeg3$K3[[i]]$centers
  k_cls <- nrow(res.all$MBCdeg3$K3[[i]]$centers)
  k_nonDEG <- res.all$MBCdeg3$K3[[i]]$nonDEGcluster
  
  #### Pattern of each gene
  pattern <- NULL
  pattern[1:k_cls] <- "DEG2"
  pattern[centers[,1] > 0] <- "DEG1"
  pattern[k_nonDEG] <- "nonDEG"
  pattern.list <- res$cluster
  for(i in 1:k_cls){
    pattern.list[pattern.list == i] <- pattern[i]
  }
  res$pattern <- pattern.list
  
  #### Identifying the k value that corresponds to the DEG cluster
  cls_DEG1<-which(pattern=='DEG1') 
  cls_DEG2<-which(pattern=='DEG2')
  cls_nonDEG<-which(pattern=='nonDEG')
  
  #### Checking the result of clustering
  hoge <- res$cluster[1:nDEG1]         # DEG1
  trueDEG1 <- sum(hoge == cls_DEG1)
  hoge <- res$cluster[(nDEG1+1):nDEG]  # DEG2
  trueDEG2 <- sum(hoge==cls_DEG2)
  hoge <- res$cluster[(nDEG+1):G]  # nonDEG
  truenonDEG<-sum(hoge==cls_nonDEG)
  true <- trueDEG1+trueDEG2+truenonDEG
  accuracy <- true/G
  accuracy.list$MBCdeg3 <- rbind(accuracy.list$MBCdeg3,accuracy)
}
# concatenating the results
matome <- data.frame(MBCdeg1 = accuracy.list$MBCdeg1,
                     MBCdeg2 = accuracy.list$MBCdeg2,
                     MBCdeg3 = accuracy.list$MBCdeg3)

# Output
out_f <- paste("../acc_2group_",name,"_",PDEG,"_",P1,"_fixed_n",n1,".txt",sep="")
write.table(matome,out_f,sep="\t",append=F,quote=F,row.names=F)
