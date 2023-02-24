#-------------------------------------
# Defining functions for data analysis
#-------------------------------------

# edgeR
my.edger <- function(counts,group,q=1){
  ptm <- proc.time()
  dge <- edgeR::DGEList(counts=counts,
                        group=group)
  dge <- edgeR::calcNormFactors(dge)               
  design <- model.matrix(~group)
  dge <- edgeR::estimateDisp(dge, design)       
  fit <- edgeR::glmQLFit(dge, design)           
  qlf <- edgeR::glmQLFTest(fit, coef = 2,)
  table <- topTags(qlf, n = nrow(counts), sort.by = "none",p.value = q)@.Data[[1]]
  ptm <- proc.time() - ptm
  table$rank <- rank(table$PValue)
  res <- list(table = table,
              time = ptm[3])
 
  return(res)                                      
}


# DESeq2
my.deseq2 <- function(counts,group){
  ptm <- proc.time()
  colData <- data.frame(condition=as.factor(group))
  dds <- DESeqDataSetFromMatrix(countData=counts,
                                colData=colData, design=~condition)
  dds <- DESeq(dds)
  tmp <- results(dds)
  ptm <- proc.time() - ptm
  tmp$pvalue[is.na(tmp$pvalue)] <- 1
  tmp$rank <- rank(tmp$pvalue)
  res <- list(table = tmp,
              time = ptm[3])
  
  return(res)
}


# TCC
my.tcc <- function(counts,group,q=1){
  ptm <- proc.time()
  tcc <- new('TCC',counts,group)
  tcc<-calcNormFactors(tcc,
                       norm.method="tmm",
                       test.method="edger", 
                       iteration=3, FDR=0.1,floorPDEG=0.05)
  tcc<-estimateDE(tcc, test.method="edger", FDR=q)　                
  df <- getResult(tcc, sort=FALSE)
  ptm <- proc.time() - ptm
  res <- list(table = df,
              norm.factors = tcc$norm.factors,
              time = ptm[3])
  
  return(res)                                                    
}


# MBCdeg
my.mbcdeg <- function(counts,group,normalizer = NULL,k=3,q=1){
  ptm <- proc.time()
  hoge<-RNASeq.Data(counts,
                    Normalizer = normalizer,
                    Treatment = group,
                    GeneID = rownames(counts))
  c0<-KmeansPlus.RNASeq(hoge,
                        nK=k,
                        model="nbinom", 
                        print.steps=F)
  cls <- Cluster.RNASeq(data=hoge, 
                        model="nbinom",
                        centers=c0$centers, 
                        method="EM")
  ptm <- proc.time() - ptm
  L2norm <- sqrt(rowSums(abs(cls$centers)^2))
  K <- which.min(L2norm)
  cluster<-cls$cluster
  pp <- cls$probability
  pp_nonDEG <- cls$probability[,K]
  centers <- cls$centers
  estimatedDEG <- rep(NA,length(cluster))
  estimatedDEG[cls$cluster!=K | cls$probability[,K]<q] <- 1
  estimatedDEG[is.na(estimatedDEG)] <- 0                       
  df <- data.frame(pp = pp,
                   pp_nonDEG = pp_nonDEG,
                   cluster = cluster,
                   rank = rank(pp_nonDEG),
                   estimatedDEG = estimatedDEG)
  rownames(df) <- rownames(counts)
  res <- list(table = df,
              centers = centers,
              nonDEGcluster = K,
              time = ptm[3])
  
  return(res)
}


# downsampling
downsampling <- function(count,n,n1,n2){  
  g1 <- count[,1:n1]
  g2 <- count[,(n1+1):(n1+n2)]
  select_g1 <- sample(colnames(g1),n)
  subsample_g1 <- g1[,select_g1]
  
  select_g2 <- sample(colnames(g2),n)
  subsample_g2 <- g2[,select_g2]
  subsumple <- cbind(subsample_g1,subsample_g2)
  
  return(subsumple)
}

