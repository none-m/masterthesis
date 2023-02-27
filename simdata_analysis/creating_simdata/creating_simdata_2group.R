#---------------------------------------------------------------------------------
# Generating simulated data for comparison between 2 groups (fixed fold change)
#---------------------------------------------------------------------------------

# Loading scripts
source("../../param.R")
source("../../package.R")


# Setting parameters
## fixed setting
G <- param$g
FC <- param$fc
data.name <- param$simdata.name

## manual setting
PDEG <- param$pdeg[[1]]
assign <- param$assign$`2group`[[1]]
P1 <- assign[1]
P2 <- assign[2]
n1 <- n2 <- param$rep[[1]]


# Generating simulated data
## preparation
simdata <- list()

## main loop
### TCC data
N_trial <- 100
for (i in 1:N_trial){
  print(i)
  set.seed(i)
  tcc <- simulateReadCounts(Ngene = G,
                            PDEG = PDEG,
                            DEG.assign= assign,
                            DEG.foldchange = c(FC, FC),
                            replicates = c(n1, n2))
  simdata[[data.name[1]]][[i]] <- tcc
}

### compcodeR data
N_trial <- 50
for (i in 1:N_trial){
  print(i)
  set.seed(i)
  compData <- generateSyntheticData(dataset = "compData", 
                                    n.vars = G,      
                                    samples.per.cond = n1, 
                                    n.diffexp = PDEG*G,
                                    fraction.upregulated = P2,
                                    filter.threshold.total = 1,
                                    filter.threshold.mediancpm = 0)
  simdata[[data.name[2]]][[i]] <- compData
}

### PROPER data
N_trial <- 50
for (i in 1:N_trial){
  print(i)
  set.seed(i)
  simOpts <- RNAseq.SimOptions.2grp(ngenes = G,
                                    p.DE = PDEG,
                                    lfc = c(rep(log(FC), G*PDEG*P1), rep(-log(FC), G*PDEG*P2)),
                                    lBaselineExpr = "bottomly",
                                    lOD = "bottomly",
                                    sim.seed = i)
  simres <- simRNAseq(simOpts, n1, n2)
  colnames(simres$counts) <- c(paste('G1_rep', 1:n1, sep = ''), paste('G2_rep', 1:n2, sep = ''))
  rownames(simres$counts) <- paste("Gene_", 1:G, sep = "")
  simdata[[data.name[3]]][[i]] <- simres
}


# Output
out_f <- paste("../simdata_2group_", PDEG, "_", P1, "_fixed_n", n1, ".obj", sep = "")
saveRDS(simdata, out_f)
