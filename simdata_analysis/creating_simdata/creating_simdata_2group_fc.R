#-------------------------------------------------------------------------------
# Generating simulated data for comparison between 2 groups (random fold change)
#-------------------------------------------------------------------------------

# Loading scripts
source("../../param.R")
source("../../package.R")

# Setting parameters
G <- param$g
PDEG <- param$pdeg[[1]]
assign <- param$assign$`2group`[[1]]
P1 <- assign[1]
P2 <- assign[2]
FC <- param$fc
n1 <- n2 <- param$rep[[1]]
data.name <- param$simdata.name
name <- data.name[1]


# Generating simulated data
## Preparation
simdata <- list()

## main loop
### TCC data
N_trial <- 100
for (i in 1:N_trial){
  print(i)
  set.seed(i)
  fc.matrix <- makeFCMatrix(Ngene = G,
                            PDEG = PDEG,
                            DEG.assign = c(P1, 1 - P1),
                            replicates = c(n1, n2))
  tcc <- simulateReadCounts(Ngene = G,
                            PDEG = PDEG,
                            DEG.assign = c(P1,1 - P1),
                            DEG.foldchange = c(FC, FC),
                            fc.matrix = fc.matrix,
                            replicates = c(n1, n2))
  simdata[[name]][[i]] <- tcc
}

# Output
out_f <- paste("../simdata_2group_fc_", name, "_", PDEG, "_", P1, "_fc_n", n1, ".obj", sep = "")
saveRDS(simdata, out_f)
