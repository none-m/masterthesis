#-------------------------------------------------------------------------------
# Generating simulated data for comparison between 3 groups (fixed fold change)
#-------------------------------------------------------------------------------

# Loading scripts
source("../../param.R")
source("../../package.R")


# Setting parameters
## fixed setting
G <- param$g
FC <- param$fc
data.name <- param$simdata.name
name <- data.name[1]

## manual setting
PDEG <- param$pdeg[[1]]
assign <- param$assign$`3group`[[1]]
P1 <- assign[1]
P2 <- assign[2]
P3 <- assign[3]
n1 <- n2 <- n3 <- param$rep[[1]]


# Generating simulated data
## preparation
simdata <- list()

## main loop
### TCC data
N_trial <- 50
for (i in 1:N_trial){
  print(i)
  set.seed(i)
  tcc<-simulateReadCounts(Ngene=G,
                          PDEG=PDEG,
                          DEG.assign=c(P1, P2, P3),
                          replicates=c(n1, n2, n3),
                          DEG.foldchange=c(FC,FC,FC))
  simdata[[name]][[i]] <- tcc
}


# Output
out_f <- paste("../simdata_3group_",name,"_",PDEG,"_",P1,"_fixed_n",n1,".obj",sep="")
saveRDS(simdata,out_f)
