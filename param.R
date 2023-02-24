#--------------------
# Setting parameters
#--------------------

param <- list()
param$g <- 10000
param$pdeg <- c(0.05,0.25,0.45,0.55,0.65,0.75)
param$assign$"2group" <- list(c(0.5,0.5),c(0.7,0.3),c(0.9,0.1),c(1.0,0))
param$assign$"3group" <- list(c(1/3,1/3,1/3),c(0.6,0.2,0.2),c(0.5,0.5,0),c(1.0,0,0))
param$fc <- 4
param$rep <- c(3,6,9,12)
param$tri <- c(50,100)
param$q <- 0.05
param$simdata.name <- c("tcc","comp","pro")
param$realdata.name <- c("liver","immune","skin","lym")
param$k <- list(3,2:4,2:5,c(2:5,seq(10,20,5)))
param$n.subsample <- c(2:10,15,20,25,30,35,40)
             
