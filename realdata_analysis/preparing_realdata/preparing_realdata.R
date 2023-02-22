#----------------------------
# Preparating real data sets
#----------------------------

# Loading scripts
source("../../package.R")

# Preparating real data sets
## preparation
realdata <- list()

## Liver samples data set
n1 <- n2 <- 4
data.all <- read.table('sample_blekhman_18.txt') 
data <- data.frame(data.all$HSF1, 
                   data.all$HSF2,
                   data.all$HSM2,
                   data.all$HSM3,
                   data.all$RMF1,
                   data.all$RMF3,
                   data.all$RMM1,
                   data.all$RMM3)
obj <- as.logical(rowSums(data) > 0) 
data <- data[obj,]
colnames(data)<-c(paste('G1_rep', 1:n1, sep = ''), paste('G2_rep', 1:n2, sep = '')) 

realdata$liver$data <- data
realdata$liver$rep <- c(n1, n2)


## Immune therapy data set
n1 <- 51
n2 <- 58
data.all <- read.csv('GSE91061_BMS038109Sample.hg19KnownGene.raw.csv', row.names = 1)  
Col_pre <- str_detect(colnames(data.all), pattern = "Pre") 
Col_on <- str_detect(colnames(data.all), pattern = "On")　 
data_pre <- data.all[, Col_pre] 
data_on <- data.all[, Col_on]   
data <- cbind(data_pre, data_on)
obj <- as.logical(rowSums(data) > 0) 
data <- data[obj,]
colnames(data) <- c(paste('G1_rep', 1:n1, sep = ''), paste('G2_rep', 1:n2, sep = '')) 

realdata$immune$data <- data
realdata$immune$rep <- c(n1, n2)


## Skin biopsies data set
n1 <- n2 <- 74
param_ID <- "SRP035988" 
download_study(param_ID, type = "rse-gene", download = T) 
load(file.path(param_ID, 'rse_gene.Rdata')) 
rse <- rse_gene                        
rse <- scale_counts(rse)               
x <- assays(rse)$counts                
x <- as.data.frame(x)                  
data <- cbind(x$SRR1146077, x$SRR1146078, x$SRR1146079,#lesional01-03 
              x$SRR1146080, x$SRR1146081, x$SRR1146082,#lesional04-06
              x$SRR1146083 + x$SRR1146084, x$SRR1146087, x$SRR1146089,#lesional07-09
              x$SRR1146090, x$SRR1146091, x$SRR1146092,#lesional10-12
              x$SRR1146093, x$SRR1146094, x$SRR1146095,#lesional13-15
              x$SRR1146097, x$SRR1146098 + x$SRR1146099, x$SRR1146100,#lesional16-18
              x$SRR1146102, x$SRR1146103, x$SRR1146104,#lesional19-21
              x$SRR1146105, x$SRR1146110, x$SRR1146112,#lesional22-24
              x$SRR1146119, x$SRR1146123, x$SRR1146124,#lesional25-27
              x$SRR1146128, x$SRR1146129, x$SRR1146132,#lesional28-30
              x$SRR1146133, x$SRR1146134, x$SRR1146135,#lesional31-33
              x$SRR1146136, x$SRR1146140, x$SRR1146141,#lesional34-36
              x$SRR1146147, x$SRR1146148, x$SRR1146149,#lesional37-39
              x$SRR1146150, x$SRR1146151, x$SRR1146152,#lesional40-42
              x$SRR1146154, x$SRR1146155, x$SRR1146156,#lesional43-45
              x$SRR1146157, x$SRR1146158, x$SRR1146159,#lesional46-48
              x$SRR1146160, x$SRR1146161, x$SRR1146162,#lesional49-51
              x$SRR1146163, x$SRR1146164, x$SRR1146165,#lesional52-54
              x$SRR1146168, x$SRR1146199, x$SRR1146202,#lesional55-57
              x$SRR1146203, x$SRR1146204, x$SRR1146205,#lesional58-60
              x$SRR1146206, x$SRR1146208, x$SRR1146209,#lesional61-63
              x$SRR1146210, x$SRR1146211, x$SRR1146212,#lesional64-66
              x$SRR1146214, x$SRR1146215, x$SRR1146216 + x$SRR1146217,#lesional67-69
              x$SRR1146218, x$SRR1146219, x$SRR1146220,#lesional70-72
              x$SRR1146221, x$SRR1146222, x$SRR1146223,#lesional73-75
              x$SRR1146224, x$SRR1146225, x$SRR1146226,#lesional76-78
              x$SRR1146227, x$SRR1146228, x$SRR1146229,#lesional79-81
              x$SRR1146233, x$SRR1146234, x$SRR1146235,#lesional82-84
              x$SRR1146237, x$SRR1146239, x$SRR1146240,#lesional85-87
              x$SRR1146241, x$SRR1146242, x$SRR1146252,#lesional88-90
              x$SRR1146253, x$SRR1146254,              #lesional91-92
              x$SRR1146076, x$SRR1146085 + x$SRR1146086, x$SRR1146088,#normal01-03
              x$SRR1146096, x$SRR1146101, x$SRR1146106,#normal04-06
              x$SRR1146107, x$SRR1146108, x$SRR1146109,#normal07-09
              x$SRR1146111, x$SRR1146113, x$SRR1146114,#normal10-12
              x$SRR1146115, x$SRR1146116, x$SRR1146117,#normal13-15
              x$SRR1146118, x$SRR1146120, x$SRR1146121,#normal16-18
              x$SRR1146122, x$SRR1146125, x$SRR1146126,#normal19-21
              x$SRR1146127, x$SRR1146130, x$SRR1146131,#normal22-24
              x$SRR1146137, x$SRR1146138, x$SRR1146139,#normal25-27
              x$SRR1146142, x$SRR1146143, x$SRR1146144,#normal28-30
              x$SRR1146145, x$SRR1146146, x$SRR1146153,#normal31-33
              x$SRR1146166, x$SRR1146167, x$SRR1146169,#normal34-36
              x$SRR1146170, x$SRR1146171, x$SRR1146172,#normal37-39
              x$SRR1146173, x$SRR1146174, x$SRR1146175,#normal40-42
              x$SRR1146176, x$SRR1146177, x$SRR1146178,#normal43-45
              x$SRR1146179, x$SRR1146180, x$SRR1146181,#normal46-48
              x$SRR1146182, x$SRR1146183, x$SRR1146184,#normal49-51
              x$SRR1146185, x$SRR1146186, x$SRR1146187,#normal52-54
              x$SRR1146188, x$SRR1146189, x$SRR1146190,#normal55-57
              x$SRR1146191, x$SRR1146192, x$SRR1146193,#normal58-60
              x$SRR1146194, x$SRR1146195, x$SRR1146196,#normal61-63
              x$SRR1146197, x$SRR1146198, x$SRR1146200,#normal64-66
              x$SRR1146201, x$SRR1146207, x$SRR1146213,#normal67-69
              x$SRR1146230, x$SRR1146231 + x$SRR1146232, x$SRR1146236,#normal70-72
              x$SRR1146238, x$SRR1146243, x$SRR1146244,#normal73-75
              x$SRR1146245, x$SRR1146246, x$SRR1146247,#normal76-78
              x$SRR1146248, x$SRR1146249, x$SRR1146250)#normal79-81
obj <- as.logical(rowSums(data) > 0) 
data <- data[obj,]
data <- as.data.frame(data)
data <- select(.data = data, -c(2, 35, 6, 10, 29, 9, 19, 13, 56, 5, 8, 1, 12, 74, 31, 87, 88, 69,  
                                145, 160, 132, 165, 95, 137, 124))
colnames(data)<-c(paste('G1_', 1:n1, sep = ''), paste('G2_', 1:n2, sep = '')) 
rownames(data) <- paste('gene', 1:nrow(data), sep = '')

realdata$skin$data <- data
realdata$skin$rep <- c(n1, n2)


## Lymphoblastoid cell lines data set
n1 <- 40
n2 <- 29
param_ID <- "SRP001540"
download_study(param_ID, type = "rse-gene", download = T) 
in_f <- "SRP001540/rse_gene.Rdata" 
load(in_f)                             
rse <- rse_gene
rse <- scale_counts(rse) 
x <- assays(rse)$counts 
x <- as.data.frame(x) 
data <- cbind(x$SRR031822, x$SRR031953 + x$SRR031873,#Female1-2 
              x$SRR031952 + x$SRR031871, x$SRR031868,#Female3-4
              x$SRR031819, x$SRR031897 + x$SRR031857,#Female5-6
              x$SRR031823, x$SRR031959, x$SRR031955,#Female7-9
              x$SRR031954, x$SRR031956, x$SRR031838,#Female10-12
              x$SRR031918, x$SRR031817, x$SRR031949 + x$SRR031852,#Female13-15
              x$SRR031841, x$SRR031865, x$SRR031896,#Female16-18
              x$SRR031853, x$SRR031820, x$SRR031874,#Female19-21
              x$SRR031895, x$SRR031870, x$SRR031839,#Female22-24
              x$SRR031958, x$SRR031867, x$SRR031848,#Female25-27
              x$SRR031847, x$SRR031818, x$SRR031919,#Female28-30
              x$SRR031866, x$SRR031849, x$SRR031877,#Female31-33
              x$SRR031814, x$SRR031914, x$SRR031812,#Female34-36
              x$SRR031842, x$SRR031843, x$SRR031860, x$SRR031837,#Female37-40
              x$SRR031917, x$SRR031821 + x$SRR031898,#Male1-2
              x$SRR031950 + x$SRR031850, x$SRR031876 + x$SRR031862,#Male3-4
              x$SRR031875, x$SRR031915, x$SRR031878 + x$SRR031863,#Male5-7
              x$SRR031869, x$SRR031864, x$SRR031845,#Male8-10
              x$SRR031951 + x$SRR031851, x$SRR031846,#Male11-12
              x$SRR031916, x$SRR031844, x$SRR031813,#Male13-15
              x$SRR031894, x$SRR031854, x$SRR031858,#Male16-18
              x$SRR031859, x$SRR031872, x$SRR031816,#Male19-21
              x$SRR031815, x$SRR031920 + x$SRR031899,#Male22-23
              x$SRR031957 + x$SRR031855, x$SRR031840,#Male24-25
              x$SRR031948, x$SRR031893, x$SRR031811, x$SRR031861)#Male26-29
obj <- as.logical(rowSums(data) > 0) 
data <- data[obj,]
data <- as.data.frame(data)
colnames(data) <- c(paste("G1_", 1:n1, sep = ""), paste("G2", 1:n2, sep = "")) 
rownames(data) <- paste('gene', 1:nrow(data), sep = '')
data <- data[obj,]
data <- as.data.frame(data)
realdata$lym$data <- data
realdata$lym$rep <- c(n1, n2)


# Output
out_f <- "../realdataset.obj"
saveRDS(realdata, out_f)
