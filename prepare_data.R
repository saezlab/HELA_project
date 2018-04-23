path <- "~/Documents/RWTH_Aachen"
source(paste(path,"/FUNCTIONS/general_functions.R", sep=""))
library(openxlsx)

x <- read.xlsx("/Users/miyang/Documents/RWTH_Aachen/HELA/ORIGINAL_DATA/merged20170531.xlsx", 1)     ## NEW
colnames(x)[2:4] <- c("Chromosome", "ChromStart", "ChromEnd")

# dup <- x[ x$gene_ %in% x$gene_[ duplicated(x$gene_) ] , ]
x <- x[ !duplicated(x$gene_) , ]
colnames(x) 
rownames(x) <- substr(x$gene_, 1, nchar(x$gene_)-1)
x <- x[ order(rownames(x)) , ] ; x<-x[,-1] 
x <- x[!is.na(x$Chromosome), ] ##  x <- x[ -grep("GL", x$gene_) , ] ; x <- x[ -grep("chrM", x$gene_) , ]
merged_v2 <- x ; merged_v2[ ,2] <- as.numeric(merged_v2[ ,2]) ; merged_v2[ ,3] <- as.numeric(merged_v2[ ,3]) 

CNV <- x[ ,grep("log2", colnames(x))] ; CNV <- CNV[ ,1:14] ; colnames(CNV) <- paste0("CNV_",sprintf("%02d", 1:14 )) ; CNV <- data.matrix(CNV)
mRNA <- x[ ,grep("mRNA", colnames(x))] ; mRNA <- mRNA[ ,1:14]; mRNA <- data.matrix(mRNA)
Prot <- x[ ,grep("Prot", colnames(x))] ; Prot <- Prot[ ,1:14]; Prot <- data.matrix(Prot)
Kloss <- x[ ,grep("Kloss", colnames(x))] ; Kloss <- Kloss[ ,3:16]; Kloss <- data.matrix(Kloss)
Let7 <- x[ ,grep("Let7", colnames(x))] ; Let7 <- Let7[ ,1:14]; Let7 <- data.matrix(Let7)
ProCopies <- x[ ,grep("ProCopies", colnames(x))] ; ProCopies <- ProCopies[ ,1:14]; ProCopies <- data.matrix(ProCopies)

Hela_data_log <- list(merged_v2, CNV, mRNA, Prot, Kloss, Let7, ProCopies)
names(Hela_data_log ) <- c("merged_v2","CNV","mRNA","Prot","Kloss","Let7","ProCopies")

save( Hela_data_log  , file="/Users/miyang/Documents/RWTH_Aachen/HELA/DATA/Hela_data_log.Rdata")

###################################

mRNA <- x[ ,grep("mRNA", colnames(x))] ; mRNA <- mRNA[ ,15:28]; mRNA <- data.matrix(mRNA)
Prot <- x[ ,grep("Prot", colnames(x))] ; Prot <- Prot[ ,15:28]; Prot <- data.matrix(Prot)
Kloss <- x[ ,grep("Kloss", colnames(x))] ; Kloss <- Kloss[ ,17:30]; Kloss <- data.matrix(Kloss)
Hela_data_AVG <- list(merged_v2, mRNA, Prot, Kloss )
names(Hela_data_AVG ) <- c("merged_v2", "mRNA","Prot","Kloss" )

save( Hela_data_AVG  , file="/Users/miyang/Documents/RWTH_Aachen/HELA/DATA/Hela_data_AVG.Rdata")

###################################
CNV <- x[ ,grep("log2", colnames(x))] ; CNV <- CNV[ ,1:14] ; colnames(CNV) <- paste0("CNV_",sprintf("%02d", 1:14 )) ; CNV <- data.matrix(CNV)
mRNA <- x[ ,grep("mRNA", colnames(x))] ; mRNA <- mRNA[ ,1:14]; mRNA <- data.matrix(mRNA)
Prot <- x[ ,grep("Prot", colnames(x))] ; Prot <- Prot[ ,1:14]; Prot <- data.matrix(Prot)
Kloss <- x[ ,grep("Kloss", colnames(x))] ; Kloss <- Kloss[ ,3:16]; Kloss <- data.matrix(Kloss)
Let7 <- x[ ,grep("Let7", colnames(x))] ; Let7 <- Let7[ ,1:14]; Let7 <- data.matrix(Let7)

CNV <-2^CNV
mRNA <- 10^mRNA
Prot <- 10^Prot
Kloss <- 10^Kloss
Let7 <- 10^Let7
ProCopies <- 10^ProCopies

Hela_data_original <- list(merged_v2, CNV, mRNA, Prot, Kloss, Let7, ProCopies)
names(Hela_data_original) <- c("merged_v2","CNV","mRNA","Prot","Kloss","Let7","ProCopies")

save(Hela_data_original, file="/Users/miyang/Documents/RWTH_Aachen/HELA/DATA/Hela_data_original.Rdata")


