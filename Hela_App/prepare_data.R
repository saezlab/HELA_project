

load("/Users/miyang/Documents/RWTH_Aachen/HELA/DATA/Hela_data_log.Rdata")

x <- Hela_data_log$merged_v2
SP_Identified <- x$SP_Identified
table <- cbind(rownames(x) , SP_Identified ) 
colnames(table) <- c("HGNC","swissprot_ID")
table <- table[complete.cases(table) , ]

write.csv(table, "/Users/miyang/Documents/RWTH_Aachen/HELA//Hela_App/Hela_dashboard/data/swissprot_ID_HGNC")


swissprot_ID_HGNC <- read.csv("~/Documents/RWTH_Aachen/HELA/Hela_App/Hela_dashboard/data/swissprot_ID_HGNC", row.names=1)
swissprot_ID_HGNC$swissprot_ID <- as.character(swissprot_ID_HGNC$swissprot_ID)
swissprot_ID_HGNC$HGNC <- as.character(swissprot_ID_HGNC$HGNC)

matrix_0 <- Hela_data_log$Prot
matrix_0 <- matrix_0[ complete.cases(matrix_0) , ]
common <- intersect( swissprot_ID_HGNC[ ,1] , rownames(matrix_0) )
matrix <- matrix_0[common, ]
x <- swissprot_ID_HGNC[which(swissprot_ID_HGNC[ ,1] %in% common), 2 ]
rownames(matrix) <- x



