# setwd("~/Documents/RWTH_Aachen/HELA/Hela_App/Hela_dashboard")
load("./data/Hela_data_log.Rdata")
Prot <- Hela_data_log$Prot ; Prot <- Prot[complete.cases(Prot) , ]
df <- Prot
df <- data.frame(df)


genes <- "A1BG  AAAS	AAGAB  AAMDC "
genes <- strsplit(genes," ")[[1]]
genes <- strsplit(genes,"\t") 
genes <- unlist(genes)

genes <- gsub('[\r\n\t]', '',genes)

genes <- trimws(genes, which = c("both"))
genes <- strsplit(genes," ")[[1]]
genes <- genes[which(genes!="")] ; genes

choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3)
label_names <- c("No.1_Kyoto","No.2_CCL2","No.3","No.4_Kyoto","No.5_Kyoto","No.6_CCL2","No.7_CCL2",
                 "No.8_Kyoto","No.9_Kyoto","No.10_Kyoto","No.11","No.12_CCL2","No.13_CCL2","No.14_CCL2")

L <- as.list(1:length(label_names)) ; names(L) <- label_names

merged_v2 <- Hela_data_log$merged_v2

############################################# put app on server #############################################

library(rsconnect)
rsconnect::deployApp("/Users/miyang/Documents/RWTH_Aachen/HELA/Hela_App/Hela_dashboard")

rsconnect::setAccountInfo(name='miyang',
                          token='41C85E8B64978827D4D680FA802AC3B5',
                          secret='/ZUMktl1JEy/z4sIsmEPcNkszKcTbzo+u/UJDPBH')

##  https://miyang.shinyapps.io/hela_dashboard/
##  https://www.shinyapps.io/admin/#/dashboard


