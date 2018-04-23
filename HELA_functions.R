
top_CV <- function(data,nb_up=500,nb_down=500) {   ### big misunderstanding, don't use this function
CV_value <- c() ; for(gene in 1:length(data[,1])) { 
  qq <- as.numeric(data[gene, ]) 
  CV_value <- c(CV_value, sd( qq-mean(qq) ))
}
x <- cbind(data, CV_value)
x <- x[complete.cases(x), ]
x <- x[ order(-x[ ,length(colnames(x))]) , ]

top <- rownames(x)[1:nb_up]
down <- rownames(x)[(length(rownames(x))-nb_down+1):length(rownames(x))]

all <- c(top,down) 
return(all)
}

###################################################### PLOTTING ######################################################

draw_by_layer <- function(folder,layer_name,layer,circle_break=0.95,min,max) {

##### add gene label 
x <- Hela_data_log$merged_v2 ; Gene <- rownames(x) ;  Gene_label <- cbind(x[,1:3], Gene) 
Gene_label <- Gene_label[ Gene_label$Gene %in% protein_names , ]
  
#### Initialize 
data(UCSC.HG19.Human.CytoBandIdeogram); 
N_track_in<-15
RCircos.Set.Core.Components(cyto.info= UCSC.HG19.Human.CytoBandIdeogram, chr.exclude="chrY",tracks.inside=N_track_in, tracks.outside=0)
#### Set parameters
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$hist.width <- 2
rcircos.params$heatmap.width <- 100
rcircos.params$sub.tracks <- 6
rcircos.params$track.height <- rcircos.params$track.height*1.2
RCircos.Reset.Plot.Parameters(rcircos.params)
cyto.band <- RCircos.Get.Plot.Ideogram()
cyto.band$BandColor <- "black" ; cyto.band$ChrColor <- "black"
RCircos.Reset.Plot.Ideogram(cyto.band)

params <- RCircos.Get.Plot.Parameters()
params$heatmap.color <- "BlueWhiteRed"
RCircos.Reset.Plot.Parameters(params)

setwd(paste0("/Users/miyang/Documents/RWTH_Aachen/HELA/PLOT/",folder) )

# png(paste0(layer_name,".png"), width = 3000, height = 3000 ,pointsize = 60)
# tiff(paste0(layer_name,".tiff"), width = 3000, height = 3000 ,pointsize = 60)
# pdf(file=paste0(layer_name,".pdf"), height=12, width=12, compress=TRUE)

##### break the circle
RCircos.ZoomOut.Chromosome(zoom.out.ratio=circle_break) 
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

##### add gene label 
x <- Hela_data_log$merged_v2 ; Gene <- rownames(x) ;  Gene_label <- cbind(x[,1:3], Gene) 

##### add all samples
track <- 1
for(i in Hela_total_order_Liu) {
  x <- layer ; common_gene <- intersect(rownames(Gene_label), rownames(x) ) 
  x <- cbind(Gene_label[ common_gene , ], x[ common_gene ,i]) ; x <- x[!is.na(x[ ,5]), ]
  x$Chromosome <- paste0("chr", x$Chromosome)
#  x <- x[ x$Chromosome == "chr2" , ]
  RCircos.Heatmap.Plot( x, data.col=5,track.num=track, side="in" ,min.value=min,max.value=max)  ; track <- track + 1 
# RCircos.Histogram.Plot( x, data.col=5,track.num=track, side="in",min.value=min,max.value=max )  ; track <- track + 1 
}

title (paste0(layer_name), line = -1.8 , cex.lab=1.1, cex.axis=1.1, cex.main=1.1, cex.sub=1.1)

RCircos.Plot.Heatmap.Color.Scale(max.value=max, min.value=min, color.type="BlueWhiteRed", scale.location=2,scale.width=1, scale.height=0.1) 
# legend(0.6, 2, legend=c("Right: GSE42081 Mouse","Left:  GSE42081 Rat"), cex=0.8)
# dev.off()
}

####################################################################################################################

draw_by_layer_CV <- function(folder,layer_name,layer,circle_break=0.96,min,max) {
  
  ##### add gene label 
  x <- Hela_data_original$merged_v2 ; Gene <- rownames(x) ;  Gene_label <- cbind(x[,1:3], Gene) 
  Gene_label <- Gene_label[ Gene_label$Gene %in% protein_names , ]
  
  #### Initialize 
  data(UCSC.HG19.Human.CytoBandIdeogram); 
  N_track_in<-3
  RCircos.Set.Core.Components(cyto.info= UCSC.HG19.Human.CytoBandIdeogram, chr.exclude="chrY",tracks.inside=N_track_in, tracks.outside=0)
  #### Set parameters
  rcircos.params <- RCircos.Get.Plot.Parameters()
  rcircos.params$track.height <- rcircos.params$track.height*2.5
  RCircos.Reset.Plot.Parameters(rcircos.params)
  
  cyto.band <- RCircos.Get.Plot.Ideogram()
  cyto.band$BandColor <- "black" ; cyto.band$ChrColor <- "black"
  RCircos.Reset.Plot.Ideogram(cyto.band)
  
  setwd(paste0("/Users/miyang/Documents/RWTH_Aachen/HELA/PLOT/",folder) )
  
  ##### break the circle
  RCircos.ZoomOut.Chromosome(zoom.out.ratio=circle_break) 
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()

  ##### add all samples
  track <- 1
  for(i in 1:length(CV_group)) {
    x <- layer ; x <- cbind(Gene_label, x[rownames(Gene_label),CV_group[[i]]]) ; x <- x[complete.cases(x), ] 
    CV_value <- c() ; for(gene in 1:length(x[,1])) { 
      qq <- as.numeric(x[gene,5:length(colnames(x))]) 
      CV_value <- c(CV_value, sd(qq) )
    }
    CV_value[CV_value>0.5] <- 0.5
    x <- cbind(x, CV_value)
    x <- x[complete.cases(x), ] ; print(min(x$CV_value)) ; print(max(x$CV_value))
    x$Chromosome <- paste0("chr", x$Chromosome)
    RCircos.Heatmap.Plot( x, data.col=length(colnames(x)),track.num=track, side="in",min.value=min,max.value=max) ; track <- track + 1
  }
  title (paste0(layer_name), line = -2 , cex.lab=1.1, cex.axis=1.1, cex.main=1.1, cex.sub=1.1)
  RCircos.Plot.Heatmap.Color.Scale(max.value=max, min.value=min,color.type="BlueWhiteRed", scale.location=2,scale.width=1, scale.height=0.1) 

}

####################################################################################################################


draw_by_layer_CV_UniGroup_RANK_range01 <- function(folder,layer_list,circle_break=0.96 ,color_choice="BlueWhiteRed") {
  
  ##### add gene label 
  x <- Hela_data_original$merged_v2 ; Gene <- rownames(x) ;  Gene_label <- cbind(x[,1:3], Gene) 
  Gene_label <- Gene_label[ Gene_label$Gene %in% protein_names , ]

  #### Initialize 
  data(UCSC.HG19.Human.CytoBandIdeogram)
  N_track_in<-6
  RCircos.Set.Core.Components(cyto.info= UCSC.HG19.Human.CytoBandIdeogram, chr.exclude="chrY",tracks.inside=N_track_in, tracks.outside=0)
  #### Set parameters
  rcircos.params <- RCircos.Get.Plot.Parameters()
  rcircos.params$heatmap.color <- color_choice
  rcircos.params$track.height <- rcircos.params$track.height*2.5
  RCircos.Reset.Plot.Parameters(rcircos.params)
  
  cyto.band <- RCircos.Get.Plot.Ideogram()
  
  cyto.band$BandColor <- "white" ; cyto.band$ChrColor <- "white"
  RCircos.Reset.Plot.Ideogram(cyto.band)
  
  setwd(paste0("/Users/miyang/Documents/RWTH_Aachen/HELA/PLOT/",folder) )
  
  ##### break the circle
  RCircos.ZoomOut.Chromosome(zoom.out.ratio=circle_break) 
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()
  
  ##### add all samples
  track <- 1
  for(i in 1:4) {
    x <- layer_list[[i]] ; x <- cbind(Gene_label, x[rownames(Gene_label), 1:14 ]) ; x <- x[complete.cases(x), ] 
    CV_value <- c() ; for(gene in 1:length(x[,1])) { 
      qq <- as.numeric(x[gene,5:length(colnames(x))]) 
      CV_value <- c(CV_value, sd(qq) )
    }
    
    CV_value <- rank(CV_value) ; CV_value <- range_01(CV_value)
 
    x <- cbind(x, CV_value)
    x <- x[complete.cases(x), ] 
    x$Chromosome <- paste0("chr", x$Chromosome)
    RCircos.Heatmap.Plot( x, data.col=length(colnames(x)),track.num=track, side="in",min.value=min(CV_value),max.value=max(CV_value)) ; track <- track + 1
    #  RCircos.Heatmap.Plot( x, data.col=length(colnames(x)),track.num=track, side="in",min.value=min(x$CV_value),max.value=max(x$CV_value)) ; track <- track + 1
  }
  title ("Covariation among all 14 cell lines", line = -3 , cex.lab=1.1, cex.axis=1.1, cex.main=1.1, cex.sub=1.1)
  
  RCircos.Plot.Heatmap.Color.Scale(max.value=1, min.value=0,color.type=color_choice, scale.location=2,scale.width=1, scale.height=0.1) 
  # legend(0.6, 2, legend=c("Right: GSE42081 Mouse","Left:  GSE42081 Rat"), cex=0.8)
  
}


####################################################################################################################

draw_by_layer_ratio <- function(folder,layer_list,circle_break=0.96,min,max) {
  
  ##### add gene label 
  x <- Hela_data_log$merged_v2 ; Gene <- rownames(x) ;  Gene_label <- cbind(x[,1:3], Gene) 
  Gene_label <- Gene_label[ Gene_label$Gene %in% protein_names , ]
  
  #### Initialize 
  data(UCSC.HG19.Human.CytoBandIdeogram); 
  N_track_in<-6
  RCircos.Set.Core.Components(cyto.info= UCSC.HG19.Human.CytoBandIdeogram, chr.exclude="chrY",tracks.inside=N_track_in, tracks.outside=0)
  #### Set parameters
  rcircos.params <- RCircos.Get.Plot.Parameters()
  rcircos.params$track.height <- rcircos.params$track.height*2.5
  RCircos.Reset.Plot.Parameters(rcircos.params)
  
  cyto.band <- RCircos.Get.Plot.Ideogram()
  cyto.band$BandColor <- "black" ; cyto.band$ChrColor <- "black"
  RCircos.Reset.Plot.Ideogram(cyto.band)
  
  setwd(paste0("/Users/miyang/Documents/RWTH_Aachen/HELA/PLOT/",folder) )
  
  ##### break the circle
  RCircos.ZoomOut.Chromosome(zoom.out.ratio=circle_break) 
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()
  
  ##### add all samples
  track <- 1
  for(i in 1:4) {  # i=4
    x <- layer_list[[i]] ; x <- cbind(Gene_label, x[rownames(Gene_label), 1:14 ]) ; x <- x[complete.cases(x), ] 
    
    mean_CCL2 <- c() ; for(gene in 1:length(x[,1])) { 
      qq <- as.numeric(x[gene,Hela_CCL2+4]) 
      mean_CCL2 <- c(mean_CCL2, mean(qq) )
    }
    mean_Kyoto <- c() ; for(gene in 1:length(x[,1])) { 
      qq <- as.numeric(x[gene,Hela_Kyoto+4]) 
      mean_Kyoto <- c(mean_Kyoto,mean(qq) )
    }
    mean_Kyoto[mean_Kyoto==0] <- 1E-10
    ratio_of_mean <- mean_CCL2/mean_Kyoto
    ratio_of_mean <- log2(ratio_of_mean)
    print(mean(ratio_of_mean))
    ratio_of_mean[ratio_of_mean > 1] <- 1
    ratio_of_mean[ratio_of_mean < -1] <- -1

    x <- cbind(x, ratio_of_mean)
    x <- x[complete.cases(x), ]  
    x$Chromosome <- paste0("chr", x$Chromosome)
    RCircos.Heatmap.Plot( x, data.col=length(colnames(x)),track.num=track, side="in",min.value=min,max.value=max) ; track <- track + 1
  }
  title ("Ratio of CCL2/Kyoto", line = -3 , cex.lab=1.1, cex.axis=1.1, cex.main=1.1, cex.sub=1.1)
  
  RCircos.Plot.Heatmap.Color.Scale(max.value=max, min.value=min,color.type="BlueWhiteRed", scale.location=2,scale.width=1, scale.height=0.1) 
  # legend(0.6, 2, legend=c("Right: GSE42081 Mouse","Left:  GSE42081 Rat"), cex=0.8)
}

draw_1_layer <- function(folder,title,data,circle_break=0.96,min,max) {
  
  ##### add gene label 
  x <- Hela_data_log$merged_v2 ; Gene <- rownames(x) ;  Gene_label <- cbind(x[,1:3], Gene) 
  Gene_label <- Gene_label[ Gene_label$Gene %in% protein_names , ]
  
  #### Initialize 
  data(UCSC.HG19.Human.CytoBandIdeogram); 
  N_track_in<-6
  RCircos.Set.Core.Components(cyto.info= UCSC.HG19.Human.CytoBandIdeogram, chr.exclude="chrY",tracks.inside=N_track_in, tracks.outside=0)
  #### Set parameters
  rcircos.params <- RCircos.Get.Plot.Parameters()
  rcircos.params$track.height <- rcircos.params$track.height*2.5
  RCircos.Reset.Plot.Parameters(rcircos.params)
  
  cyto.band <- RCircos.Get.Plot.Ideogram()
  cyto.band$BandColor <- "black" ; cyto.band$ChrColor <- "black"
  RCircos.Reset.Plot.Ideogram(cyto.band)
  setwd(paste0("/Users/miyang/Documents/RWTH_Aachen/HELA/PLOT/",folder) )
  
  pdf(file=paste0(title,".pdf"), height=12, width=12, compress=TRUE)
  
  ##### break the circle
  RCircos.ZoomOut.Chromosome(zoom.out.ratio=circle_break) 
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()
  
  ##### add all samples
  track <- 1
  common_gene <- intersect(names(data),rownames(Gene_label))
  Gene_label <- Gene_label[rownames(Gene_label) %in% common_gene, ]
  data <- data[common_gene]
  x <- cbind(Gene_label, data) ; x <- x[complete.cases(x), ] 
  
  x$Chromosome <- paste0("chr", x$Chromosome)
  RCircos.Heatmap.Plot( x, data.col=length(colnames(x)),track.num=track, side="in",min.value=min,max.value=max)  

  title (title, line = -3 , cex.lab=1.1, cex.axis=1.1, cex.main=1.1, cex.sub=1.1)
  
  RCircos.Plot.Heatmap.Color.Scale(max.value=max, min.value=min,color.type="BlueWhiteRed", scale.location=2,scale.width=1, scale.height=0.1) 
  
  dev.off()
}

####################################################################################################################
library(DMwR)
Hela_heatmap_overlap_Prot <- function(x,folder,row_names,col_names,name_file,title,fontsize=24,fontsize_row=22, fontsize_col=22) {
  setwd(paste0("/Users/miyang/Documents/RWTH_Aachen/HELA/PLOT/",folder) )

  x <- x[apply(x, 1, function(y) !all(is.na(y))),] ; x <- x[ rowSums(x)!=0, ] 
  x <- x[complete.cases(x), ] ; x <- x[ rownames(x) %in% protein_names , ] ; # x <- t(knnImputation(t(x), 10)) ; 
  x <- quantile_normalisation(x)
  colnames(x) <- paste0("sample_",as.character(1:14))
 
  pdf(name_file, width = 10 , height = 14, onefile = F )
  plot_pheatmap(x, row_names=row_names, col_names=col_names, title=title ,fontsize, fontsize_row,fontsize_col, cl = T)
  dev.off()
}

