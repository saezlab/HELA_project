`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

convert_ID <- function(ID=1,data) {
swissprot_ID_HGNC <- read.csv("./data/swissprot_ID_HGNC", row.names=1)
if(ID==2) { 
  common <- intersect( swissprot_ID_HGNC[ ,1] , rownames(data) )
  matrix <- data[common, ]
  x <- swissprot_ID_HGNC[which(swissprot_ID_HGNC[ ,1] %in% common), 2 ]
  rownames(matrix) <- x
  output <- matrix
} else { output <- data }
 return(output)
}

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

Hela_heatmap_overlap_Prot <- function(x,row_names,col_names,title,fontsize=14,fontsize_row=8, fontsize_col=12 ) {
  x <- x[apply(x, 1, function(y) !all(is.na(y))),] ; x <- x[ rowSums(x)!=0, ] 
#  x <- x[complete.cases(x), ] ; x <- x[ rownames(x) %in% protein_names , ] ; # x <- t(knnImputation(t(x), 10)) ; 
  x <- quantile_normalisation(x)
  colnames(x) <- as.character(1:14)

  plot_pheatmap(x, row_names=row_names, col_names=col_names, title=title ,fontsize, fontsize_row,fontsize_col, cl = T )
}

######## Label rotation pheatmap
## Edit body of pheatmap:::draw_colnames, customizing it to your liking
draw_colnames <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = 1, 
            hjust = 0.5, rot = 0, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}
## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames",
                  ns=asNamespace("pheatmap"))

plot_pheatmap <- function(mat, row_names, col_names , title , cl = T,fontsize=18,fontsize_row=18, fontsize_col=18, scale="none" ) {
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.95, name="vp", just=c("right","top"))), action="prepend")
  pheatmap(mat,  main=title, fontsize=fontsize, fontsize_row=fontsize_row,fontsize_col=fontsize_col,cluster_rows = cl, cluster_cols = cl,scale=scale,display_numbers=F)
  setHook("grid.newpage", NULL, "replace")
  grid.text(row_names, x=-0.03, rot=90, gp=gpar(fontsize=18))
  grid.text(col_names, y=0.01, gp=gpar(fontsize=18)  )
}

scatter_plot <- function(gene,real_gene,df,label_X,label_Y,text_size=26,title_size=2.5,label_names) {
  if(real_gene %not in% c(gene,"")){
    ggplot(df ) + ggtitle(paste0(real_gene," is not in common for those 2 layers"))
    
  } else {
    
  # equation, correlation and p value
  out <- cor.test(df$a,df$b) ; r <- out$estimate ; p <- out$p.value
  lm_eqn <- function(df){
    m <- lm(b ~ a, df);
    eq <- substitute(~~italic("r")~"="~r*","~~italic("p")~"="~p,
                     list(a = format(coef(m)[1], digits = 2), 
                          b = format(coef(m)[2], digits = 2), 
                          r = format(r, digits = 2),
                          p = format(p, digits=2)))
    as.character(as.expression(eq));                 
  }
  
  g <- ggplot(df, aes(a, b, color = "#000033" )) + 
    geom_point(shape = 16, size = 3, show.legend = FALSE, alpha = .8 ) + geom_smooth(method=lm,se=F,show.legend=F) + 
    labs(x = label_X, y=label_Y ) + ggtitle(gene) + 
    theme(legend.position="none",axis.text=element_text(size= text_size) , axis.title= element_text(size=text_size), plot.title = element_text(size =rel(title_size), hjust = 0.5 ),
    panel.background = element_rect(fill='white'),panel.grid.major = element_line(colour = "grey90") ) + 
    geom_text(x = min(df$a) + 0.5*(max(df$a)-min(df$a)), y = min(df$b), label = lm_eqn(df), parse = TRUE,show.legend=F,color="black",size = 7) +
    geom_text_repel(data=df ,aes(label=label_names ), size=6 ) 
  g
  }
}

Add_scale_to_layer <- function(x) {
  if (x == "Proteomics") {
    out <- "Proteomics (log10)"
  } else if (x == "mRNA") {
    out <- "mRNA (log10)"
  } else if (x == "CNV") {
    out <- "CNV (log2)"
  } else if (x == "Kloss") {
    out <- "Kloss"
  } else if (x == "Let7/Control") {
    out <- "Let7/Control"
  }
  return(out)
}

density_plot <- function(df,gene,cell=1,layer="Proteomics",text_size=20,title_size=2,label_names) {
  if(gene %not in% c(rownames(df),"")){
    df <- cbind(df[ ,cell], rep(colnames(df)[1], length(rownames(df)))) ; df <- data.frame(df) ; df$X1 <- as.numeric(as.character(df$X1))
    df <- df[complete.cases(df) , ];
    ggplot(df ) + ggtitle(paste0("The qualified ",layer," data is lacking for ",gene))
  
  } else {
    if(gene==""){gene <- rownames(df)[1]}
    gene_abundance <- df[gene,cell] 
    df <- cbind(df[ ,cell], rep(colnames(df)[1], length(rownames(df)))) ; df <- data.frame(df) ; df$X1 <- as.numeric(as.character(df$X1))
    df <- df[complete.cases(df) , ]; df_ordered <- df[ order(-df$X1) , ] ; rank <- which(df_ordered$X1 == gene_abundance)[1]
    ggplot(df, aes(x=X1)) + geom_density(fill="grey95") + ggtitle(paste0(gene," in Hela ",label_names[cell])) + xlab(Add_scale_to_layer(layer)) + 
      theme(legend.position="none",axis.text=element_text(size= text_size) , axis.title= element_text(size=text_size), 
            plot.title = element_text(size =rel(title_size), hjust = 0.5 ),
            panel.background = element_rect(fill='white'),panel.grid.major = element_line(colour ='white') )  + 
    geom_vline(xintercept=gene_abundance, color="red") + geom_text(aes(x=gene_abundance, label=paste0("Rank ",rank,"/",length(df$X1)), y=0.05), size=7 ) 
  }
  
}

# label_names <- c("No.1_Kyoto","No.2_CCL2","No.3","No.4_Kyoto","No.5_Kyoto","No.6_CCL2","No.7_CCL2","No.8_Kyoto","No.9_Kyoto","No.10_Kyoto","No.11","No.12_CCL2","No.13_CCL2","No.14_CCL2")
# density_plot(df=df, gene="",label_names=label_names)
# df <- Hela_data_log$Prot
# density_plot(df=df,gene="AAAS",cell=1,layer="CNV",label_names=label_names)


