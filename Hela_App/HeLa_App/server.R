# server.R - Dashboard
library(shiny)
library(DMwR)
library(grid)
library(pheatmap)
library(ggplot2)
library(ggrepel)
# setwd("~/Documents/RWTH_Aachen/HELA/Hela_App/Hela_dashboard")
source("./functions.R")
load("./data/Hela_data_log.Rdata")
label_names <- c("No.1_Kyoto","No.2_CCL2","No.3_Kyoto","No.4_Kyoto","No.5_S3","No.6_CCL2","No.7_CCL2",
                 "No.8_Kyoto","No.9_Kyoto","No.10_Kyoto","No.11","No.12_CCL2","No.13_CCL2","No.14_CCL2")
List_of_cell_lines <- as.list(1:length(label_names)) ; names(List_of_cell_lines) <- label_names

options(shiny.maxRequestSize=30*1024^2)

shinyServer(
  function(input, output, session) {

    # Reactive function regarding the omics layer selection
    F1_1_layer_0 = reactive({
      switch(input$F1_1_layer,
             "Proteomics" = Hela_data_log$Prot[complete.cases(Hela_data_log$Prot), ],
             "mRNA" = Hela_data_log$mRNA[complete.cases(Hela_data_log$mRNA), ],
             "CNV" = Hela_data_log$CNV[complete.cases(Hela_data_log$CNV), ],
             "Kloss" = Hela_data_log$Kloss[complete.cases(Hela_data_log$Kloss), ],
             "Let7/Control" = Hela_data_log$Let7[complete.cases(Hela_data_log$Let7), ])
    })
    F1_1_layer = reactive({ F1_1_layer <- convert_ID(input$select_Identifier,F1_1_layer_0() ) })

    # Reactive function regarding the omics layer selection
    F2_1_layer_0 = reactive({
      switch(input$F2_1_layer,
             "Proteomics" = Hela_data_log$Prot[complete.cases(Hela_data_log$Prot), ],
             "mRNA" = Hela_data_log$mRNA[complete.cases(Hela_data_log$mRNA), ],
             "CNV" = Hela_data_log$CNV[complete.cases(Hela_data_log$CNV), ],
             "Kloss" = Hela_data_log$Kloss[complete.cases(Hela_data_log$Kloss), ],
             "Let7/Control" = Hela_data_log$Let7[complete.cases(Hela_data_log$Let7), ])
    })
    F2_1_layer = reactive({ 
      F2_1_layer <- convert_ID(input$select_Identifier,F2_1_layer_0() ) 
      if(length(F2_list_genes())<=1) { F2_1_layer } else { F2_1_layer <- F2_1_layer[which(rownames(F2_1_layer) %in% F2_list_genes()), ] ; F2_1_layer }
      })
    
    F3_layer_X_0 = reactive({
      switch(input$F3_layer_X,
             "Proteomics" = Hela_data_log$Prot[complete.cases(Hela_data_log$Prot), ],
             "mRNA" = Hela_data_log$mRNA[complete.cases(Hela_data_log$mRNA), ],
             "CNV" = Hela_data_log$CNV[complete.cases(Hela_data_log$CNV), ],
             "Kloss" = Hela_data_log$Kloss[complete.cases(Hela_data_log$Kloss), ],
             "Let7/Control" = Hela_data_log$Let7[complete.cases(Hela_data_log$Let7), ])
    })
    F3_layer_X = reactive({ F3_layer_X <- convert_ID(input$select_Identifier,F3_layer_X_0() ) })
    
    F3_layer_Y_0 = reactive({
      switch(input$F3_layer_Y,
             "Proteomics" = Hela_data_log$Prot[complete.cases(Hela_data_log$Prot), ],
             "mRNA" = Hela_data_log$mRNA[complete.cases(Hela_data_log$mRNA), ],
             "CNV" = Hela_data_log$CNV[complete.cases(Hela_data_log$CNV), ],
             "Kloss" = Hela_data_log$Kloss[complete.cases(Hela_data_log$Kloss), ],
             "Let7/Control" = Hela_data_log$Let7[complete.cases(Hela_data_log$Let7), ])
    })
    F3_layer_Y = reactive({ F3_layer_Y <- convert_ID(input$select_Identifier,F3_layer_Y_0() ) })
    
    # Display Omics layer
    output$F1_1_layer = DT::renderDataTable(
      F1_1_layer(),
      options = list(scrollX = TRUE)
    )
    output$F2_1_layer = DT::renderDataTable(
      F2_1_layer(),
      options = list(scrollX = TRUE)
    )
    
    F3_layer_X_common = reactive({ 
      F3_layer_X()[ intersect(rownames(F3_layer_X()), rownames(F3_layer_Y())) , ] 
    })
    
    # Display Omics layer
    output$F3_layer_X_common = DT::renderDataTable(
      F3_layer_X_common(),
      options = list(scrollX = TRUE)
    )
    
    ########## Input genes and cell lines
    gene_density  = reactive({ input$F1_1_gene })
    F2_list_genes = reactive({ 
      genes <- strsplit(input$F2_list_genes," ")[[1]]
      genes <- strsplit(genes,"\t") 
      genes <- unlist(genes)
      genes <- unique(genes) ; genes
      })
    
    gene_Scatter_plot = reactive({ input$F3_1_gene })
   
    F4_list_genes = reactive({ 
      genes <- strsplit(input$F4_list_genes," ")[[1]]
      genes <- strsplit(genes,"\t") 
      genes <- unlist(genes)
      genes <- unique(genes) ; genes
    })
    ##########################################

    F4_merged_v2 = reactive({
      Prot <- convert_ID(input$select_Identifier,Hela_data_log$Prot )
      mRNA <- convert_ID(input$select_Identifier,Hela_data_log$mRNA )
      CNV <- convert_ID(input$select_Identifier,Hela_data_log$CNV )
      Kloss <- convert_ID(input$select_Identifier,Hela_data_log$Kloss )
      Let7 <- convert_ID(input$select_Identifier,Hela_data_log$Let7 )
      m <- cbind(Prot,mRNA,CNV,Kloss,Let7) ; m
    
    })
    
    F4_1_layer_0 = reactive({
      switch(input$F4_1_layer,
             "Proteomics" = Hela_data_log$Prot[complete.cases(Hela_data_log$Prot), ],
             "mRNA" = Hela_data_log$mRNA[complete.cases(Hela_data_log$mRNA), ],
             "CNV" = Hela_data_log$CNV[complete.cases(Hela_data_log$CNV), ],
             "Kloss" = Hela_data_log$Kloss[complete.cases(Hela_data_log$Kloss), ],
             "Let7/Control" = Hela_data_log$Let7[complete.cases(Hela_data_log$Let7), ],
             "ALL layers" = F4_merged_v2() )
 
    })
    F4_1_layer = reactive({ 
      F4_1_layer <- convert_ID(input$select_Identifier,F4_1_layer_0() ) 
      if(length(F4_list_genes())<=1) { F4_1_layer } else { F4_1_layer <- F4_1_layer[which(rownames(F4_1_layer) %in% F4_list_genes()), ] ; F4_1_layer }
      })
    
    # generate table of genes output
    output$F4_1_layer = DT::renderDataTable(
      F4_1_layer(),
      options = list(scrollX = TRUE)
    )
    F4_list_genes
    
    # Plot Density
    output$density = renderPlot({
      withProgress(message="Draw Density plot",value=1,{
        df <- as.data.frame(F1_1_layer())
        density_plot(df=df,gene=gene_density(),cell=as.numeric(input$F1_1_cell),layer=input$F1_1_layer,label_names=label_names)
      })
    })
    # switch from update to results menu as soon as the go button is pushed
    observeEvent(input$go_density, {
      new.menu = switch(input$menu,"Data_Navigation"="Results")
      updateTabItems(session,"menu",new.menu)
    })
    
    # Plot Heatmap
    output$heatmap = renderPlot({
      withProgress(message="Draw Heatmap",value=1,{
        genes <- rownames(F2_1_layer())[which(rownames(F2_1_layer()) %in% F2_list_genes())]
        if(length(genes)<1) { mat_2<-F2_1_layer() } else {  mat_2 <- F2_1_layer()[ genes,  ] }
        Hela_heatmap_overlap_Prot(mat_2,row_names="Genes",col_names="Hela Cells",title=input$F2_1_layer)
      })
    })
    # switch from update to results menu as soon as the go button is pushed
    observeEvent(input$go_heatmap, {
      new.menu = switch(input$menu,"Data_Navigation"="Results")
      updateTabItems(session,"menu",new.menu)
    })
    

    # Plot scatter plot
    output$scatter_plot = renderPlot({
      withProgress(message="Draw Scatter plot",value=1,{
        if(gene_Scatter_plot() %not in% rownames(F3_layer_X_common())) {gene <- rownames(F3_layer_X_common())[1]} else {gene <- gene_Scatter_plot()}
        
        df <- cbind( F3_layer_X()[gene, ] , F3_layer_Y()[gene, ] ) 
        colnames(df) <- c("a","b") ; df <- data.frame(df)
        label_X <- Add_scale_to_layer(input$F3_layer_X)
        label_Y <- Add_scale_to_layer(input$F3_layer_Y)
        scatter_plot(gene=gene,real_gene=gene_Scatter_plot(),df=df,label_X,label_Y,text_size=22,title_size=2,label_names=label_names) 

      })
    })
    # switch from update to results menu as soon as the go button is pushed
    observeEvent(input$go_scatter_plot, {
      new.menu = switch(input$menu,"Data_Navigation"="Results")
      updateTabItems(session,"menu",new.menu)
    })
    
    # Download selected genes
    fileext = reactive({
      switch(input$datatype,
             "Excel (CSV)" = "csv", "Text (TSV)" = "txt",
             "Text (Space Separated)" = "txt", "Doc" = "doc") 
    })
    output$downloadData = downloadHandler(
      filename = function() {
        paste0("Hela_selected_genes_",input$F4_1_layer,".", fileext(), sep=".")
      },
      content = function(file) {
        sep = switch(input$datatype,
                     "Excel (CSV)" = ";", "Text (TSV)" = "\t",
                     "Text (Space Separated)" = " ", "Doc" = " " )
        write.table(F4_1_layer(), file, sep=sep )
      }
    )
    
    outputOptions(output,"heatmap",suspendWhenHidden=F)
    outputOptions(output,"scatter_plot",suspendWhenHidden=F)
  }
)

