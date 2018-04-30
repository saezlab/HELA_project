## ui.R ##
library(shinydashboard)
library(shinyBS)  

label_names <- c("No.1_Kyoto","No.2_CCL2","No.3_Kyoto","No.4_Kyoto","No.5_S3","No.6_CCL2","No.7_CCL2",
                 "No.8_Kyoto","No.9_Kyoto","No.10_Kyoto","No.11","No.12_CCL2","No.13_CCL2","No.14_CCL2")
List_of_cell_lines <- as.list(1:length(label_names)) ; names(List_of_cell_lines) <- label_names


#################### HEADER ####################
header = dashboardHeader(
  title = "Hela App",
  # Button pointing to the Saezlab Homepage
  tags$li(class = "dropdown",id = "saezlab",
          tags$a(href="http://saezlab.org", target="_blank", 
                 tags$img(src="logo_saezlab_removed_transparent_background.png",
                          height = 18, width = 18)
          ),
          bsTooltip("saezlab", "Go to the Saezlab Homepage",
                    "bottom", options = list(container = "body"))),
  
  # Button pointing to the orginal HeLa publiation
  tags$li(class = "dropdown", id="publication",
          tags$a(href="https://www.biorxiv.org/content/early/2018/04/30/307421",
                 target="_blank", tags$img(icon("file-text"))
          ),
          bsTooltip("publication", "Go to the HeLa publication",
                    "bottom", options = list(container = "body"))),
  
  # Button pointing to the GitHub respository
  tags$li(class = "dropdown", id="github",
          tags$a(href="https://github.com/saezlab/HELA_project/tree/master/Hela_App/HeLa_App",
                 target="_blank", tags$img(icon("github"))
          ),
          bsTooltip("github", "View the code of this webside on GitHub",
                    "bottom", options = list(container = "body")))
)

#################### SIDEBAR ####################
sidebar = dashboardSidebar(
  sidebarMenu(id = "menu",
              menuItem("Home", tabName = "Home", icon=icon("home")),
              menuItem("Data Navigation", tabName = "Data_Navigation", icon=icon("upload")),
              menuItem("Results", tabName = "Results", icon=icon("bar-chart")),
              menuItem("Contact", tabName = "Contact", icon=icon("group"))
  )
)

#################### BODY ####################
body = dashboardBody(heigth="auto",width=1, 
                     tabItems(
                       tabItem(tabName = "Home",
                              h1("Welcome to the HeLa App"),
                               h3(tags$i("A web application to visualize multi omics Hela cell lines")),
                               br(),
                               h4("Description:"),
                               p("There is a growing concern over the issue of reproducibility in biology research as discussed in Nature. We hypothesize that at least some of this lack of reproducibility comes from the research reagents that scientists use, particularly also cell lines and research reagents. For example, cells might have the same name but be actually different (e.g. HeLa). With arrayCGH, mRNA-Seq, our newly developed data independent acquisition method SWATH mass spectrometry, and pulsed SILAC- SWATH method, we consistently measured the gene expression at each levels for 14 different Hela cell lines across different laboratories. Collectively, the quantitative data of the CNV (whole genome-wide), transcripts (11,365 high-quality mRNAs), proteins (5030), protein turnover rates (2083), and the responsive proteomic events under Let7 mimics treatment were acquired and included in this database, demonstrating a heterogeneity of gene expression control of Hela cells from different labs."),
                               br(),
                               h4("The original Hela publication can be cited as:"),
                               p("Liu et al", 
                                 tags$a("Genomic, Proteomic and Phenotypic Heterogeneity in HeLa Cells across Laboratories: Implications for Reproducibility of Research Results",
                                        href="https://www.biorxiv.org/content/early/2018/04/30/307421",
                                        target="_blank"))
                       ),
                       tabItem(tabName = "Data_Navigation",
                               fluidRow(
                                 box(width = 10 ,
                                     title = tagList(icon("photo"),"User inputs"), status="primary", solidHeader=TRUE,
                                     
                                     tabsetPanel(
                                       tabPanel(
                                         title = tagList(icon("edit"), "Settings"),
                                         status = "primary", solidHeader = TRUE,
                                         
                                         helpText(strong("1. Select gene identifier")),
                                         selectInput("select_Identifier", label=NULL, 
                                                     choices = list("HGNC" = 1, "SwissProt ID" = 2), selected = 1)
                                         
                                       ),
                                       
                                       tabPanel(
                                         title = tagList(icon("edit"), "Abundance rank in 1 layer"),
                                         status = "primary", solidHeader = TRUE,
                                         
                                         # Select 1 gene
                                         helpText(strong("1. Enter 1 gene name")),
                                         textInput("F1_1_gene", label=NULL),
                                         
                                         # select omics
                                         helpText(strong("2. Select an omics layer")),
                                         selectInput("F1_1_layer", label=NULL,
                                                     choices = c("Proteomics","mRNA","CNV","Kloss","Let7/Control"),
                                                     selected = "Proteomics", selectize = FALSE),
                                         
                                         # select a cell line
                                         helpText(strong("3. Select a cell line")),
                                         selectInput("F1_1_cell", label = NULL, choices = List_of_cell_lines, selected = 1),
                                         
                                         actionButton("go_density", label="Density plot",icon=icon("bar-chart"))
                                         
                                       ),
                                       
                                       tabPanel(
                                         title = tagList(icon("edit"), "Heatmap in 1 layer"),
                                         status = "primary", solidHeader = TRUE,
                                         
                                         # Select genes
                                         helpText(strong("1. Enter list of genes"),("(Space Separated)")),
                                         textInput("F2_list_genes", label=NULL),
                                         
                                         # select omics
                                         helpText(strong("2. Select an omics layer")),
                                         selectInput("F2_1_layer", label=NULL,
                                                     choices = c("Proteomics","mRNA","CNV","Kloss","Let7/Control"),
                                                     selected = "Proteomics", selectize = FALSE),
                                         
                                         actionButton("go_heatmap", label="Heatmap",icon=icon("bar-chart"))
                                         
                                       ),
                                       
                                       tabPanel(
                                         title = tagList(icon("edit"), "Correlation plot between 2 layers"),
                                         status = "primary", solidHeader = TRUE,
                                         
                                         helpText(strong("1. Enter 1 gene name")), 
                                         textInput("F3_1_gene", label=NULL),
                                         
                                         helpText(strong("2. Select omics layer X")), 
                                         selectInput("F3_layer_X", label=NULL,
                                                     choices = c("Proteomics","mRNA","CNV","Kloss","Let7/Control"),
                                                     selected = "Proteomics", selectize = FALSE),
                                         
                                         helpText(strong("3. Select omics layer Y")),
                                         selectInput("F3_layer_Y", label=NULL,
                                                     choices = c("Proteomics","mRNA","CNV","Kloss","Let7/Control"),
                                                     selected = "Proteomics", selectize = FALSE),
                                         
                                         actionButton("go_scatter_plot", label="Scatter Plot",icon=icon("bar-chart"))
                                         
                                       ),
                                       
                                       tabPanel(
                                         title = tagList(icon("edit"), "Data navigation and Download"),
                                         status = "primary", solidHeader = TRUE,
                                         
                                         # Select genes
                                         helpText(strong("1. Enter list of genes"),("(Space Separated)")),
                                         textInput("F4_list_genes", label=NULL),
                                         
                                         # select omics
                                         helpText(strong("2. Select an omics layer")),
                                         selectInput("F4_1_layer", label=NULL,
                                                     choices = c("Proteomics","mRNA","CNV","Kloss","Let7/Control", "ALL layers"),
                                                     selected = "Proteomics", selectize = FALSE)
                                         
                                       )
                                     ) 
                                     
                                 ),  ## end of box
                                 
                                 tabBox(
                                   width = 10,
                                   title = tagList(icon("table"), "Data available for each input module"),
                                   id = "matrices",
                                   tabPanel("Abundance rank in 1 layer", DT::dataTableOutput("F1_1_layer")),
                                   tabPanel("Heatmap in 1 layer", DT::dataTableOutput("F2_1_layer")),
                                   tabPanel("Correlation plot between 2 layers", DT::dataTableOutput("F3_layer_X_common")),
                                   tabPanel("Data navigation and Download", DT::dataTableOutput("F4_1_layer"),
                                            helpText("Select a format and click on \'Download selected genes\'
                                                     to download the genes across all layers and for all cell lines."),
                                            radioButtons("datatype", "Select the file type:",
                                                         choices = c("Excel (CSV)",
                                                                     "Text (TSV)",
                                                                     "Text (Space Separated)",
                                                                     "Doc")),
                                            downloadButton("downloadData", label="Download selected genes")
                                            )
                               )
                               
                       )
                       ),
                       
                       tabItem(tabName = "Results",
                               
                               box(width = 6 ,
                                   title = tagList(icon("photo"),"Abundance rank in 1 layer"), status="primary",
                                   solidHeader=TRUE, plotOutput("density", height="500px")
                               ),
                               
                               box(width = 6 ,
                                   title = tagList(icon("photo"),"Heatmap in 1 layer"), status="primary",
                                   solidHeader=TRUE, plotOutput("heatmap", height="500px")
                               ),
                               
                               box(width = 6 ,
                                   title = tagList(icon("photo"),"Correlation plot between 2 layers"), status="primary",
                                   solidHeader=TRUE, plotOutput("scatter_plot", height="500px")
                               )
                               
                       ),
                       
                       tabItem(tabName = "Contact",
                               h3("Contact Us"),
                               p("Please do not hesitate to", 
                                 a("contact us", href="mailto:mi.yang@rwth-aachen.de"),
                                 " for feedback, questions or suggestions."),
                               br(),
                               h3("On Twitter"),
                               p("You can follow institute news and updates on",
                                 a("@sysbiomed", 
                                   href="https://twitter.com/sysbiomed",
                                   target="_blank"),"."),
                               br(),
                               br(),
                               br(),
                               br(),
                               a(img(src = "logo_saezlab_removed_transparent_background.png", height = 72, width = 72),
                                 href="http://saezlab.org", target="_blank"
                               ), "MTZ Pauwelstrasse 19, 52074 Aachen, Germany"
                       )
                     )
                       )

#################### JOINING ####################
tagList(dashboardPage(header, sidebar, body),
        tags$footer("HeLa App, version 0.1 (2018)", align = "center", style = "
              position:absolute;
              bottom:0;
              width:100%;
              height:50px;   /* Height of the footer */
              color: black;
              padding: 10px;
              #background-color: black;
              z-index: 1000;"))
