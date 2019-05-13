suppressMessages(library(shiny))
suppressMessages(library(shinyjs))
suppressMessages(library(shinyBS))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinyFiles))
suppressMessages(library(DT))
suppressMessages(library(rhandsontable))
suppressMessages(library(shinycssloaders))
library(plotly)
library(xtable)

appCSS <- "
.multicol { 
                                   height: 150px;
-webkit-column-count: 5; /* Chrome, Safari, Opera */ 
-moz-column-count: 5;    /* Firefox */ 
column-count: 5; 
-moz-column-fill: auto;
-column-fill: auto;
} 
#loading-content {
position: absolute;
background: white !important;
opacity: 0.8;
z-index: 1000000;
left: 0;
right: 0;
top: 0;
bottom: 0;
font-size: 50px;
text-align: center;
color: #black;
}
#loading-gif { 
opacity: 0.8; 
display: block;
margin-left: auto;
margin-right: auto;
vertical-align: middle;
z-index: 1000000;
}
"

fluidPage(
  useShinyjs(),
  useShinyjs(),
  inlineCSS(appCSS),
  hidden(div(id="loading-content",
             img(id="loading-gif", src="screen-loading.gif"),
             p(id="loadingText", "WORKING"),
             p("...")
  )),
  dashboardPage(
    dashboardHeader(title="BMDx: dose response analysis for expression data", titleWidth="35%"),
    dashboardSidebar(disable=FALSE,
                     bsCollapse(id="bsSidebar1", open="LOAD PHENOTYPE DATA",
                                bsCollapsePanel("LOAD PHENOTYPE DATA", style="warning",
                                                fluidRow(
                                                  column(12, align="center",
                                                         shinyBS::bsButton("import_pheno_submit", label="Import Phenotype Data", style="danger",icon=icon("exclamation-circle")),
                                                         shinyBS::bsTooltip("import_pheno_submit", "Launch a graphical window, to configure import of phenotype data from a file!", placement="bottom")
                                                  )
                                                )
                                ),
                                bsCollapsePanel("LOAD EXPRESSION MATRIX", style="danger",
                                                fluidRow(
                                                  column(12, align="center",
                                                         shinyBS::bsButton("import_expr_submit", label="Import Expression Matrix", style="danger", icon=icon("exclamation-circle")),
                                                         shinyBS::bsTooltip("import_expr_submit", "Launch a graphical window, to configure import of expression matrix from a file!", placement="bottom")
                                                         
                                                  )
                                                )
                                ),
                                bsCollapsePanel("ANOVA FILTERING", style="danger",
                                                fluidRow(
                                                  column(12, align="center",
                                                         shinyBS::bsButton("anova_filtering_button", label="Anova Filtering", style="danger", icon=icon("exclamation-circle")),
                                                         shinyBS::bsTooltip("anova_filtering_button", "Launch a graphical window, to configure anova parameters!", placement="bottom")
                                                  )
                                                ),
                                                
                                                fluidRow(
                                                  column(12, align="center",
                                                         shinyBS::bsButton("skip_anova_filtering_button", label="Skip Anova", style="danger", icon=icon("exclamation-circle")),
                                                         shinyBS::bsTooltip("skip_anova_filtering_button", "Skip Anova Filtering!", placement="bottom")
                                                         
                                                  )
                                                )       
                                                
                                                # fluidRow(
                                                #   column(12, align="center",
                                                #          shinyBS::bsButton("skip_anova_filtering_button", label="Skip Anova Filtering", style="danger", icon=icon("exclamation-circle")),
                                                #          shinyBS::bsTooltip("skip_anova_filtering_button", "Skip the anova filtering!", placement="bottom")
                                                #   )
                                                # )
                                ),
                                bsCollapsePanel("COMPUTE BMD", style="danger",
                                                fluidRow(
                                                  column(12, align="center",
                                                         shinyBS::bsButton("bmd_button", label="Compute BMD", style="danger", icon=icon("exclamation-circle")),
                                                         shinyBS::bsTooltip("bmd_button", "Launch a graphical window to configure BMD analysis!", placement="bottom")
                                                  )
                                                )
                                ),
                                bsCollapsePanel("PATHWAY ENRICHMENT", style="danger",
                                                fluidRow(
                                                  column(12, align="center",
                                                         shinyBS::bsButton("enrich_button", label="Enrichment", style="danger", icon=icon("exclamation-circle")),
                                                         shinyBS::bsTooltip("enrich_button", "Launch a graphical window to configure enrichment analysis!", placement="bottom")
                                                  )
                                                )
                                )
                                # bsCollapsePanel("TIMEPOINTS", style="danger",
                                #                 fluidRow(
                                #                   column(12, align="center",
                                #                          shinyBS::bsButton("tp_button", label="TP Comparison", style="danger", icon=icon("exclamation-circle")),
                                #                          shinyBS::bsTooltip("tp_button", "Launch a graphical to window to compare timepoints!", placement="bottom")
                                #                   )
                                #                 )
                                # )
                     )

    ),
    dashboardBody(
      shinyBS::bsModal("enrichPathways", "Enrichment", "enrich_button", size="large",
                       tags$h5("1. Input gene lists"),
                       wellPanel(
                         fluidRow(
                           column(6,
                                  radioButtons("organism","1) Organisms",
                                               choices = c(human = "Human", mouse = "Mouse", rat = "Rat"),selected = "Human"),
                                  shinyBS::bsTooltip(id = "organism",title = "Note: select organism and gene ID before uploading the file",placement = "bottom")
                                  
                           ) ,
                           column(6,radioButtons("idtype","2) GeneID",
                                                 choices = c(symbols = "SYMBOL", ensemble = "ENSEMBL",entrez = "ENTREZID"),
                                                 selected = "SYMBOL"),
                                  shinyBS::bsTooltip(id = "idtype",title = "Note: select organism and gene ID before uploading the file",placement = "bottom")
                                  
                                  
                           )
                           # column(4,
                           #        fileInput("file1", "3) Choose Excel File",
                           #                  multiple = FALSE,
                           #                  accept = c("text/csv/xlsx",
                           #                             "text/comma-separated-values,text/plain/excel",
                           #                             ".xlsx")))
                           
                         )
                       ),
                       tags$h5("2. Functional annotation parameters"),
                       wellPanel(
                         fluidRow(
                           column(6,radioButtons("EnrichType","Select Functional Annotation",
                                                 choices = c(KEGG = "KEGG", REACTOME="REACTOME",GO = "GO"),
                                                 selected = "KEGG")
                           ),
                           column(6,radioButtons("GOType","Select GO",
                                                 choices = c(BP = "BP", CC="CC",MF = "MF"),
                                                 selected = "BP")
                           )),
                         fluidRow(
                           column(6,radioButtons("pcorrection","Correction Method",
                                                 choices = c(gSCS = "analytical", fdr = "fdr", bonferroni = "bonferroni", Nominal = "none"),
                                                 selected = "analytical"),
                                  shinyBS::bsTooltip(id = "pcorrection",
                                                     title = "Default is g:SCS. Check g:Profiler web page for more info",
                                                     placement = "top")
                           ),
                         column(6,selectInput(inputId = "pvalueTh", label = "P-value threshold:",choices = list(0.001,0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09),selected = 0.05)
                                
                                
                         )),
                         fluidRow(
                           column(6,checkboxInput("only_annotated", "Annotated genes only", value = TRUE)),
                           column(6, sliderInput("min_intersection", "Minumum number of genes in the intersection:",
                                                          min = 0, max = 100, value = 0))
                         )
                       ),
                       tags$h5("3. Display parameters"),
                       wellPanel(
                         # Horizontal line ----
                         fluidRow(
                           
                           column(6,radioButtons("aggregation","Aggregation Function",
                                                 choices = c(min = "min", max = "max",mean = "mean",
                                                             median = "median"),
                                                 selected = "mean")  
                           ),
                          
                           column(6,
                                  radioButtons("MapValueType","Choose Values Type",
                                               choices = c(Pvalue = "PVAL", GenesModifications="FC",GenesModifications_PValue  = "FCPV"),
                                               selected = "FC")
                           )
                           ),
                         
                         fluidRow(column(6,radioButtons("continuous","Plot modification",
                                                        choices = c(value = "continuous",
                                                                    sign = "discrete"),
                                                        selected = "continuous")),
                                 # column(6,actionButton("computePathways","Generate Map"))
                                  column(6,  shinyBS::bsButton("enrichment_analysis", label="Run Enrichment", style="info", icon=icon("hand-o-right")))
                                 )
                        
                       )
      ),
      shinyBS::bsModal("computeBMD", "Compute BMD Value", "bmd_button", size="large",
                       fluidRow(
                         column(3,
                                uiOutput("MaxDose")
                         ),column(3,
                                  selectInput("LOOF", "Lack-of-fit PValue Th:", choices=c(0.3,0.2,0.1,0.05),selected=0.1)
                        ),
                        # column(3,
                        #        selectInput("Interval", "Interval type:", choices=c("none","delta","fls"),selected="delta")
                        #  ),
                        #column(3, selectInput("RespLev","Response Level",choices = seq(10,90,10), selected = 10))
                        column(3, textInput("RespLev", label = "Response Level", value =1.349))
                        
                       ),
                       fluidRow(
                         column(12,selectInput("BMDSettings", "Select the BMD analysis setting", choices=c("All","Regulatory","Degree of Freedom", "Custom"),selected="Custom"))
                       ),
                       fluidRow(
                         column(12,
                                wellPanel(
                                    uiOutput("bmd_checkbox")
                                )
                                )
                       ),
                       # fluidRow(
                       #   column(3, textInput("Power", label = "Power value", value =2)),
                       #   column(3, textInput("Hill_pow", label = "Hill power", value =2)),
                       #   column(3, textInput("Hill_kd", label = "Hill Kd", value =10))
                       # ),
                       
                       fluidRow(
                         column(12, align="right",
                                shinyBS::bsButton("bmd_analysis", label="Run BMD", style="info", icon=icon("hand-o-right"))
                         )
                       ),
                       fluidRow(
                         column(12,
                         h5("Models Description"),                     
                                  bsCollapse(id="BMDModHelp",
                                      bsCollapsePanel("Linear Model", style="warning",
                                        withMathJax(),
                                        helpText("The formula for the linear model is
                                                  $$ f(dose) = \\beta_0 + \\beta_1 dose $$
                                                  The linear model is a special case of the polynomial model, with $$n=1$$")
                                        ),
                                        bsCollapsePanel("Polynomial Model (Quadratic/Cubic)", style="warning",
                                                  withMathJax(),
                                                  helpText("The formula for the polynomial model is
                                                           $$ f(dose) = \\beta_0 + \\beta_1 dose + \\beta_2 dose^2 + \\ldots + \\beta_n dose^n $$
                                                           Here n is the degree of the polynomial. The user can choose between $$n = 2, 3$$")
                                        ),
                                        bsCollapsePanel("Power Model", style="warning",
                                                      withMathJax(),
                                                      helpText("The formula for the power model is
                                                               $$ f(dose) = \\beta_0 + (dose)^\\delta $$
                                                               The user can choose between $$\\delta = 2, 3, 4$$")                                           
                                        ),
                                        bsCollapsePanel("Exponential Model", style="warning",
                                                      withMathJax(),
                                                      helpText("The formula for the exponential model is
                                                               $$ f(dose) = \\beta_0 + exp(dose) $$")                                           
                                                      ),
                                         bsCollapsePanel("Hill Model", style="warning",
                                                      withMathJax(),
                                                      helpText("The formula for the hill model is 
                                                               $$ f(dose) = \\beta_0 + \\dfrac{dose^n}{Kd + dose^n} $$
                                                               The user can choose between $$n = 0.5,1,2,3,4,5$$ while Kd is fixed to 10.")
                                                              
                                                      ),
                                        bsCollapsePanel("Asymptotic Regression", style = "warning",
                                                        withMathJax(),
                                                        helpText("The formula for the asymptotic regression model is the following:
                                                              $$f(dose) = c + (d-c) \\times (1-exp(-dose/e)) $$
                                                               The parameter c is the lower limit (at x=0), the parameter d is the upper limit and 
                                                               the parameter e>0 is determining the steepness of the increase of dose.
                                                              The AR.3 model is the one depending from c, d and e parameters. The AR.2 model depends only on d and e parameters, while c is set to zero")
                                                      ),
                                      bsCollapsePanel("Michaelis-Menten Model", style = "warning",
                                                      withMathJax(),
                                                      helpText("The model is defined by the three-parameter model (MM.3) function
                                                               $$f(dose, (c, d, e)) = c + \\dfrac{d-c}{1+(e/dose)}$$
                                                               It is increasing as a function of the dose, attaining the lower limit c at dose 0 (x=0) and the upper limit d for infinitely large doses. 
                                                               The parameter e corresponds to the dose yielding a response halfway between c and d. 
                                                               The common two-parameter Michaelis-Menten model (MM.2) is obtained by setting c equal to 0.")
                                      )
                                    )
                         )
                         )
                       
      ),
      shinyBS::bsModal("computeAnova", "Filter Genes by Anova", "anova_filtering_button", size="large",
                       fluidRow(
                        column(3,
                                  uiOutput("timePoint")
                         ),column(3,
                                  radioButtons("adjBool", "Pvalue:",
                                               c("FDR" = TRUE,
                                                 "Nominal" = FALSE))
                         ),
                        column(3,
                               selectInput("anovaPvalTh", "Anova PValue Th:", choices=c(0.05,0.04,0.03,0.02,0.01))
                        )
                       ),fluidRow(
                         column(12, align="right",
                                shinyBS::bsButton("anova_analysis", label="Run Anova", style="info", icon=icon("hand-o-right"))
                         )
                       )                       
      ),
      shinyBS::bsModal("importGxModal", "Import Gene Expression Table", "import_expr_submit", size="large",
                       fluidRow(
                         column(3,
                                fileInput("gx", label="File")
                         ),column(3,
                                  uiOutput("selGxSep")
                         ),column(3,
                                  textInput("gxSepT", "Other Seperator", value=":")
                         ),column(3,
                                  uiOutput("selGxQuote")
                         )
                       ),fluidRow(
                         column(12, align="right",
                                shinyBS::bsButton("upload_gx_submit", label="Import", style="info", icon=icon("hand-o-right"))
                         )
                       )
      ),
      shinyBS::bsModal("importPhenoModal", "Import Phenotype Data", "import_pheno_submit", size="large",
                       fluidRow(
                         column(3,
                                fileInput("fPheno", label="Phenotype File")
                         ),column(3,
                                  uiOutput("selSep")
                         ),column(3,
                                  textInput("sepT", "Other Seperator", value=":")
                         ),column(3,
                                  uiOutput("selQuote")
                         )
                       ),fluidRow(
                         column(1,
                                actionButton("load_pheno_submit", "Preview")
                         ),column(2, align="left",
                                  textOutput("phRowsText"),
                                  textOutput("phColsText")
                         )
                       ),hr(),
                       fluidRow(
                         column(12,
                                hidden(div(id="phenoPreviewDiv",
                                           fluidRow(
                                             column(12,
                                                    rhandsontable::rHandsontableOutput("phenoTypesRH")
                                             )
                                           ),hr(),
                                           fluidRow(
                                             column(4,
                                                      uiOutput("selSampleIDCol")
                                             ),column(4,
                                                      uiOutput("selDoseCol")
                                             ),column(4,
                                                      uiOutput("selTPCol")
                                             )
                                           ),fluidRow(
                                             column(12, align="right",
                                                    shinyBS::bsButton("upload_pheno_submit", label="Import", style="info", icon=icon("hand-o-right"))
                                             )
                                           )
                                ))
                         )
                       )
      ),
      fluidRow(
        tabBox(id="display", title="", width=12,
               tabPanel(value="pdTab", title="Phenotype Data",
                        fluidRow(
                          column(12,
                                 DT::dataTableOutput("filtered")
                          )
                        )
               ),
               tabPanel(value="gExpTab", title="Gene Expression Matrix",
                        fluidRow(
                          column(12,
                                 DT::dataTableOutput("gExpMat")
                          )
                        )
               ),
               tabPanel(value="AnovaTab", title="Anova",
                        fluidRow(
                          column(6,uiOutput("timePointSel")),
                          column(6,downloadButton("downloadAnovaData", "Download"))
                        ),
                        fluidRow(
                          column(6,
                             #fluidRow(downloadButton("downloadAnovaData", "Download")),
                             fluidRow(DT::dataTableOutput("Anova_table"))
                          ),
                          column(6,
                                 fluidRow(plotOutput("anovaPlot"))
                                 #fluidRow(plotOutput("AnovaVenn"))
                          )
                        )
               ),
               tabPanel(value="BMDTab", title="BMD",
                        tabsetPanel(
                          tabPanel("Gene Level", 
                                   fluidRow(
                                     column(6,uiOutput("timePointSel2")),
                                     column(6,fluidRow(downloadButton("downloadBMDData", "Download")))
                                   ),
                                   fluidRow(
                                     column(6,
                                            fluidRow(DT::dataTableOutput("BMD_table"))
                                     ),
                                     column(6,
                                            #checkboxInput("xlog", label = "Log scale x axis", value = TRUE),
                                            #selectInput("genePlotType", "Plot type", choices=c("average","none","obs","all","bars","confidence"),selected = "bars"),
                                            plotOutput("bmd_fitting")
                                   ))         
                          ),
                          # tabPanel("Pathway Level", 
                          #          fluidRow(
                          #            column(6, uiOutput("timePointSelPat"))
                          #          ),
                          #          fluidRow(
                          #            column(12,
                          #                   DT::dataTableOutput("PAT_table")
                          #            )),
                          #          fluidRow(
                          #            column(12,
                          #                   plotOutput("path_bmd_dist"))
                          #          )
                          #          ),
                          tabPanel("Compare TP",
                                   tabPanel("Compare Time points",
                                            fluidRow(
                                              bsCollapse(id="TPSidebar", open="BMD",
                                                         bsCollapsePanel("BMD Values", style="warning",
                                                                         #"Here i will print the BMD value distribution",
                                                                         plotlyOutput("BMD_dist_TP")
                                                         ),
                                                         bsCollapsePanel("Lack of fit Pvalues", style="warning",
                                                                         #"Here i will print the BMD value distribution",
                                                                         plotlyOutput("BMD_pval_fitting")
                                                                         
                                                         ),
                                                         bsCollapsePanel("BMD/BMDL", style="warning",
                                                                         #"Here i will print the BMD value distribution",
                                                                         plotlyOutput("BMD_BMDL")

                                                         ),
                                                         bsCollapsePanel("Fitted models", style="warning",
                                                                         #"Here i will print the BMD value distribution",
                                                                         fluidRow(column(12,plotOutput("BMD_dist_models"))),
                                                                         fluidRow(column(12, plotlyOutput("BMD_BMDL_BMDU_by_model")))
                                                         ),
                                                         bsCollapsePanel("GENE X TP", style="warning",
                                                                         #"Here i will plot the BMD for every gene at different TP",
                                                                         plotlyOutput("NGTime")
                                                         ),
                                                         bsCollapsePanel("Venn diagram responsive genes", style="warning",
                                                                         #"Here i will plot the BMD for every gene at different TP",
                                                                         # uiOutput("geneList"),
                                                                         # plotlyOutput("gene_bmd_plot")
                                                                         fluidRow(
                                                                           column(6,plotOutput("NGVenn")),
                                                                           column(6, DT::dataTableOutput("VennDF"))
                                                                         ),
                                                                         hr(),
                                                                        # "Genes that show dose-response behaviour at all timepoints",
                                                                         fluidRow(column(6,uiOutput("intersectionNameUI"))),
                                                                         fluidRow(
                                                                             column(6,uiOutput("nclustvenn")),
                                                                             column(6,selectInput(inputId = "geneClustMeth", label = "Agglomerative Method",
                                                                                                  choices = list("ward","average","complete","single"), selected = "ward"))
                                                                          ),
                                                                          fluidRow(column(12, plotlyOutput("gene_bmd_plot"))
                                                                          )
                                                         )

                                              ) 
                                            )
                                            
                                   )     
                                   
                          )
                        )
               ),
               tabPanel(value = "enrTab",title = "Enrichment",
                #         fluidRow(
                #           column(6, uiOutput("timePointSel3"))
                #         ),    
                # fluidRow(column(12, DT::dataTableOutput("enrich_table")))         
        
                sidebarLayout(
                  sidebarPanel(
                    wellPanel(
                      tags$h5("1. Data Selection"),
                      fluidRow(                      
                        selectInput(inputId = "level", label = "Browse hierarchy: choose a level",choices = list(1,2,3))
                      ),
                      fluidRow(
                        column(4,uiOutput("chose_lev1")),
                        shinyBS::bsTooltip(id = "chose_lev1",title = "Note: remove ALL from the list for specific selection.",placement = "top"),
                        column(4,uiOutput("chose_lev2")),
                        shinyBS::bsTooltip(id = "chose_lev2",title = "Note: remove ALL from the list for specific selection.",placement = "top"),
                        column(3,uiOutput("chose_lev3")),
                        shinyBS::bsTooltip(id = "chose_lev3",title = "Note: remove ALL from the list for specific selection.",placement = "top")
                        
                      ),
                      fluidRow(
                        uiOutput("selectColumn"),
                        shinyBS::bsTooltip(id = "selectColumn",title = "Note: remove ALL from the list for specific selection.",placement = "top")
                        
                      )
                    ),
                    wellPanel(
                      tags$h5("2. Plot section"),
                      
                      fluidRow(
                        column(4,checkboxInput("doGrouping", "Show categories", value = TRUE)),
                        column(4,checkboxInput("aspectRatio", "Keep aspect ratio", value = TRUE)),
                        column(4,actionButton("do", "Plot Map")),
                        shinyBS::bsTooltip(id = "do",title ="NOTE: press the Plot Mat button every time you update the map!",placement = "bottom")
                        
                      )
                    ),
                    wellPanel(
                      tags$h5("3. Download Selection"),
                      fluidRow(
                        column(4,textInput(inputId ="img_width", value = 15,label = "Width")), #width
                        column(4,textInput(inputId ="img_height", value = 30,label = "Height")),
                        column(4,downloadButton('downloadData')),
                        shinyBS::bsTooltip(id = "downloadData",title ="NOTE: when downloading, specify image size in inches ",placement = "bottom")
                      )
                    )
                    # wellPanel(
                    #   tags$h5("4. Clustering Selection"),
                    #   fluidRow(
                    #     column(4,uiOutput("nClust")),
                    #     column(4,selectInput("ClusterMethod","Select aggregation method",list("ward","complete","single"),selected = "complete")),
                    #     column(4,selectInput("Distance","Select distance",list("euclidean","jaccard","jaccard+euclidean"),selected = "jaccard"))
                    #   ),
                    #   fluidRow(
                    #     column(6,actionButton("doCluster", "Cluster Samples")),
                    #     column(6,actionButton("resetCluster","Reset Cluster"))
                    #   )
                    # )
                    
                  ),
                  mainPanel(
                    tags$h5("Use scrollbars to navigate and see the whole map"),
                    
                    tabsetPanel(
                      
                      tabPanel("Heatmap",fluidRow(column(12,align="left",shinycssloaders::withSpinner(plotOutput(outputId="heatmap"), type=6)))),
                      #tabPanel("Clustering",plotOutput(outputId="hclust_plot", width = "100%")),
                      tabPanel("Cluster Bubble Plot", plotlyOutput("pathway_bubble")),
                      tabPanel("Mean BMD for Time Point", plotlyOutput("meanBMD_timepoint")),
                      tabPanel("Gene BMD in pathway",
                               fluidRow(
                                 column(6, uiOutput("timePointSelPat"))
                               ),
                               fluidRow(
                                 column(12,
                                        DT::dataTableOutput("PAT_table")
                                 )),
                               fluidRow(
                                 column(12,
                                        plotlyOutput("path_bmd_dist"))
                               )
                               ),
                      tabPanel("Pathways Table", 
                               "Here I will display the enriched pathway in tabular form",
                                fluidRow(uiOutput("timePointSelTab")),
                                fluidRow(downloadButton("downloadEnrichedPathwayTables", "Download")),
                                fluidRow(DT::dataTableOutput("PatTable"))
                      ),
                      #tabPanel("Genes Heatmap", "Here I will display the gene heatmap")
                      tabPanel("Heatmap Genes",
                               fluidRow(column(4,
                                               selectInput(inputId="levelGene", label="Choose a hierarchy level", choices=list(1,2,3))
                               ),column(4,
                                        uiOutput("choosePath")
                               ),column(4,
                                        selectInput(inputId="selScoreType", label="Show Values", choices=list("BMD"="lfc", "P-Value"="pval", "Combined"="comb"), selected="lfc")
                               )),fluidRow(column(4,
                                                  #shinyBS::bsButton("doGeneHeatMap", label="Plot", style="danger", icon=icon("exclamation-circle"))
                                                  actionButton("doGeneHeatMap", "Plot")
                               )),fluidRow(column(12,align="center",
                                                  shinycssloaders::withSpinner(plotOutput(outputId="heatmapGenes"), type=6)
                               ))
                      )
                    )
                    #fluidRow(plotOutput(outputId="heatmap",height = 900))
                  )
                )
               
                
               )
        )
      )
    )
  )
)

