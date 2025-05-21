ui <- function(request){
  fluidPage(
    theme = bslib::bs_theme(bootswatch = "united"),
    useShinyjs(),
    tags$head(tags$link(rel = "icon", href = "icon.png")),
    tags$head(tags$style(HTML(".shiny-output-error-validation {color: red; font-weight: bold;}"))),
    list(tags$style(HTML(".navbar-default .navbar-nav { font-weight: bold; font-size: 16px; }"))),
    # navbar_options = navbar_options(
    #   bg = phs_colours(colourname = "phs-purple"), # background navbar colour
    #   collapsible = TRUE # collapse tabs on smaller screens
    # ),
    #   tags$style("
    #   .sidebar {
    #     background-color: #f2f2f2;
    #   }
    # "),
    # imageOutput("headimage",height = "100px"),
    
    ### Page title 
    fluidRow(
      br(),
      # column(3, titlePanel("Data set:")),
      column(3,imageOutput("headimage", height = 90)),
      column(2),
      column(
        4, br(),
        selectInput(inputId = "dn",label = NULL,dataname,width = "300px")
      ),
      column(
        3, br(),
        downloadButton(
          outputId = "SeuratResults",
          label = "Download SeuratObject"
          # class = "btn-lg btn-success"
        )
      )
    ),
    tabsetPanel(
      id = "Tab",
      # navbarPage(
      NULL,
      # Dim Plot-----
      ### cellInfo vs geneExpr on dimRed
      tabPanel(
        value = "DimPlot",
        HTML("CellInfo vs GeneExpr"),
        h4("Cell information vs gene expression on reduced dimensions"),
        "In this tab, users can visualise both cell information and gene ",
        "expression side-by-side on low-dimensional representions.",
        br(),br(),
        layout_sidebar(
          sidebar = sidebar(
            width = 300,
            selectInput("coord", "Coordinates:", choices = coordinates),
            accordion(
              ## Subset-----
              accordion_panel(
                title = "Subset", value = "dimsubset",
                subsetUI(id = "dimPlotSubset")
              ),
              ## graphics controls-----
              accordion_panel(
                title = "Graphics controls",value = "dimctrl",
                sliderInput(
                  "sc1a1siz", "Point size:",
                  min = 0, max = 4, value = 1.25, step = 0.25
                ),
                PlotFontSizeUI(id = "dimPlotPFsize"),
                radioButtons("sc1a1asp", "Aspect ratio:",
                             choices = c("Square", "Fixed", "Free"),
                             selected = "Square", inline = TRUE),
                checkboxInput("sc1a1txt", "Show axis text", value = FALSE)
              )
            )
          ),
          layout_columns(
            card(
              scatterUI("left", whichinfo = "CellInfo", cellLabels = defaultcellmeta)
            ),
            card(
              scatterUI("right", whichinfo = "GeneExpr", cellLabels = defaultcellmeta)
            )
          )
        ),
      ),     # End of tab (2 space)
      
      # Gene Coexpression-----
      tabPanel(
        value = "Coex",
        HTML("Gene coexpression"),
        h4("Coexpression of two genes on reduced dimensions"),
        "In this tab, users can visualise the coexpression of two genes ",
        "on low-dimensional representions.",
        br(),br(),
        layout_sidebar(
          sidebar = sidebar(
            width = 300,
            selectGeneUI(id = "CoExprGene1", label = "Gene 1:"),
            selectGeneUI(id = "CoExprGene2", label = "Gene 2:"),
            selectInput("coexpr_coord", "Coordinates:", choices = coordinates),
            accordion(
              accordion_panel(
                title = "Subset", value = "coexsubset",
                subsetUI(id = "CoExprSubset")
              ),
              accordion_panel(
                title = "Graphics controls",value = "coexctrl",
                sliderInput(
                  "CoExprPointSize", "Point size:",
                  min = 0, max = 4, value = 1.25, step = 0.25
                ),
                PlotFontSizeUI(id = "CoExprPFsize"),
                radioButtons("CoExprAsp", "Aspect ratio:",
                             choices = c("Square", "Fixed", "Free"),
                             selected = "Square", inline = TRUE),
                checkboxInput("CoExprText", "Show axis text", value = FALSE),
                
                radioButtons(
                  "CoExprColor", "Colour:",
                  choices = c("Red (Gene1); Blue (Gene2)",
                              "Orange (Gene1); Blue (Gene2)",
                              "Red (Gene1); Green (Gene2)",
                              "Green (Gene1); Blue (Gene2)"),
                  selected = "Red (Gene1); Blue (Gene2)"
                ),
                radioButtons(
                  "CoExprOrder", "Plot order:",
                  choices = c("Max-1st", "Min-1st", "Original", "Random"),
                  selected = "Max-1st", inline = TRUE
                )
              )
            )
          ),
          layout_columns(
            card(
              uiOutput("CoexPlot.ui"),
              fluidRow(downloadPlotUI("download_coex"))
            ),
            card(
              uiOutput("coexLeg.ui"),
              # downloadPlotUI("download_coex_leg"),
              br(), h4("Cell numbers"),
              card(dataTableOutput("coex.dt"))
            ),
            col_widths = c(7,5)
          )
        )
        
        # fluidRow(
        #   column(
        #     3, # h4("Dimension Reduction"),
        #     fluidRow(
        #       column(
        #         12, selectInput("coexpr_coord", "Coordinates:", choices = coordinates)
        #       )
        #     )
        #   ),
        #   column(3, subsetUI(id = "CoExprSubset")),
        #   ## graph control---
        #   column(
        #     6,
        #     actionButton("CoExprGraphToggle", "Toggle graphics controls"),
        #     conditionalPanel(
        #       condition = "input.CoExprGraphToggle % 2 == 1",
        #       fluidRow(
        #         column(
        #           6,
        #           sliderInput(
        #             "CoExprPointSize", "Point size:",
        #             min = 0, max = 4, value = 1.25, step = 0.25
        #           ),
        #           PlotFontSizeUI(id = "CoExprPFsize")
        #         ),
        #         column(
        #           6, radioButtons("CoExprAsp", "Aspect ratio:",
        #                           choices = c("Square", "Fixed", "Free"),
        #                           selected = "Square", inline = TRUE),
        #           checkboxInput("CoExprText", "Show axis text", value = FALSE)
        #         )
        #       )
        #     )
        #   )
        # ),   # End of fluidRow (4 space)
        ## CoExpr Plot-----
        # sidebarLayout(
        #   sidebarPanel(
        #     width = 3,
        #     # style="border-right: 2px solid black",
        #     h4("Gene Expression"),
        #     selectGeneUI(id = "CoExprGene1", label = "Gene 1:"),
        #     selectGeneUI(id = "CoExprGene2", label = "Gene 2:"),
        #     actionButton("CoExprPlotToggle", "Toggle plot controls"),
        #     conditionalPanel(
        #       condition = "input.CoExprPlotToggle % 2 == 1",
        #       radioButtons("CoExprColor", "Colour:",
        #                    choices = c("Red (Gene1); Blue (Gene2)",
        #                                "Orange (Gene1); Blue (Gene2)",
        #                                "Red (Gene1); Green (Gene2)",
        #                                "Green (Gene1); Blue (Gene2)"),
        #                    selected = "Red (Gene1); Blue (Gene2)"),
        #       radioButtons("CoExprOrder", "Plot order:",
        #                    choices = c("Max-1st", "Min-1st", "Original", "Random"),
        #                    selected = "Max-1st", inline = TRUE)
        #     )
        #   ), # End of column (6 space)
        #   mainPanel(
        #     width = 9,
        #     fluidRow(
        #       column(
        #         8, #style="border-right: 2px solid black",
        #         card(
        #           uiOutput("CoexPlot.ui")
        #         ),
        #         downloadPlotUI("download_coex")
        #       ), # End of column (6 space)
        #       column(
        #         4,
        #         card(
        #           uiOutput("coexLeg.ui"),
        #           # downloadPlotUI("download_coex_leg"),
        #           br(), h4("Cell numbers"),
        #           dataTableOutput("coex.dt")
        #         )
        #       )  # End of column (6 space)
        #     )
        #   )
        # )
      ), # fluidRow

      # violin / box-----
      tabPanel(
        value = "viobox",
        HTML("Violinplot / Boxplot"),
        h4("Cell information / gene expression violin plot / box plot"),
        "In this tab, users can visualise the gene expression or continuous cell information ",
        "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).",
        br(),br(),
        layout_sidebar(
          sidebar = sidebar(
            width = 400,
            hidden(textInput("hiddenVioMeta",label = "hidden vio meta")),
            selectInput("viobox_metaname", "Cell information (X-axis):",
                        choices = defaultcellmeta,
                        selected = defaultcellmeta[2]) |>
              helper(type = "inline", size = "m", fade = TRUE,
                     title = "Cell information to group cells by",
                     content = c("Select categorical cell information to group cells by",
                                 "- Single cells are grouped by this categorical covariate",
                                 "- Plotted as the X-axis of the violin plot / box plot")),
            radioButtons(
              "SelectCellGene",
              label = NULL,
              choices = c("Cell", "Gene"),
              selected = "Cell",
              inline = TRUE
            ),
            tabsetPanel(
              id = "CellGeneUI",
              type = "hidden",
              tabPanel(
                "Cell",
                selectInput("CellNum", "Cell Info (Y-axis):", choices=NULL) |>
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell Info / Gene to plot",
                    content = c(
                      "Select cell info / gene to plot on Y-axis",
                      "- Can be continuous cell information (e.g. nUMIs / scores)",
                      "- Can also be gene expression"
                    )
                  )
              ),
              tabPanel(
                "Gene",
                selectGeneUI(id = "vioGene", label = "Gene name (Y-axis):"),
              )
            ),
            # radioButtons(
            #   "vioboxType", "Plot type:",
            #   choices = c("violin", "boxplot"),
            #   selected = "violin", inline = TRUE
            # ),
            checkboxInput("vioboxpts", "Show data points", value = FALSE),
            accordion(
              accordion_panel(
                title = "Subset", value = "viosubset",
                subsetUI(id = "violinSubset") # ,br(),
              ),
              accordion_panel(
                title = "Graphics controls",value = "vio_graph_ctrl",
                sliderInput(
                  "vioboxPointSize", "Point size:",
                  min = 0, max = 4, value = 1.25, step = 0.25
                ),
                PlotFontSizeUI(id = "vioboxPFsize")
                # actionButton("violinGraphToggle", "Toggle graphics controls"),
                # conditionalPanel(
                #   condition = "input.violinGraphToggle % 2 == 1",
                #   
                # )
              )
            )
          ), # End of column (6 space)
          layout_columns(
            # uiOutput("viobox.ui"),
            navset_pill(
              nav_panel(
                "Violin Plot", br(),
                fluidRow(
                  uiOutput("vio_plot.ui"),
                  downloadPlotUI("download_violin")
                )
              ), 
              nav_panel(
                "Box Plot", br(),
                fluidRow(
                  uiOutput("box_plot.ui"),
                  downloadPlotUI("download_box")
                )
              )
            ),
            
            # downloadPlotUI("download_viobox")
            # fluidRow(downloadPlotUI("download_viobox"))
            
            # downloadPlotUI("download_viobox")
            # tabsetPanel(
            #   id = "vioboxPlotText",
            #   type = "hidden",selected = "plot",
            #   tabPanel("plot", uiOutput("viobox.ui")),
            #   tabPanel("text", HTML("Enter a gene name"))
            # )
          )  # End of column (6 space)
        )    # End of fluidRow (4 space)
      ),     # End of tab (2 space)

      # Proportion-----
      tabPanel(value = "prop",
        HTML("Proportion plot"),
        h4("Proportion / cell numbers across different cell information"),
        "In this tab, users can visualise the composition of single cells based on one discrete ",
        "cell information across another discrete cell information. ",
        "Usage examples include the library or cellcycle composition across clusters.",
        br(),br(),
        layout_sidebar(
          sidebar = sidebar(
            width = 400,
            # 3, style="border-right: 2px solid black",
            selectInput("prop_metaname", "Cell information to plot (X-axis):",
                        choices = defaultcellmeta,
                        selected = defaultcellmeta[2]) |>
              helper(type = "inline", size = "m", fade = TRUE,
                     title = "Cell information to plot cells by",
                     content = c("Select categorical cell information to plot cells by",
                                 "- Plotted as the X-axis of the proportion plot")),
            selectInput("prop_group", "Cell information to group / colour by:",
                        choices = defaultcellmeta,
                        selected = "seurat_clusters") |>
              helper(type = "inline", size = "m", fade = TRUE,
                     title = "Cell information to group / colour cells by",
                     content = c("Select categorical cell information to group / colour cells by",
                                 "- Proportion / cell numbers are shown in different colours")),
            radioButtons("prop_value", "Plot value:",
                         choices = c("Proportion", "CellNumbers"),
                         selected = "Proportion", inline = TRUE),
            checkboxInput("prop_filp", "Flip X/Y", value = FALSE),
            accordion(
              accordion_panel(
                title = "Subset", value = "prosubset",
                subsetUI(id = "propSubset") # ,br(),
              ),
              accordion_panel(
                title = "Graphics controls",value = "proP_graph_ctrl",
                PlotFontSizeUI(id = "propPFsize")
              )
            )#,
            # actionButton("propGraphToggle", "Toggle graphics controls"),
            # conditionalPanel(
            #   condition = "input.propGraphToggle % 2 == 1",
            #   
            # )
          ), # End of column (6 space)
          layout_columns(
            fluidRow(
              uiOutput("prop.ui"),
              downloadPlotUI("download_prop")
            ),
          )  # End of column (6 space)
        )    # End of fluidRow (4 space)
      ),     # End of tab (2 space)

      # Bubble/Heatmap-----
      tabPanel(
        value = "BH",
        HTML("Bubbleplot / Heatmap"),
        h4("Gene expression bubbleplot / heatmap"),
        "In this tab, users can visualise the gene expression patterns of ",
        "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(),
        "The normalised expression are averaged, log-transformed and then plotted.",
        br(),br(),
        layout_sidebar(
          sidebar = sidebar(
            width = 400,
      #   # sidebarLayout(
      #   #   sidebarPanel(
      #       # width = 3,
      #       # style="border-right: 2px solid black",
            textAreaInput("BHgenes", HTML("List of gene names <br />
                                          (Max 50 genes, separated <br />
                                           by , or ; or newline):"),
                          height = "200px") |>
              helper(type = "inline", size = "m", fade = TRUE,
                     title = "List of genes to plot on bubbleplot / heatmap",
                     content = c("Input genes to plot",
                                 "- Maximum 50 genes (due to ploting space limitations)",
                                 "- Genes should be separated by comma, semicolon or newline")),
      h5(htmlOutput("BHtxt")),      
      selectInput("BHgroup", "Group by:",
                        choices = defaultcellmeta,
                        selected = defaultcellmeta[2]) |>
              helper(type = "inline", size = "m", fade = TRUE,
                     title = "Cell information to group cells by",
                     content = c("Select categorical cell information to group cells by",
                                 "- Single cells are grouped by this categorical covariate",
                                 "- Plotted as the X-axis of the bubbleplot / heatmap")),
            # radioButtons("BHplt", "Plot type:",
            #              choices = c("Bubbleplot", "Heatmap"),
            #              selected = "Bubbleplot", inline = TRUE),
            checkboxInput("BHscl", "Scale gene expression", value = TRUE),
            checkboxInput("BHrow", "Cluster rows (genes)", value = TRUE),
            checkboxInput("BHcol", "Cluster columns (samples)", value = FALSE),
            br(),
            accordion(
              accordion_panel(
                title = "Subset", value = "BHsubset",
                subsetUI(id = "BHSubset") # , br(), br(),
              ),
              accordion_panel(
                title = "Graphics controls",value = "BH_graph_ctrl",
                radioButtons(
                  "BHcols",
                  "Colour scheme:",
                  choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                  selected = "Blue-Yellow-Red"
                ),
                radioButtons(
                  "BHpsz",
                  "Plot size:",
                  choices = c("Small", "Medium", "Large"),
                  selected = "Medium",
                  inline = TRUE
                ),
                radioButtons(
                  "BHfsz",
                  "Font size:",
                  choices = c("Small", "Medium", "Large"),
                  selected = "Medium",
                  inline = TRUE
                )
              ),
            ),
            # actionButton("BHGraphToggle", "Toggle graphics controls"),
            # conditionalPanel(
            #   condition = "input.BHGraphToggle % 2 == 1",
            #
            # )
          ), # End of column (6 space)
          layout_columns(
            # 9,
            # card(
            #   h4(htmlOutput("BHtxt")),
            #   uiOutput("BHplot.ui"),
            #   downloadPlotUI("download_BH")
            # )
            navset_pill(
              nav_panel(
                "Bubble Plot",br(),
                fluidRow(
                  # ,br(),
                  uiOutput("Bubplot.ui"),br(),br(),
                  downloadPlotUI("download_bubble")
                )
              ), 
              nav_panel(
                "Heatmap",br(),
                fluidRow(
                  # h4(htmlOutput("BHtxt")),br(),
                  uiOutput("HeatmapPlot.ui"),
                  downloadPlotUI("download_heatmap")
                )
              )
            )
          )  # End of column
        )    # End of layout_sidebar
      ),      # End of tab (2 space)
      # DEG-----
      # nav_panel(
      tabPanel(
        value = "DEG",
        HTML("Differential Expression"),
        h4("Volcano plot for differential gene expression"),
        "This tab enables visualization of differentially expressed genes (log2 fold change) across cell clusters.",
        br(),br(),
        layout_sidebar(
          sidebar = sidebar(
            width = 400,
            selectInput("deg_metaname", "Metaname:",choices = NULL),
            selectInput("deg_level", "Cluster:", choices = NULL),
            fluidRow(
              column(
                6,
                numericInput(
                  "p_adj", "Adjust p value:",
                  value = 0.05, min = 0, max = 0.05, step = 0.01
                )
              ),
              column(
                6,
                numericInput(
                  "log2fc", "Log2 Fold Change:",
                  value = 1, min = 0.2, max = 10, step = 0.1
                )
              )
            ),
            accordion(
              accordion_panel(
                title = "Control 2",value = "control 2",
                sliderInput(
                  "deg_label_num", "labeled cells Number",
                  value = 5, min = 0, max = 10
                )
                
              )
            )
          ),
          layout_columns(
            navset_pill(
              nav_panel(
                "Volcano", br(),
                fluidRow(
                  uiOutput("volcano_plot_ui"),br(),
                  downloadPlotUI("download_volcano")
                )
              ), 
              nav_panel(
                "FC-Pct", br(),
                fluidRow(
                  uiOutput("fc_pct_plot_ui"),
                  downloadPlotUI("download_fcpct")
                )
              ),
              nav_panel("Table", br(), dataTableOutput("degTable"))
            )
          )
        )
      ) # End of DEG
    )
  )
}
