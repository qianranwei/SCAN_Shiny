# Subset------------------------------------------------------------------------
subsetUI <- function(id) {
  tagList(
    selectInput(
      NS(id,"metaname"), "Cell information to subset:", 
      choices = defaultcellmeta, 
      selected = NULL # defaultcellmeta[1]
    ),
    checkboxGroupInput(
      NS(id,"cellLevels"), "Select which cells to show", inline = TRUE
    ),
    layout_columns(
      actionButton(NS(id, "allButton"), "Select all groups", class = "btn btn-primary"), 
      actionButton(NS(id, "nonButton"), "Deselect all groups", class = "btn btn-primary")
    )
  )
}

subsetServer <- function(id, metanames, meta_levels) {
  moduleServer(id, function(input, output, session) {
    # toggle button-----
    # observe({print(paste0("actionButton: ",input$ToggleButton))})
    
    # toggle_status <- ifelse(button.status(), "close", "open") |> reactive()
    # observeEvent(toggle_status(),{
    #   updateTabsetPanel(session, inputId = "SubsetPanel",selected = toggle_status())
    # })
    
    # metanames <- conf()[type == "label"]$ID |> reactive() # 筛选出因子型metaname用于取子集
    # observe({print(paste0("metanames: ", paste0(metanames())))})
    # ns <- session$ns  # 关键：获取模块的命名空间函数
    # select metaname-----
    observeEvent(metanames(),{
      updateSelectInput(
        session, "metaname", choices = metanames(), selected = metanames()[1])
    })
    # observe({print(paste0("input$metaname: ",input$metaname))})
    
    # select cell levels-----
    cellLevels <- meta_levels()[[input$metaname]] |> reactive()
    # observe({print(paste0("cellLevels: ",paste0(cellLevels())))})
    # observe({print(paste0("input$cellLevels: ",input$cellLevels))})
    observeEvent(cellLevels(), {
      updateCheckboxGroupInput(
        session, inputId = "cellLevels",
        label = "Select which cells to show",
        choices = cellLevels(), selected = cellLevels(), inline = TRUE
      )
    })
    observeEvent(input$nonButton, {
      updateCheckboxGroupInput(
        session, inputId = "cellLevels",
        label = "Select which cells to show",
        choices = cellLevels(), selected = NULL, inline = TRUE
      )
    })
    observeEvent(input$allButton, {
      updateCheckboxGroupInput(
        session, inputId = "cellLevels",
        label = "Select which cells to show",
        choices = cellLevels(), selected = cellLevels(), inline = TRUE
      )
    })
    
    
    # observe({message("input$metaname: ",input$metaname)})
    # observe({message("input$cellLevels: ",input$cellLevels)})
    # observe({message("selected.metaname: ",selected.metaname())})
    # observe({message("selected.levels: ",selected.levels())})
    
    list(
      metaname = input$metaname |> reactive(),
      cellLevels = input$cellLevels |> reactive()
    )
  })
}

# select gene-----
selectGeneUI <- function(id, label = "Gene name:") {
  tagList(
    selectizeInput(
      NS(id,"genename"), label = label, choices = NULL,
      multiple = TRUE,
      options = list(
        placeholder = 'Enter a gene name',
        # maxOptions = 6,
        maxItems = 2,
        closeAfterSelect = TRUE,
        openOnFocus = FALSE,    # 禁止点击输入框时自动展开下拉框
        create = FALSE, # 禁止用户新增选项
        persist = TRUE # 清空已选项后不保留历史
      )
    ) |>
      helper(
        type = "inline",
        size = "m",
        fade = TRUE,
        title = "Gene expression to colour cells by",
        content = c(
          "Select gene to colour cells by gene expression",
          paste0(
            "- Gene expression are coloured in a ",
            "White-Red colour scheme which can be ",
            "changed in the plot controls"
          )
        )
      ),
    hidden(textInput(NS(id,"hiddenGeneName"),label = "hidden gene name")) # 用于实现selectizeInput的bookmarking功能
  )
}

selectGeneServer <- function(id, genenames) {
  moduleServer(
    id, function(input, output, session) {
      observeEvent(
        genenames(),{
          updateSelectizeInput(
            session,
            "genename",
            choices = genenames(),
            selected = NULL,
            server = TRUE
          ) # updateSelectizeInput
        }
      ) # oberveEvent
      # observe({print(input$genename)})
      observeEvent(
        input$hiddenGeneName, {
          updateSelectizeInput(
            session,
            "genename",
            choices = genenames(),
            selected = input$hiddenGeneName,
            server = TRUE
          )
        }
      )
      reactive(input$genename[1])
    }
  ) # moduleServer
}

# download plot-----
downloadPlotUI <- function(id) {
  tagList(
    fluidRow(
    # layout_columns(
      column(width = 3,downloadButton(NS(id,"pdf"), "Download PDF")),
      column(width = 3, downloadButton(NS(id,"png"), "Download PNG")),
      column(
        width = 3,
        numericInput(
          NS(id, "h"),
          "PDF / PNG height:",
          width = "138px",
          min = 4,
          max = 20,
          value = 6,
          step = 0.5
        )
      ),
      column(
        width = 3,
        numericInput(
          NS(id, "w"), "PDF / PNG width:",
          width = "138px",
          min = 4, max = 20, value = 8, step = 0.5
        )
      )
      # div(
      #   style = "display:inline-block",
      #   numericInput(
      #     NS(id, "h"), "PDF / PNG height:",
      #     width = "138px",
      #     min = 4, max = 20, value = 6, step = 0.5
      #   )
      # ),
      # div(
      #   style = "display:inline-block",
      #   numericInput(
      #     NS(id, "w"), "PDF / PNG width:",
      #     width = "138px",
      #     min = 4, max = 20, value = 8, step = 0.5
      #   )
      # ),
      # col_widths = c(3, 3, 3, 3)
    )
    
    
  )
}

downloadPlotServer <- function(id, plotname, plot) {
  moduleServer(id, function(input, output, session) {
    # download cell pdf-----
    output$pdf <- downloadHandler(
      filename = function() {
        paste0(plotname(), ".pdf")
      },
      content = function(file) {
        ggsave(
          file, device = "pdf",
          height = input$h, width = input$w,
          useDingbats = FALSE,
          plot = plot()
        )
      }
    )
    # download png-----
    output$png <- downloadHandler(
      filename = function() {
        paste0(plotname(), ".png")
      },
      content = function(file) {
        ggsave(
          file, device = "png",
          height = input$h, width = input$w,
          plot = plot()
        )
      }
    )
  })
}

# Dim Plot----------------------------------------------------------------------
scatterUI <- function(id, whichinfo, cellLabels = NULL) {
  tagList(
    # CellInfo-----
    navset_pill(
      id = NS(id,"cellgene"),
      selected = whichinfo,
      # type = "hidden",
      nav_panel(
        "CellInfo",
        br(),
        fluidRow(
          column(
            6,
            # 把除UMAP和tSNE以外的metanames放在这里
            selectInput(NS(id,"cellinfo"), label = "Cell information:",
                           choices = cellLabels) |>
              helper(
                type = "inline",
                size = "m",
                fade = TRUE,
                title = "Cell information to colour cells by",
                content = c(
                  "Select cell information to colour cells",
                  "- Categorical covariates have a fixed colour palette",
                  paste0(
                    "- Continuous covariates are coloured in a ",
                    "Blue-Yellow-Red colour scheme, which can be ",
                    "changed in the plot controls"
                  ),
                  "",
                  "The meaning of the variable name",
                  "- orig.ident: The initial cell label",
                  "- seurat_clusters: The cell classification labels from Seurat",
                  "- singler_by_cell: Use SingleR to annotate cell types for each cell individually",
                  "- singler_by_cluster_main: Use SingleR to annotate each cell cluster as the main cell type",
                  "- singler_by_cluster_fine: Use SingleR to annotate each cell cluster as the fine cell type"
                )
              )
          ),
          column(
            6,
            # actionButton(NS(id, "sc1a1tog1"),  # 确保按钮ID也使用命名空间
            actionButton(NS(id, "sc1a1tog1"), "Toggle plot controls"),
            conditionalPanel( 
              condition = paste0("input['", NS(id, "sc1a1tog1"), "'] % 2 == 1"),
              radioButtons(NS(id,"col1"), "Colour (Continuous data):",
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"),
                           selected = "Blue-Yellow-Red"),
              radioButtons(NS(id,"ord1"), "Plot order:",
                           choices = c("Max-1st", "Min-1st", "Original", "Random"),
                           selected = "Original", inline = TRUE),
              checkboxInput(NS(id,"lab"), "Show cell info labels", value = TRUE)
            )
          )
        ),
        fluidRow(column(12, uiOutput(NS(id,"DimPlot.ui")))),
        # plotOutput(NS(id,"DimPlot") , click = "plot_click",width = "600px",height = "700px")
        # plotOutput(NS(id,"testplot") , click = "plot_click")
        
        downloadPlotUI(NS(id,"download_cell"))
        
        # actionButton("sc1a1tog9", "Toggle to show cell numbers / statistics"),
        # conditionalPanel(
        #   condition = "input.sc1a1tog9 % 2 == 1",
        #   h4("Cell numbers / statistics"),
        #   radioButtons("sc1a1splt", "Split continuous cell info into:",
        #                choices = c("Quartile", "Decile"),
        #                selected = "Decile", inline = TRUE),
        #   dataTableOutput("sc1a1.dt")
        # )
      ),
      # GeneExpr-----
      nav_panel(
        title = "GeneExpr",
        br(),
        fluidRow(
          column(
            6, selectGeneUI(id = NS(id,"GeneExpr"))
          ), # end of column
          column(
            6, actionButton(NS(id,"sc1a1tog2"), "Toggle plot controls"),
            conditionalPanel(
              condition = paste0("input['", NS(id, "sc1a1tog2"), "'] % 2 == 1"),
              radioButtons(NS(id,"col2"), "Colour:",
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"),
                           selected = "White-Red"),
              radioButtons(NS(id,"ord2"), "Plot order:",
                           choices = c("Max-1st", "Min-1st", "Original", "Random"),
                           selected = "Max-1st", inline = TRUE)
            )
          )
        ),
        # textOutput(NS(id,"genes"))
        uiOutput(NS(id,"DimPlotGene.ui")),
        downloadPlotUI(NS(id,"download_gene"))
      )
    )
  )
}

scatterserver <- function(id, dn, cellLabel, gene_loc, genenames, conf, meta, coord,
                          metaname_sub, levels_sub, height,
                          point_size, font_size, aspect_ratio, axis_text) {
  moduleServer(
    id,
    function(input, output, session) {
      # update cell metanames-----
      observeEvent(cellLabel(),{
        updateSelectInput(
          session, "cellinfo", choices = cellLabel(), selected = "seurat_clusters") # cellLabel()[1]
      })
      # dim plot CellInfo-----
      cell_mataname <- input$cellinfo |> reactive()
      # observe({print(paste0("subsetname: ",paste0(metaname_sub())))})
      # observe({print(paste0("subsetlevels: ",paste0(levels_sub())))})
      DRcell <- scDRcell(
        inpConf = conf(),
        inpMeta = meta(),
        inpdrX = paste0(coord(),"1"),
        inpdrY = paste0(coord(),"2"), # "UMAP2",
        inp1 = cell_mataname(), #
        inpsub1 = metaname_sub(),
        inpsub2 = levels_sub(),
        inpsiz = point_size(),
        inpcol = input$col1, # "Blue-Yellow-Red",
        inpord = input$ord1,
        inpfsz = font_size(),
        inpasp = aspect_ratio(),
        inptxt = axis_text(),
        inplab = input$lab
      ) |> reactive()
      output$DimPlot <- renderPlot({DRcell()})
      ns <- session$ns  # 关键：获取模块的命名空间函数
      output$DimPlot.ui <- renderUI({ 
        plotOutput(ns("DimPlot"), height = height()) 
      })
      cellUMAPname <- paste(coord(),cell_mataname(),sep = "_") |> reactive()
      downloadPlotServer(
        id = "download_cell", plotname = cellUMAPname, plot = DRcell
      )
      
      # output table
      # output$sc1a1.dt <- renderDataTable({
      #   ggData = scDRnum(
      #     conf(), meta(), input$cellinfo, input$sc1a1inp2,
      #     input$sc1a1sub1, input$sc1a1sub2,
      #     "sc1gexpr.h5", sc1gene, input$sc1a1splt
      #   )
      #   datatable(ggData, rownames = FALSE, extensions = "Buttons",
      #             options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>%
      #     formatRound(columns = c("pctExpress"), digits = 2)
      # })
      # output$testplot <- renderPlot({
      #   hist(meta()$nCount_RNA)
      # })
      
      # dim plot GeneExpr-----
      gene_name <- selectGeneServer(id = "GeneExpr", genenames = genenames)
      h5path <- file.path(dn(),"sc1gexpr.h5") |> reactive()
      DRgene <- scDRgene(
        inpConf = conf(),
        inpMeta = meta(),
        inpdrX = paste0(coord(),"1"),
        inpdrY = paste0(coord(),"2"),
        inp1 = gene_name(),
        inpsub1 = metaname_sub(),
        inpsub2 = levels_sub(),
        inpH5 = h5path(),
        inpGene = gene_loc(),
        inpsiz = point_size(),
        inpcol = input$col2, # "Blue-Yellow-Red",
        inpord = input$ord2,
        inpfsz = font_size(),
        inpasp = aspect_ratio(),
        inptxt = axis_text()
      ) |> reactive()
      output$DimPlotGene <- renderPlot({DRgene()}) 
      output$DimPlotGene.ui <- renderUI({ 
        plotOutput(ns("DimPlotGene"), height = height())
      })
      geneUMAPname <- paste(coord(), gene_name(), sep = "_") |> reactive()
      downloadPlotServer(
        id = "download_gene", plotname = geneUMAPname, plot = DRgene
      )
    }
  )
}

# Plot Font size-----
PlotFontSizeUI <- function(id) {
  tagList(
    radioButtons(
      NS(id,"plot_size"),
      "Plot size:",
      choices = c("Small", "Medium", "Large"),
      selected = "Medium",
      inline = TRUE
    ),
    radioButtons(
      NS(id,"font_size"),
      "Font size:",
      choices = c("Small", "Medium", "Large"),
      selected = "Medium",
      inline = TRUE
    )
  )
}

PlotFontSizeServer <- function(id) {
  moduleServer(
    id, function(input, output, session) {
      # observe({print(paste0("plot_size: ",pList2[input$plot_size]))})
      # observe({print(paste0("plot_size: ",input$font_size))})
      list(
        plot_size = reactive(pList2[input$plot_size]),
        font_size = reactive(input$font_size)
      )
    }
  )
}