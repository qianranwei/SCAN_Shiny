### Start server code 
server <- function(input, output, session) {
  observe_helpers()
  output$headimage <- renderImage({
    list(
      src = "www/webicon.png",
      # contentType = "image/jpeg",
      width = "100%",
      height = 80
    )
  }, deleteFile = FALSE)
  
  # output$webhead <- renderImage({
  #   list(
  #     src = file.path("www/webhead.png"),
  #     contentType = "image/jpeg",
  #     width = 1400
  #     # height = 120
  #   )
  # }, deleteFile = FALSE)
  # observe({print(file.path(dn(), "sc1conf.rds"))})
  # observe({print(getwd())})
  # observe({message("dn: ",input$dn)})
  # read data-----
  dn <- input$dn |> reactive()
  conf <- readRDS(file.path(dn(), "sc1conf.rds")) |> reactive()
  meta <- readRDS(file.path(dn(), "sc1meta.rds")) |> reactive()
  gene_loc <- readRDS(file.path(dn(), "sc1gene.rds")) |> reactive()
  def  <- readRDS(file.path(dn(), "sc1def.rds")) |> reactive()
  deg  <- readRDS(file.path(dn(), "sc1deg.rds")) |> reactive()
  # sc1info  <- readRDS(file.path(input$dn,"sc1info.rds"))|>reactive()
  # observe({message("dn: ",dn())})
  # observe({print(head(conf()))})
  
  cell_metanames <- conf()[type == "label"]$ID |> reactive() # classified variable
  label_levels <- conf()[match(cell_metanames(),ID)]$fID |>
    strsplit(split = "\\|") |> setNames(cell_metanames()) |> reactive()
  cell_num <- conf()[type == "num"]$ID |> reactive() # numberic variable
  genenames <- reactive(names(gene_loc()))
  h5path <- file.path(dn(),"sc1gexpr.h5") |> reactive()
  
  # DimPlot-----
  ## update subset-----
  dimPlotSubsetSelected <- subsetServer(id = "dimPlotSubset", metanames = cell_metanames, meta_levels = label_levels)
  dimPlotPFsize <- PlotFontSizeServer(id = "dimPlotPFsize")
  
  ## Plot-----
  dim_coord <- reactive(input$coord)
  dimPlotPointSize <- reactive(input$sc1a1siz)
  dimPlotAspectRatio <- reactive(input$sc1a1asp)
  dimPlotAxisText <- reactive(input$sc1a1txt)
  
  scatterserver(
    "left",
    dn = dn,
    cellLabel = cell_metanames,
    gene_loc = gene_loc,
    genenames = genenames,
    conf = conf,
    meta = meta,
    coord = dim_coord,
    metaname_sub = dimPlotSubsetSelected$metaname, # sub_name,
    levels_sub = dimPlotSubsetSelected$cellLevels, # sub_levels,
    height = dimPlotPFsize$plot_size,
    point_size = dimPlotPointSize,
    font_size = dimPlotPFsize$font_size,
    aspect_ratio = dimPlotAspectRatio,
    axis_text = dimPlotAxisText
  )
  scatterserver(
    "right",
    dn = dn,
    cellLabel = cell_metanames,
    gene_loc = gene_loc,
    genenames = genenames,
    conf = conf,
    meta = meta,
    coord = dim_coord,
    metaname_sub = dimPlotSubsetSelected$metaname, # 用于取子集的meta列名
    levels_sub = dimPlotSubsetSelected$cellLevels, # 取子集的levels
    height = dimPlotPFsize$plot_size,
    point_size = dimPlotPointSize,
    font_size = dimPlotPFsize$font_size,
    aspect_ratio = dimPlotAspectRatio,
    axis_text = dimPlotAxisText
  )
  
  # Gene Coexpr-----
  CoExprSubsetSelected <- subsetServer(id = "CoExprSubset", metanames = cell_metanames, meta_levels = label_levels)
  CoExprMetaName <- CoExprSubsetSelected$metaname
  CoExprLevels <- CoExprSubsetSelected$cellLevels

  CoExprPFsize <- PlotFontSizeServer(id = "CoExprPFsize")
  CoexPlotSize <- CoExprPFsize$plot_size
  CoexFontSize <- CoExprPFsize$font_size
  # observe({print(paste0("plot_size: ",CoexPlotSize()))})

  CoExprGene1 <- selectGeneServer(id = "CoExprGene1", genenames = genenames)
  CoExprGene2 <- selectGeneServer(id = "CoExprGene2", genenames = genenames)
  # observe({cat(paste0("Gene1: \n",CoExprGene1()|>is.null(),"\n"))})
  # observe({cat(paste0("Gene2: \n",CoExprGene2()|>length(),"\n"))})

  coex_coord_1 <- paste0(input$coexpr_coord,"1") |> reactive()
  coex_coord_2 <- paste0(input$coexpr_coord,"2") |> reactive()

  ## DRcoex-----
  Coex <- scDRcoex(
    inpConf = conf(),
    inpMeta = meta(),
    inpdrX = coex_coord_1(),
    inpdrY = coex_coord_2(),
    inp1 = CoExprGene1(),
    inp2 = CoExprGene2(),
    inpsub1 = CoExprMetaName(),
    inpsub2 = CoExprLevels(),
    inpH5 = h5path(),
    inpGene = gene_loc(),
    inpsiz = input$CoExprPointSize,
    inpcol = input$CoExprColor,
    inpord = input$CoExprOrder,
    inpfsz = CoexFontSize(),
    inpasp = input$CoExprAsp,
    inptxt = input$CoExprText
  ) |> reactive()
  output$CoexPlot <- renderPlot({Coex()})
  output$CoexPlot.ui <- renderUI({
    plotOutput("CoexPlot", height = CoexPlotSize()) # CoExprPFsize$plot_size
  })
  ## download Coex-----
  coexPlotName <- paste(
    input$coexpr_coord,CoExprGene1(),CoExprGene2(),sep = "_"
  ) |> reactive()
  downloadPlotServer(
    id = "download_coex", plotname = coexPlotName, plot = Coex
  )
  ## Leg-----
  CoexLeg <- scDRcoexLeg(
    CoExprGene1(),
    CoExprGene2(),
    input$CoExprColor,
    CoexFontSize()
  ) |> reactive()
  output$coexLeg <- renderPlot({CoexLeg()})
  output$coexLeg.ui <- renderUI({
    plotOutput("coexLeg", height = "300px")
  })
  ## download Leg-----
  coexPlotName <- paste(
    "Coex_Leg",CoExprGene1(),CoExprGene2(),sep = "_"
  ) |> reactive()
  downloadPlotServer(
    id = "download_coex_leg", plotname = coexPlotName, plot = CoexLeg
  )
  ## Table-----
  output$coex.dt <- renderDataTable({
    ggData = scDRcoexNum(
      conf(), meta(), CoExprGene1(), CoExprGene2(),
      CoExprMetaName(), CoExprLevels(), h5path(), gene_loc()
    )
    datatable(
      ggData, rownames = FALSE, extensions = "Buttons",
      options = list(
        pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel")
      )
    ) |> formatRound(columns = c("percent"), digits = 2)
  })
  
  # Violin/Boxplot-----
  ## update cell metanames-----
  observeEvent(cell_metanames(),{
    updateSelectInput(
      session, "viobox_metaname", choices = cell_metanames(), selected = NULL) # cell_metanames()[1]
  })
  
  observeEvent(input$hiddenVioMeta,{
    updateSelectInput(
      session, "viobox_metaname", choices = cell_metanames(),
      selected = ifelse(input$hiddenVioMeta == "","orig.ident",input$hiddenVioMeta)
    ) # cell_metanames()[1]
  })

  observeEvent(input$SelectCellGene,{
    updateTabsetPanel(session, inputId = "CellGeneUI",selected = input$SelectCellGene)
  })
  observeEvent(cell_num(),{
    updateSelectInput(
      session, "CellNum", choices = cell_num(), selected = cell_num()[1])
  })
  vioboxGene <- selectGeneServer(id = "vioGene", genenames = genenames)
  vioboxSubsetSelected <- subsetServer(id = "violinSubset", metanames = cell_metanames, meta_levels = label_levels)
  vioboxMetaName <- vioboxSubsetSelected$metaname
  vioboxLevels <- vioboxSubsetSelected$cellLevels

  vioboxPFsize <- PlotFontSizeServer(id = "vioboxPFsize")
  vioboxPlotSize <- vioboxPFsize$plot_size
  vioboxFontSize <- vioboxPFsize$font_size

  vioboxY <- switch(
    input$SelectCellGene,
    Cell = input$CellNum,
    Gene = vioboxGene()
  ) |> reactive()
  
  # violin plot
  VioPlot <- scVioBox(
    conf(),
    meta(),
    input$viobox_metaname,
    vioboxY(),
    vioboxMetaName(),
    vioboxLevels(),
    h5path(),
    gene_loc(),
    "violin",
    input$vioboxpts,
    input$vioboxPointSize,
    vioboxFontSize()
  ) |> reactive()
  output$vio_plot <- renderPlot({VioPlot()})
  output$vio_plot.ui <- renderUI({
    plotOutput("vio_plot", height = vioboxPlotSize())
  })
  # download violin
  violinPlotName <- paste(
    "violin",input$viobox_metaname,vioboxY(),sep = "_"
  ) |> reactive()
  downloadPlotServer(
    id = "download_violin", plotname = violinPlotName, plot = VioPlot
  )
  
  # box plot
  BoxPlot <- scVioBox(
    conf(),
    meta(),
    input$viobox_metaname,
    vioboxY(),
    vioboxMetaName(),
    vioboxLevels(),
    h5path(),
    gene_loc(),
    "boxplot",
    input$vioboxpts,
    input$vioboxPointSize,
    vioboxFontSize()
  ) |> reactive()
  output$box_plot <- renderPlot({BoxPlot()})
  output$box_plot.ui <- renderUI({
    plotOutput("box_plot", height = vioboxPlotSize())
  })
  # download box
  boxPlotName <- paste(
    "Box",input$viobox_metaname,vioboxY(),sep = "_"
  ) |> reactive()
  downloadPlotServer(
    id = "download_box", plotname = boxPlotName, plot = BoxPlot
  )
  
  
  # observeEvent(vioboxY(),{
  #   if(is.null(vioboxY())) {
  #     vioboxPT <- reactive("text")
  #   } else {
  #     vioboxPT <- reactive("plot")
  #   }
  #   updateTabsetPanel(session,inputId = "vioboxPlotText",selected = vioboxPT())
  # })
  ## download-----
  # vioboxPlotName <- paste(
  #   input$vioboxType,input$viobox_metaname,vioboxY(),sep = "_"
  # ) |> reactive()
  # downloadPlotServer(
  #   id = "download_viobox", plotname = vioboxPlotName, plot = VioBox
  # )
  
  # proportion-----
  ## select metaname-----
  observeEvent(cell_metanames(),{
    updateSelectInput(
      session, "prop_metaname", choices = cell_metanames(), selected = cell_metanames()[1])
  })
  observeEvent(cell_metanames(),{
    updateSelectInput(
      session, "prop_group", choices = cell_metanames(), selected = "seurat_clusters")
  })
  propSubsetSelected <- subsetServer(id = "propSubset", metanames = cell_metanames, meta_levels = label_levels)
  propMetaName <- propSubsetSelected$metaname
  propLevels <- propSubsetSelected$cellLevels
  
  propPFsize <- PlotFontSizeServer(id = "propPFsize")
  propPlotSize <- propPFsize$plot_size
  propFontSize <- propPFsize$font_size
  ## plot-----
  Prop <- scProp(
    conf(),
    meta(),
    input$prop_metaname,
    input$prop_group,
    propMetaName(),
    propLevels(),
    input$prop_value,
    input$prop_filp,
    propFontSize()
  ) |> reactive()
  output$prop <- renderPlot({Prop()}) 
  output$prop.ui <- renderUI({ 
    plotOutput("prop", height = propPlotSize()) 
  })
  ## download-----
  # propPlotName <- paste(
  #   input$prop_value,input$prop_metaname,input$prop_group,sep = "_"
  # ) |> reactive()
  # downloadPlotServer(
  #   id = "download_prop", plotname = propPlotName, plot = Prop
  # )
  
  # Bubble-----
  observeEvent(def(), {
    updateTextAreaInput(
      session,
      inputId = "BHgenes",
      value = paste0(def()$genes, collapse = ", ")
    )
  })
  observeEvent(cell_metanames(),{
    updateSelectInput(
      session, "BHgroup", choices = cell_metanames(), selected = cell_metanames()[1])
  })
  BHSubsetSelected <- subsetServer(id = "BHSubset", metanames = cell_metanames, meta_levels = label_levels)
  BHMetaName <- BHSubsetSelected$metaname
  BHLevels <- BHSubsetSelected$cellLevels
  ## plot-----

  output$BHtxt <- renderUI({
    geneList <- scGeneList(input$BHgenes, gene_loc()) |> reactive()
    # observe(print(input$BHgenes))
    if (nrow(geneList()) > 50) {
      HTML("More than 50 input genes! Please reduce the gene list!")
    } else {
      oup <- ifelse(
        test = nrow(geneList()[present == FALSE]) == 0,
        yes = paste0(nrow(geneList()[present == TRUE]), " genes OK and will be plotted"),
        no = paste0(
          paste0(nrow(geneList()[present == TRUE]), " genes OK and will be plotted"),
          "<br/>",
          nrow(geneList()[present == FALSE]),
          " genes not found (",
          paste0(geneList()[present == FALSE]$gene, collapse = ", "),
          ")"
        )
      ) |> reactive()
      HTML(oup())
    }
  })
  # Bubble plot
  Bubplot <- scBubbHeat(
    conf(),
    meta(),
    input$BHgenes,
    input$BHgroup,
    "Bubbleplot",
    BHMetaName(),
    BHLevels(),
    h5path(),
    gene_loc(),
    input$BHscl,
    input$BHrow,
    input$BHcol,
    input$BHcols,
    input$BHfsz
  ) |> reactive()
  output$Bubplot <- renderPlot({Bubplot()})
  output$Bubplot.ui <- renderUI({
    plotOutput(
      "Bubplot"#, height = pList3[input$BHpsz]
    )
  })
  ## bubble download
  BubPlotName <- paste(
    "Bubble",input$BHgroup,sep = "_"
  ) |> reactive()
  downloadPlotServer(
    id = "download_bubble", plotname = BubPlotName, plot = Bubplot
  )
  
  # heatmap plot
  HeatmapPlot <- scBubbHeat(
    conf(),
    meta(),
    input$BHgenes,
    input$BHgroup,
    "Heatmap",
    BHMetaName(),
    BHLevels(),
    h5path(),
    gene_loc(),
    input$BHscl,
    input$BHrow,
    input$BHcol,
    input$BHcols,
    input$BHfsz
  ) |> reactive()
  output$HeatmapPlot <- renderPlot({HeatmapPlot()})
  output$HeatmapPlot.ui <- renderUI({
    plotOutput(
      "HeatmapPlot"#, height = pList3[input$BHpsz]
    )
  })
  ## download-----
  HeatmapPlotName <- paste(
    "Heatmap",input$BHgroup,sep = "_"
  ) |> reactive()
  downloadPlotServer(
    id = "download_heatmap", plotname = HeatmapPlotName, plot = HeatmapPlot
  )
  
  # DEG-----
  deg_metanames <- names(deg()) |> reactive()
  deg_levels <- deg()[[input$deg_metaname]]$cluster |> levels() |> reactive()
  observeEvent(cell_metanames(),{
    updateSelectInput(session,"deg_metaname", "Metaname:",choices = deg_metanames())
  })
  observeEvent(deg_levels(),{
    updateSelectInput(
      session,
      "deg_level",
      "Cluster:",
      choices = deg_levels()
    )
  })

  deg_table <- deg()[[input$deg_metaname]] |> data.table() |> reactive()
  deg_cluster_table <- deg_table()[cluster == input$deg_level] |> reactive()
  deg_filter_table <- filterDegTable(
    df = deg_cluster_table(), p = input$p_adj, fc = input$log2fc
  ) |> reactive()
  
  volcano_plot <- plot_volcano(
    data = deg_cluster_table(),
    FC = input$log2fc,
    PValue = input$p_adj,
    label_num = input$deg_label_num,
    volcano_title = NULL
  ) |> reactive()
  output$volcano_plot <- renderPlot({volcano_plot()})
  output$volcano_plot_ui <- renderUI({
    plotOutput("volcano_plot",height = "500px",width = "800px") # , height = pList3[input$BHpsz]
  })
  ## download
  volcanoPlotName <- paste(
    "Volcano",input$deg_metaname,input$deg_level,sep = "_"
  ) |> reactive()
  downloadPlotServer(
    id = "download_volcano", plotname = volcanoPlotName, plot = volcano_plot
  )
  
  fc_pct_plot <- plot_fc_pctdiff(
    data = deg_cluster_table(),
    FC = input$log2fc,
    PValue = input$p_adj,
    top_n_diff = input$deg_label_num,
    diff_title = NULL
  ) |> reactive()
  output$fc_pct_plot <- renderPlot({fc_pct_plot()})
  output$fc_pct_plot_ui <- renderUI({
    plotOutput("fc_pct_plot",height = "500px",width = "800px") # , height = pList3[input$BHpsz]
  })
  ## download fc_pct
  fcpctPlotName <- paste(
    "Volcano",input$deg_metaname,input$deg_level,sep = "_"
  ) |> reactive()
  downloadPlotServer(
    id = "download_fcpct", plotname = fcpctPlotName, plot = fc_pct_plot
  )

  output$degTable <- renderDataTable(
    datatable(
      rownames = FALSE,
      options = list(
        columnDefs = list(
          list(
            className = "dt-center",
            targets = "_all"
          )
        )
      ), {deg_filter_table()[,c("gene","cluster","pct.1","pct.2","avg_log2FC","p_val_adj")]}
    ) |>
      formatSignif("p_val_adj", 3) |>
      formatRound(c("pct.1","pct.2","avg_log2FC"),digits = 3) |>
      formatStyle(
        c("gene","cluster","pct.1","pct.2","avg_log2FC","p_val_adj"),
        textAlign = "center"
      )
  )
}
