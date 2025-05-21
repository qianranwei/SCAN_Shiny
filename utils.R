# functions=====================================================================
# cellcluster <- lapply(
#   sc1deg,
#   function(x) unique(as.character(x$cluster))
# )
# 
# tidy conf & meta
glimpse_meta <- function(data_dir) {
  conf <- readRDS(file.path(data_dir,"sc1conf.rds"))
  meta <- readRDS(file.path(data_dir,"sc1meta.rds"))
  
  res_id <- grep("RNA_snn_res",colnames(meta))
  
  dplyr::glimpse(conf)
  dplyr::glimpse(meta[-res_id,])
}

tidyconf <- function(x, num_variable = NULL) {
  if(is.null(num_variable)) {
    num_variable <- c("nCount_RNA","nFeature_RNA","percent.mito","percent.ribo")
  }
  
  dim_id <- which(x$dimred) # 维度坐标
  num_id <- which(x$ID %in% num_variable) # 连续变量
  lab_id <- setdiff(seq(nrow(x)),union(dim_id,num_id)) # 离散变量
  
  res_id <- grep("RNA_snn_res",x$ID)
  main_id <- setdiff(lab_id,res_id) # main + res = label
  
  # type variable-----
  type <- rep(NA,nrow(x))
  type[dim_id] <- "dimred"
  type[num_id] <- "num"
  type[lab_id] <- "label"
  type <- factor(type)
  
  x <- data.table(x,type)
  
  # sort-----
  x <- x[c(main_id,num_id,dim_id),] # res_id,
  return(x)
}

tidymeta <- function(m,c) {
  label_names <- c[type == "label"]$ID # 筛选标签变量名
  for(l in label_names) {
    if(!is.factor(m[[l]])) m[[l]] <- factor(m[[l]]) # 转换为因子型
  }
  return(m)
}

adjust_conf_files <- function(data_dir, num_variable = NULL) {
  conf <- readRDS(file.path(data_dir,"sc1conf.rds"))
  meta <- readRDS(file.path(data_dir,"sc1meta.rds"))
  
  conf <- tidyconf(x = conf,num_variable = num_variable)
  meta <- tidymeta(m = meta, c = conf)
  
  metanames <- conf$ID[conf$type == "label" & is.na(conf$fID)]
  for(m in metanames) {
    nl <- nlevels(meta[[m]])
    conf[ID == m]$fID <- paste0(levels(meta[[m]]), collapse = "|")
    conf[ID == m]$fCL <- paste0(
      colorRampPalette(
        RColorBrewer::brewer.pal(12, "Paired")
      )(nl),
      collapse = "|"
    )
    conf[ID == m]$fRow <- ceiling(nl / 4)
    if(nl >= 2) conf[ID == m]$grp <- TRUE
  }
  
  saveRDS(conf, file.path(data_dir,"sc1conf.rds"))
  saveRDS(meta, file.path(data_dir,"sc1meta.rds"))
  
  vlist <- list(
    label = conf[type == "label"]$ID,
    num = conf[type == "num"]$ID,
    dimred = conf[type == "dimred"]$ID
  )
  return(vlist)
}

# Function to extract legend
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

# Plot theme
sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5){
  oupTheme = theme(
    text =             element_text(size = base_size), # , family = "Helvetica"
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line =   element_line(colour = "black"),
    axis.ticks =  element_line(colour = "black", linewidth = base_size / 20),
    axis.title =  element_text(face = "bold"),
    axis.text =   element_text(size = base_size),
    axis.text.x = element_text(angle = Xang, hjust = XjusH),
    legend.position = "bottom",
    legend.key =      element_rect(colour = NA, fill = NA)
  )
  if(!XYval){
    oupTheme = oupTheme + theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  return(oupTheme)
}

### Common plotting functions
# Plot cell information on dimred
scDRcell <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2,
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt, inplab){
  # message("start================================")
  # message("inpsub1: ",inpsub1)
  # message("inpsub2: ",inpsub2)
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]}
  # Prepare ggData
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID,
                       inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID),
                   with = FALSE]
  colnames(ggData) = c("X", "Y", "val", "sub")
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y))
  bgCells = FALSE
  # message("127: ",nrow(ggData))
  # message("inpsub2: ",inpsub2)
  # message("levels(ggData$sub): ",levels(ggData$sub))
  # message(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub))
  # if(is.null(inpsub2)) inpsub2 <- levels(ggData$sub)
  if(length(inpsub2) != nlevels(ggData$sub)){
    bgCells = TRUE
    ggData2 = ggData[!sub %in% inpsub2]
    ggData = ggData[sub %in% inpsub2]
  }
  # message("132: ",nrow(ggData))
  if(inpord == "Max-1st"){
    ggData = ggData[order(val)]
  } else if(inpord == "Min-1st"){
    ggData = ggData[order(-val)]
  } else if(inpord == "Random"){
    ggData = ggData[sample(nrow(ggData))]
  }
  # Do factoring if required
  if(!is.na(inpConf[UI == inp1]$fCL)){
    ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]]
    names(ggCol) = levels(ggData$val)
    ggLvl = levels(ggData$val)[levels(ggData$val) %in% unique(ggData$val)]
    # message(levels(ggData$val))
    # message(levels(ggData$val) %in% unique(ggData$val))
    ggData$val = factor(ggData$val, levels = ggLvl)
    # message(levels(ggData$val))
    ggCol = ggCol[ggLvl]
  }
  
  # Actual ggplot
  ggOut = ggplot(ggData, aes(X, Y, color = val))
  if(bgCells){
    ggOut = ggOut +
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16)
  }
  ggOut = ggOut +
    geom_point(size = inpsiz, shape = 16) + xlab(inpdrX) + ylab(inpdrY) +
    sctheme(base_size = sList[inpfsz], XYval = inptxt)
  if(is.na(inpConf[UI == inp1]$fCL)){
    ggOut = ggOut + scale_color_gradientn("", colours = cList[[inpcol]]) +
      guides(color = guide_colorbar(barwidth = 15))
  } else {
    sListX = min(nchar(paste0(levels(ggData$val), collapse = "")), 200)
    sListX = 0.75 * (sList - (1.5 * floor(sListX/50)))
    ggOut = ggOut + scale_color_manual("", values = ggCol) +
      guides(color = guide_legend(override.aes = list(size = 5),
                                  nrow = inpConf[UI == inp1]$fRow)) +
      theme(legend.text = element_text(size = sListX[inpfsz]))
    if(inplab){
      ggData3 = ggData[, .(X = mean(X), Y = mean(Y)), by = "val"]
      lListX = min(nchar(paste0(ggData3$val, collapse = "")), 200)
      lListX = lList - (0.25 * floor(lListX/50))
      ggOut = ggOut +
        geom_text_repel(data = ggData3, aes(X, Y, label = val),
                        color = "grey10", bg.color = "grey95", bg.r = 0.15,
                        size = lListX[inpfsz], seed = 42)
    }
  }
  if(inpasp == "Square") {
    ggOut = ggOut + coord_fixed(ratio = rat)
  } else if(inpasp == "Fixed") {
    ggOut = ggOut + coord_fixed()
  }
  # message("end================================")
  return(ggOut)
}

# scDRnum <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, 
#                     inpH5, inpGene, inpsplt){ 
#   if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
#   # Prepare ggData 
#   ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), 
#                    with = FALSE] 
#   colnames(ggData) = c("group", "sub") 
#   h5file <- H5File$new(inpH5, mode = "r") 
#   h5data <- h5file[["grp"]][["data"]] 
#   ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
#   ggData[val2 < 0]$val2 = 0 
#   h5file$close_all() 
#   if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
#     ggData = ggData[sub %in% inpsub2] 
#   } 
#   
#   # Split inp1 if necessary 
#   if(is.na(inpConf[UI == inp1]$fCL)){ 
#     if(inpsplt == "Quartile"){nBk = 4} 
#     if(inpsplt == "Decile"){nBk = 10} 
#     ggData$group = cut(ggData$group, breaks = nBk) 
#   } 
#   
#   # Actual data.table 
#   ggData$express = FALSE 
#   ggData[val2 > 0]$express = TRUE 
#   ggData1 = ggData[express == TRUE, .(nExpress = .N), by = "group"] 
#   ggData = ggData[, .(nCells = .N), by = "group"] 
#   ggData = ggData1[ggData, on = "group"] 
#   ggData = ggData[, c("group", "nCells", "nExpress"), with = FALSE] 
#   ggData[is.na(nExpress)]$nExpress = 0 
#   ggData$pctExpress = 100 * ggData$nExpress / ggData$nCells 
#   ggData = ggData[order(group)] 
#   colnames(ggData)[3] = paste0(colnames(ggData)[3], "_", inp2) 
#   return(ggData) 
# } 
# Plot gene expression on dimred
scDRgene <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2,
                     inpH5, inpGene,
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt){
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]}
  # Prepare ggData
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID,
                       inpConf[UI == inpsub1]$ID),
                   with = FALSE]
  colnames(ggData) = c("X", "Y", "sub")
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y))

  # message(inp1)
  # message("inp1 class is ", class(inp1))
  # message("inp1 length is ", length(inp1))
  # if(is.null(inp1)) {
  #   # message("inp1 is null")
  #   ggData$val = 0
  # } else {
  #   # message("inp1 is not null")
  #   
  # }
  if(is.null(inp1)) {
    ggData$val = NA  # 设置无效值避免影响颜色映射
    fixed_color <- "snow2"  # 固定颜色标志
  } else {
    h5file <- H5File$new(inpH5, mode = "r")
    h5data <- h5file[["grp"]][["data"]]
    ggData$val = h5data$read(args = list(inpGene[inp1], quote(expr=)))
    ggData[val < 0]$val = 0
    h5file$close_all()
    fixed_color <- NULL  # 保持正常颜色映射
  }
  bgCells = FALSE
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){
    bgCells = TRUE
    ggData2 = ggData[!sub %in% inpsub2]
    ggData = ggData[sub %in% inpsub2]
  }
  if(inpord == "Max-1st"){
    ggData = ggData[order(val)]
  } else if(inpord == "Min-1st"){
    ggData = ggData[order(-val)]
  } else if(inpord == "Random"){
    ggData = ggData[sample(nrow(ggData))]
  }

  # Actual ggplot
  ggOut = ggplot(ggData, aes(X, Y))
  if(bgCells){
    ggOut = ggOut +
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16)
  }
  if(is.null(inp1)) {
    # 当inp1为NULL时，所有点强制为snow2
    ggOut = ggOut +
      geom_point(color = "snow2", size = inpsiz, shape = 16)
  } else {
    # 正常颜色映射逻辑
    ggOut = ggOut +
      geom_point(aes(color = val), size = inpsiz, shape = 16) +
      scale_color_gradientn(inp1, colours = cList[[inpcol]]) +
      guides(color = guide_colorbar(barwidth = 15))
  }
  # 添加公共图形元素
  ggOut = ggOut +
    xlab(inpdrX) + ylab(inpdrY) +
    sctheme(base_size = sList[inpfsz], XYval = inptxt)
  
  if(inpasp == "Square") {
    ggOut = ggOut + coord_fixed(ratio = rat)
  } else if(inpasp == "Fixed") {
    ggOut = ggOut + coord_fixed()
  }
  return(ggOut)
}

# Plot gene coexpression on dimred
bilinear <- function(x,y,xy,Q11,Q21,Q12,Q22){
  oup = (xy-x)*(xy-y)*Q11 + x*(xy-y)*Q21 + (xy-x)*y*Q12 + x*y*Q22
  oup = oup / (xy*xy)
  return(oup)
}

scDRcoex <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inp2,
                     inpsub1, inpsub2, inpH5, inpGene,
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt){
  
  # 处理输入参数
  has1 <- !is.null(inp1) && inp1 %in% names(inpGene)
  has2 <- !is.null(inp2) && inp2 %in% names(inpGene)
  
  # 准备基础数据
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, 
                       inpConf[UI == inpdrY]$ID,
                       inpConf[UI == inpsub1]$ID),
                   with = FALSE]
  colnames(ggData) = c("X", "Y", "sub")
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y))
  
  # 读取基因表达数据
  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  
  if(has1) {
    ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=)))
    ggData[val1 < 0]$val1 = 0
  }
  
  if(has2) {
    ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=)))
    ggData[val2 < 0]$val2 = 0
  }
  h5file$close_all()
  
  # 处理分组筛选
  bgCells = FALSE
  if(length(inpsub2) != 0 && length(inpsub2) != nlevels(ggData$sub)){
    bgCells = TRUE
    ggData2 = ggData[!sub %in% inpsub2]
    ggData = ggData[sub %in% inpsub2]
  }
  
  # 颜色处理逻辑
  if(has1 && has2) {
    # 双基因共表达颜色混合
    cInp = strsplit(inpcol, "; ")[[1]]
    if(cInp[1] == "Red (Gene1)"){
      c10 = c(255,0,0)
    } else if(cInp[1] == "Orange (Gene1)"){
      c10 = c(255,140,0)
    } else {
      c10 = c(0,255,0)
    }
    if(cInp[2] == "Green (Gene2)"){
      c01 = c(0,255,0)
    } else {
      c01 = c(0,0,255)
    }
    c00 = c(217,217,217) ; c11 = c10 + c01
    nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2
    gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1)))
    gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid
    gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid
    gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1])
    gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2])
    gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3])
    gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255)
    gg = gg[, c("v1", "v2", "cMix")]
    
    # Map colours
    ggData$v1 = round(nTot * ggData$val1 / max(ggData$val1))
    ggData$v2 = round(nTot * ggData$val2 / max(ggData$val2))
    ggData$v0 = ggData$v1 + ggData$v2
    ggData = gg[ggData, on = c("v1", "v2")]
    if(inpord == "Max-1st"){
      ggData = ggData[order(v0)]
    } else if(inpord == "Min-1st"){
      ggData = ggData[order(-v0)]
    } else if(inpord == "Random"){
      ggData = ggData[sample(nrow(ggData))]
    }
  } else if(has1 || has2) {
    # 单基因颜色梯度
    geneCol <- if(has1) c(255,0,0) else c(0,0,255) # 默认红/蓝
    ggData$cMix <- rgb(
      geneCol[1] * ggData[[if(has1) "val1" else "val2"]] / 
        max(ggData[[if(has1) "val1" else "val2"]]),
      geneCol[2],
      geneCol[3],
      maxColorValue = 255
    )
  } else {
    # 无基因时全灰
    ggData$cMix <- "snow2"
  }
  
  # 绘图逻辑
  ggOut <- ggplot(ggData, aes(X, Y)) # + 
    # theme_minimal()
  
  if(bgCells) {
    ggOut <- ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16)
  }
  
  if(has1 || has2) {
    # 有基因时的主图层
    ggOut <- ggOut +
      geom_point(size = inpsiz, shape = 16, aes(color = cMix)) +
      scale_color_identity()
  } else {
    # 无基因时全灰
    ggOut <- ggOut +
      geom_point(size = inpsiz, shape = 16, color = "snow2")
  }
  
  # 坐标轴和主题设置
  ggOut <- ggOut +
    xlab(inpdrX) + ylab(inpdrY) +
    sctheme(base_size = sList[inpfsz], XYval = inptxt)
    # theme(axis.text = element_text(size = inpfsz),
    #       axis.title = element_text(size = inpfsz + 2))
  
  # 长宽比处理
  if(inpasp == "Square") {
    ggOut <- ggOut + coord_fixed(ratio = rat)
  } else if(inpasp == "Fixed") {
    ggOut <- ggOut + coord_fixed()
  }
  
  return(ggOut)
}

scDRcoexLeg <- function(inp1, inp2, inpcol, inpfsz){
  # Generate coex color palette
  cInp = strsplit(inpcol, "; ")[[1]]
  if(cInp[1] == "Red (Gene1)"){
    c10 = c(255,0,0)
  } else if(cInp[1] == "Orange (Gene1)"){
    c10 = c(255,140,0)
  } else {
    c10 = c(0,255,0)
  }
  if(cInp[2] == "Green (Gene2)"){
    c01 = c(0,255,0)
  } else {
    c01 = c(0,0,255)
  }
  c00 = c(217,217,217) ; c11 = c10 + c01
  nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2
  gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1)))
  gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid
  gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid
  gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1])
  gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2])
  gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3])
  gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255)
  gg = gg[, c("v1", "v2", "cMix")]

  # Actual ggplot
  ggOut = ggplot(gg, aes(v1, v2)) +
    geom_tile(fill = gg$cMix) +
    xlab(inp1) + ylab(inp2) + coord_fixed(ratio = 1) +
    scale_x_continuous(breaks = c(0, nTot), label = c("low", "high")) +
    scale_y_continuous(breaks = c(0, nTot), label = c("low", "high")) +
    sctheme(base_size = sList[inpfsz], XYval = TRUE)
  return(ggOut)
}

scDRcoexNum <- function(inpConf, inpMeta, inp1, inp2,
                        inpsub1, inpsub2, inpH5, inpGene){
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]}
  
  has1 <- !is.null(inp1) && inp1 %in% names(inpGene)
  has2 <- !is.null(inp2) && inp2 %in% names(inpGene)
  
  # Prepare ggData
  ggData = inpMeta[, c(inpConf[UI == inpsub1]$ID), with = FALSE]
  colnames(ggData) = c("sub")
  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  if(has1) {
    ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=)))
    ggData[val1 < 0]$val1 = 0
  }
  if(has2) {
    ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=)))
    ggData[val2 < 0]$val2 = 0
  }
  h5file$close_all()
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){
    ggData = ggData[sub %in% inpsub2]
  }

  # Actual data.table
  ggData$express = "none"
  if(has1 && !has2) {
    ggData[val1 > 0]$express = inp1
    ggData$express = factor(ggData$express, levels = unique(c("both", inp1, "none")))
  }
  if(has2 && !has1) {
    ggData[val2 > 0]$express = inp2
    ggData$express = factor(ggData$express, levels = unique(c("both", inp2, "none")))
  }
  if(has1 && has2) {
    ggData[val1 > 0]$express = inp1
    ggData[val2 > 0]$express = inp2
    ggData[val1 > 0 & val2 > 0]$express = "both"
    ggData$express = factor(ggData$express, levels = unique(c("both", inp1, inp2, "none")))
  }
  ggData = ggData[, .(nCells = .N), by = "express"]
  ggData$percent = 100 * ggData$nCells / sum(ggData$nCells)
  ggData = ggData[order(express)]
  colnames(ggData)[1] = "expression > 0"
  return(ggData)
}

# Plot violin / boxplot
scVioBox <- function(inpConf, inpMeta, inp1, inp2,
                     inpsub1, inpsub2, inpH5, inpGene,
                     inptyp, inppts, inpsiz, inpfsz){
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]}
  # Prepare ggData
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID),
                   with = FALSE]
  colnames(ggData) = c("X", "sub")

  # Load in either cell meta or gene expr
  if(is.null(inp2)) inp2 <- ""
  if(inp2 %in% inpConf$UI){
    ggData$val = inpMeta[[inpConf[UI == inp2]$ID]]
  } else {
    h5file <- H5File$new(inpH5, mode = "r")
    h5data <- h5file[["grp"]][["data"]]
    if(inp2 == "") {
      ggData$val = 0
    } else {
      ggData$val = h5data$read(args = list(inpGene[inp2], quote(expr=)))
      ggData[val < 0]$val = 0
    }
    set.seed(42)
    tmpNoise = rnorm(length(ggData$val)) * diff(range(ggData$val)) / 1000
    ggData$val = ggData$val + tmpNoise
    h5file$close_all()
  }
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){
    ggData = ggData[sub %in% inpsub2]
  }

  # Do factoring
  ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]]
  names(ggCol) = levels(ggData$X)
  ggLvl = levels(ggData$X)[levels(ggData$X) %in% unique(ggData$X)]
  ggData$X = factor(ggData$X, levels = ggLvl)
  ggCol = ggCol[ggLvl]

  # Actual ggplot
  if(inptyp == "violin"){
    ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_violin(scale = "width")
  } else {
    ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_boxplot()
  }
  if(inppts){
    ggOut = ggOut + geom_jitter(size = inpsiz, shape = 16)
  }
  ggOut = ggOut + xlab(inp1) + ylab(inp2) +
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
    scale_fill_manual("", values = ggCol) +
    theme(legend.position = "none")
  return(ggOut)
}

# # Plot proportion plot 
scProp <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2,
                   inptyp, inpflp, inpfsz){
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]}
  # Prepare ggData
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inp2]$ID,
                       inpConf[UI == inpsub1]$ID),
                   with = FALSE]
  colnames(ggData) = c("X", "grp", "sub")
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){
    ggData = ggData[sub %in% inpsub2]
  }
  ggData = ggData[, .(nCells = .N), by = c("X", "grp")]
  ggData = ggData[, {tot = sum(nCells)
  .SD[,.(pctCells = 100 * sum(nCells) / tot,
         nCells = nCells), by = "grp"]}, by = "X"]

  # Do factoring
  ggCol = strsplit(inpConf[UI == inp2]$fCL, "\\|")[[1]]
  names(ggCol) = levels(ggData$grp)
  ggLvl = levels(ggData$grp)[levels(ggData$grp) %in% unique(ggData$grp)]
  ggData$grp = factor(ggData$grp, levels = ggLvl)
  ggCol = ggCol[ggLvl]

  # Actual ggplot
  if(inptyp == "Proportion"){
    ggOut = ggplot(ggData, aes(X, pctCells, fill = grp)) +
      geom_col() + ylab("Cell Proportion (%)")
  } else {
    ggOut = ggplot(ggData, aes(X, nCells, fill = grp)) +
      geom_col() + ylab("Number of Cells")
  }
  if(inpflp){
    ggOut = ggOut + coord_flip()
  }
  ggOut = ggOut + xlab(inp1) +
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
    scale_fill_manual("", values = ggCol) +
    theme(legend.position = "right")
  return(ggOut)
}

# Get gene list
scGeneList <- function(inp, inpGene) {
  geneList = data.table(
    gene = unique(trimws(strsplit(inp, ",|;|\n")[[1]])),
    present = TRUE
  )
  geneList[!gene %in% names(inpGene)]$present = FALSE
  return(geneList)
}

# Plot gene expression bubbleplot / heatmap
scBubbHeat <- function(inpConf, inpMeta, inp, inpGrp, inpPlt,
                       inpsub1, inpsub2, inpH5, inpGene, inpScl, inpRow, inpCol,
                       inpcols, inpfsz, save = FALSE){
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]}
  # Identify genes that are in our dataset
  geneList = scGeneList(inp, inpGene)
  geneList = geneList[present == TRUE]
  shiny::validate(need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!"))
  shiny::validate(need(nrow(geneList) > 1, "Please input at least 2 genes to plot!"))

  # Prepare ggData
  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData = data.table()
  for(iGene in geneList$gene){
    tmp = inpMeta[, c("sampleID", inpConf[UI == inpsub1]$ID), with = FALSE]
    colnames(tmp) = c("sampleID", "sub")
    tmp$grpBy = inpMeta[[inpConf[UI == inpGrp]$ID]]
    tmp$geneName = iGene
    tmp$val = h5data$read(args = list(inpGene[iGene], quote(expr=)))
    ggData = rbindlist(list(ggData, tmp))
  }
  h5file$close_all()
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){
    ggData = ggData[sub %in% inpsub2]
  }
  shiny::validate(need(uniqueN(ggData$grpBy) > 1, "Only 1 group present, unable to plot!"))

  # Aggregate
  ggData$val = expm1(ggData$val)
  ggData = ggData[, .(val = mean(val), prop = sum(val>0) / length(sampleID)),
                  by = c("geneName", "grpBy")]
  ggData$val = log1p(ggData$val)

  # Scale if required
  colRange = range(ggData$val)
  if(inpScl){
    ggData[, val:= scale(val), keyby = "geneName"]
    colRange = c(-max(abs(range(ggData$val))), max(abs(range(ggData$val))))
  }

  # hclust row/col if necessary
  ggMat = dcast.data.table(ggData, geneName~grpBy, value.var = "val")
  tmp = ggMat$geneName
  ggMat = as.matrix(ggMat[, -1])
  rownames(ggMat) = tmp
  if(inpRow){
    hcRow = ggdendro::dendro_data(as.dendrogram(hclust(dist(ggMat))))
    ggRow = ggplot() + coord_flip() +
      geom_segment(data = hcRow$segments, aes(x=x,y=y,xend=xend,yend=yend)) +
      scale_y_continuous(breaks = rep(0, uniqueN(ggData$grpBy)),
                         labels = unique(ggData$grpBy), expand = c(0, 0)) +
      scale_x_continuous(breaks = seq_along(hcRow$labels$label),
                         labels = hcRow$labels$label, expand = c(0, 0.5)) +
      sctheme(base_size = sList[inpfsz]) +
      theme(axis.title = element_blank(), axis.line = element_blank(),
            axis.ticks = element_blank(), axis.text.y = element_blank(),
            axis.text.x = element_text(color="white", angle = 45, hjust = 1))
    ggData$geneName = factor(ggData$geneName, levels = hcRow$labels$label)
  } else {
    ggData$geneName = factor(ggData$geneName, levels = rev(geneList$gene))
  }
  if(inpCol){
    hcCol = ggdendro::dendro_data(as.dendrogram(hclust(dist(t(ggMat)))))
    ggCol = ggplot() +
      geom_segment(data = hcCol$segments, aes(x=x,y=y,xend=xend,yend=yend)) +
      scale_x_continuous(breaks = seq_along(hcCol$labels$label),
                         labels = hcCol$labels$label, expand = c(0.05, 0)) +
      scale_y_continuous(breaks = rep(0, uniqueN(ggData$geneName)),
                         labels = unique(ggData$geneName), expand=c(0,0)) +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      theme(axis.title = element_blank(), axis.line = element_blank(),
            axis.ticks = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_text(color = "white"))
    ggData$grpBy = factor(ggData$grpBy, levels = hcCol$labels$label)
  }

  # Actual plot according to plottype
  if(inpPlt == "Bubbleplot"){
    # Bubbleplot
    ggOut = ggplot(ggData, aes(grpBy, geneName, color = val, size = prop)) +
      geom_point() +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      scale_x_discrete(expand = c(0.05, 0)) +
      scale_y_discrete(expand = c(0, 0.5)) +
      scale_size_continuous("proportion", range = c(0, 8),
                            limits = c(0, 1), breaks = c(0.00,0.25,0.50,0.75,1.00)) +
      scale_color_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) +
      guides(color = guide_colorbar(barwidth = 15)) +
      theme(axis.title = element_blank(), legend.box = "vertical")
  } else {
    # Heatmap
    ggOut = ggplot(ggData, aes(grpBy, geneName, fill = val)) +
      geom_tile() +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      scale_x_discrete(expand = c(0.05, 0)) +
      scale_y_discrete(expand = c(0, 0.5)) +
      scale_fill_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) +
      guides(fill = guide_colorbar(barwidth = 15)) +
      theme(axis.title = element_blank())
  }

  # Final tidy
  ggLeg = g_legend(ggOut)
  ggOut = ggOut + theme(legend.position = "none")
  if(!save){
    if(inpRow & inpCol){ggOut =
      gridExtra::grid.arrange(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),
                   layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))
    } else if(inpRow){ggOut =
      gridExtra::grid.arrange(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),
                   layout_matrix = rbind(c(1,3),c(2,NA)))
    } else if(inpCol){ggOut =
      gridExtra::grid.arrange(ggOut, ggLeg, ggCol, heights = c(1,7,2),
                   layout_matrix = rbind(c(3),c(1),c(2)))
    } else {ggOut =
      gridExtra::grid.arrange(ggOut, ggLeg, heights = c(7,2),
                   layout_matrix = rbind(c(1),c(2)))
    }
  } else {
    if(inpRow & inpCol){ggOut =
      arrangeGrob(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),
                  layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))
    } else if(inpRow){ggOut =
      arrangeGrob(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),
                  layout_matrix = rbind(c(1,3),c(2,NA)))
    } else if(inpCol){ggOut =
      arrangeGrob(ggOut, ggLeg, ggCol, heights = c(1,7,2),
                  layout_matrix = rbind(c(3),c(1),c(2)))
    } else {ggOut =
      arrangeGrob(ggOut, ggLeg, heights = c(7,2),
                  layout_matrix = rbind(c(1),c(2)))
    }
  }
  return(ggOut)
}

plot_volcano <- function(data, FC = 1, PValue = 0.05, label_num = 5, volcano_title = "Volcano") {
  # data$gene <- rownames(data)
  if (!all(c("p_val_adj", "avg_log2FC", "gene") %in% names(data)))
    stop("colnames must contain p_val_adj, avg_log2FC, gene")
  
  data <- data |> filter(p_val_adj < PValue)
  # 判断每个基因的上下调
  data$sig[-1*log10(data$p_val_adj) >= -1*log10(PValue) & data$avg_log2FC >= FC] <- "Up"
  data$sig[-1*log10(data$p_val_adj) >= -1*log10(PValue) & data$avg_log2FC <= -FC] <- "Down"
  data$sig[is.na(data$sig)] <- "NotSig"
  
  # label df
  data_filter <- data[which(data$p_val_adj < PValue & abs(data$avg_log2FC) > FC),]
  data_up <- data_filter[which(data_filter$sig == "Up"),]
  data_up_p_min <- data_up[order(data_up$p_val_adj),] |> head(label_num)
  data_up_log2fc_max <- data_up[order(data_up$avg_log2FC,decreasing = TRUE),] |> head(label_num)
  data_down <- data_filter[which(data_filter$sig == "Down"),]
  data_down_p_min <- data_down[order(data_down$p_val_adj),] |> head(label_num)
  data_down_log2fc_max <- data_down[order(data_down$avg_log2FC),] |> head(label_num)
  data_label <- rbind.data.frame(data_up_p_min,data_down_p_min,data_up_log2fc_max,data_down_log2fc_max) |>
    unique()
  
  # 绘制火山图
  volcano_plot <- ggplot(data, aes(avg_log2FC, -1*log10(p_val_adj))) +
    geom_point(aes(color = sig)) +
    labs(title="Volcano plot", x="log2 (FC)", y="-log10 (PValue)") +
    geom_hline(yintercept=-log10(PValue), linetype=2) +
    geom_vline(xintercept=c(-FC, FC), linetype=2) +
    scale_color_manual(values=c("Up" = "#ff4757", "Down" = "#546de5", "NotSig" = "#d2dae2")) +
    # 基因标签
    geom_label_repel(
      data = data_label,
      aes(label = gene)#,
      # size = label_font_si
    ) +
    ggtitle(volcano_title) +
    theme_bw() +
    theme(legend.position = "right")
  return(volcano_plot)
}

plot_fc_pctdiff <- function(data, FC = 1, PValue = 0.05, top_n_diff = 5, diff_title = NULL) {
  # data$gene <- rownames(data)
  if (!all(c("p_val_adj", "avg_log2FC", "gene") %in% names(data)))
    stop("colnames must contain p_val_adj, avg_log2FC, gene")
  # data <- data|>mutate(diff = pct.1-pct)
  
  data <- data |> filter(p_val_adj < PValue)
  # 判断每个基因的上下调
  # data$sig[ (-1*log10(data$p_val_adj) < -1*log10(PValue)|data$p_val_adj=="NA")|(data$avg_log2FC < FC)& data$avg_log2FC > -FC] <- "NotSig"
  data$sig[-1*log10(data$p_val_adj) >= -1*log10(PValue) & data$avg_log2FC >= FC] <- "Up"
  data$sig[-1*log10(data$p_val_adj) >= -1*log10(PValue) & data$avg_log2FC <= -FC] <- "Down"
  data$sig[is.na(data$sig)] <- "NotSig"
  
  # 绘图
  p <- ggplot(data,aes(x = pct.1 - pct.2, y = avg_log2FC)) +
    geom_point(aes(color = sig)) +
    labs(title = diff_title, x="Percentage Difference", y="Log2 Fold Change") +
    scale_color_manual(values=c("Up" = "#ff4757", "Down" = "#546de5", "NotSig" = "#d2dae2"))
  # 基因标签
  p <- p + ggrepel::geom_label_repel(
    data = data |> filter(sig %in% c("Up","Down")) |> filter(abs(avg_log2FC) > FC) |>
      mutate(diff = pct.1 - pct.2) |> group_by(sig) |>
      top_n(top_n_diff,abs(diff)),
    aes(label = gene)#,label.size = 3
  )
  p <- p + theme_bw() + geom_hline(yintercept = 0,linetype = "dashed") +
    geom_vline(xintercept = 0,linetype = "dashed")
  return(p)
}

plot_deg <- function(data, FC = 1, Pvalue = 0.05, label_num = 5, plot_type) {
  if(plot_type == "Volcano") {
    p <- plot_volcano(
      data = data,
      FC = FC,
      PValue = Pvalue,
      label_num = label_num,
      volcano_title = NULL
    )
  } else {
    p <- plot_fc_pctdiff(
      data = data,
      FC = FC,
      PValue = Pvalue,
      top_n_diff = label_num,
      diff_title = NULL
    )
  }
  p <- p + sctheme(base_size = 18) # base_size = sList[inpfsz]
  return(p)
}

filterDegTable <- function(df, p = 0.05, fc = 1, absolute = TRUE) {
  p_row <- df$p_val_adj < p
  fc_row <- ifelse(absolute,
                   df$avg_log2FC > fc | df$avg_log2FC < -fc,
                   df$avg_log2FC > fc)
  if (absolute) {
    fc_row <- df$avg_log2FC > fc | df$avg_log2FC < -fc
  } else{
    fc_row <- df$avg_log2FC > fc
  }
  df <- df[p_row & fc_row, ]
  return(df)
}