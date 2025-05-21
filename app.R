suppressWarnings(library(shiny))
suppressWarnings(library(shinyjs)) # 隐藏UI
suppressWarnings(library(shinyhelper))
library(bslib) # theme
# library(phsstyles) # for phs colour palette

suppressWarnings(library(data.table))
suppressWarnings(library(dplyr))

suppressWarnings(library(ggplot2))
suppressWarnings(library(ggrepel))
suppressWarnings(library(hdf5r))

# library(Matrix)
suppressWarnings(library(DT))
# library(magrittr)
# library(ggdendro)
# library(gridExtra)

source("global.R")
source("utils.R")
source("main_ui.R")
source("main_server.R")
source("modules.R")
shinyApp(ui = ui,server = server,enableBookmarking = "url")


# test-----
# conf <- readRDS("case_shiny_cell_bookmark/case/sc1conf.rds")
# conf <- tidyconf(conf)
# meta <- readRDS("case_shiny_cell_bookmark/case/sc1meta.rds")
# meta <- tidymeta(m = meta,c = conf)
# gene <- readRDS("case_shiny_cell_bookmark/case/sc1gene.rds")
# sc1def  <- readRDS("case_shiny_cell_bookmark/case/sc1def.rds")
# sc1deg  <- readRDS(file.path(input$dataname,"sc1deg.rds"))|>reactive()
# 调整conf文件-----
# glimpse_meta(data_dir = "case_shiny_cell_bookmark/pbmc3k/")
# adjust_conf_files(data_dir = "case_shiny_cell_bookmark/pbmc3k")
