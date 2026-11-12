library(bs4Dash)
library(shinyWidgets)
library(shinyscreenshot)
library(shinyalert)
# library(waiter)
library(ggdendro)
library(ggplot2)
library(plotly)
library(DT)
library(dplyr)
# remotes::install_github("deepanshu88/summaryBox")
# library(summaryBox)
# devtools::install_github("glin/reactable")
library(reactable)
# library(circlize)
## library(colorRamp2)
# library(shinyjs)
library(conflicted)
library(ComplexHeatmap)
library(htmltools)
library(ggprism)

# database
# library(AnnotationDbi)

# WGCNA 
# library(WGCNA)

WGCNA::enableWGCNAThreads(nThreads = 20)

# Draw sft plot
# library(ggrepel)
# library(patchwork)

# annotation
# library(ClusterGVis)
# library(clusterProfiler)

# network
# library(purrr)
# library(tidygraph)
library(ggraph)
# library(igraph)

#####################################
# conflicted
conflicts_prefer(bs4Dash::box)
conflicts_prefer(bs4Dash::tabsetPanel)
conflicts_prefer(bs4Dash::column)
conflicts_prefer(bs4Dash::actionButton)
conflicts_prefer(plotly::layout)
conflicts_prefer(shinyscreenshot::runExample)
conflicts_prefer(WGCNA::cor)

# set the max upload file size
options(shiny.maxRequestSize = 100 * 1024^2) 

options(
  shiny.useragg = FALSE,
  shiny.usecairo = FALSE,
  bitmapType = "cairo"
)

# Define UI for application
ui <- dashboardPage(
  # preloader = list(html = tagList(spin_1(), "Loading ..."), color = "#adb5bd"),
  fullscreen = TRUE,
  dark = NULL,
  help = NULL,
  header = dashboardHeader(
    title = "",
    center = tags$h2(style="color:white;font-family:Comic Sans MS;","iWGCNA：An integrative Web Server for Weighted Gene Co-expression Network Analysis"), 
    status = "info"
  ),
  sidebar = bs4DashSidebar(
    skin = "light",
    status = "info",
    elevation = 3,
    collapsed = TRUE,
    sidebarMenu(
      id = "steps",
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("1. Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("2. Preprocessing", tabName = "preprocessing", icon = icon("filter")),
      menuItem("3. Network Construction and Gene Group (Module) Detection", tabName = "create", icon = icon("check-to-slot")),
      menuItem("4. Module-Trait Association", tabName = "interested", icon = icon("location-pin")),
      menuItem("5. Hub Gene Identification", tabName = "hub", icon = icon("star", class = "fa-solid")),
      menuItem("6. Functional Enrichment Analysis", tabName = "function", icon = icon("map")),
      menuItem("Tutorial", tabName = "tutorial", icon = icon("book")),
      menuItem("News", tabName = "news", icon = icon("newspaper"))
      # menuItem("GEO Data Retrieval", tabName = "geo", icon = icon("database"))
    )
  ),
  body = dashboardBody(
    
    shinyjs::useShinyjs(),
    
    #here's where you throw the css into the header
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    ),
    tags$head(
      tags$style(HTML("
    /* 滑鼠懸停時的顏色 */
    .nav-item > .nav-link:hover {
      background-color: #17a2b8 !important;
      color: #fff !important;
    }

    /* 被選中的頁面title顏色 */
    .nav-item > .nav-link.active,
    .nav-item > .nav-link.active:hover,
    .nav-item > .nav-link.active:focus {
      background-color: #17a2b8 !important;
      color: #fff !important;
      border: none !important;
    }
  "))
    ),
    tags$style("
      .bold {
        font-weight: bold;
      }
      .highlight {
        font-weight: bold;
        border: 2px solid #17a2b8;    /* 外框 */
      }
    "),
    waiter::autoWaiter(html = waiter::spin_3circles()),
    tabItems(
      tabItem(tabName = "home", source("UI/0_home.R", local = TRUE)$value),
      tabItem(tabName = "upload", source("UI/1_upload.R", local = TRUE)$value),
      tabItem(tabName = "preprocessing", source("UI/2_preprocessing.R", local = TRUE)$value),
      tabItem(tabName = "create", source("UI/3_create.R", local = TRUE)$value),
      tabItem(tabName = "interested", source("UI/4_interested.R", local = TRUE)$value),
      tabItem(tabName = "hub", source("UI/5_hub.R", local = TRUE)$value),
      tabItem(tabName = "function", source("UI/6_function.R", local = TRUE)$value),
      tabItem(tabName = "tutorial", source("UI/0_tutorial.R", local = TRUE)$value),
      tabItem(tabName = "news", source("UI/0_news.R", local = TRUE)$value)
      # tabItem(tabName = "geo", source("UI/#_geo.R", local = TRUE)$value)
    )
  )
)

server <- function(input, output, session) {
  
  source("Server/0_home.R", local = TRUE)$value
  source("Server/1_upload.R", local = TRUE)$value
  source("Server/2_preprocessing.R", local = TRUE)$value
  source("Server/3_create.R", local = TRUE)$value
  source("Server/4_interested.R", local = TRUE)$value
  source("Server/5_hub.R", local = TRUE)$value
  source("Server/6_function.R", local = TRUE)$value
  source("Server/0_tutorial.R", local = TRUE)$value
  source("Server/0_news.R", local = TRUE)$value
  # source("Server/#_geo.R", local = TRUE)$value
}


shinyApp(ui = ui, server = server)
