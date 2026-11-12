# UI
tabItem(
  tabName = 'create',
  fluidRow(
    box(
      id = 'create_network',
      height = "auto",
      width = 12,
      elevation = 4,
      title = div(
        style = "font-size: 20px; font-weight: bold;", 
        "3. Network Construction and Gene Group (Module) Detection"
      ),
      collapsible = TRUE,
      solidHeader = FALSE,
      maximizable = TRUE,
      status = "info",
      tabsetPanel(
        tabPanel(
          title = "3a. Soft-Threshold Power Evaluation (β)",
          value = "power_tab",
          icon = icon("table"),
          fluidRow(
            column(
              width = 3,
              div(style = "margin-bottom: 20px;"),
              # singed and unsinged
              prettyRadioButtons(
                inputId = "Networktype",
                label = tags$span(
                  style = "font-size: 18px;", "Network Type: ", " ", 
                  actionButton("network",label = "",icon = icon("info"), class = "btn-info btn-xs")
                ),
                choices = c("Unsigned" = "unsigned", "Signed" = "signed","Signed Hybrid" = "signed hybrid"),
                selected = "signed",
                inline = TRUE, status = "info", fill = TRUE
              ),
              # choose R2 cutoff
              sliderInput(
                inputId = 'cutoff',
                label = tags$span(style = "font-size: 18px;", HTML('R<sup>2</sup> Cutoff: '),
                                  " ",
                                  actionButton("Cutoff",label = "",icon = icon("info"), class = "btn-info btn-xs")), 
                min = 0.0, max = 1.0, value = 0.9, step = 0.05
              ),
              uiOutput("recommend"),
              hr(),
              bs4Dash::actionButton(
                "cre_run",em("Evaluate Soft-Threshold Power"), status = "primary",
                style = "box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.4); border-radius: 8px;"
              ),
              div(style = "margin-bottom: 20px;"),
              shinyjs::hidden(div(
                id = "step3b_ui", 
                bs4Dash::actionButton("step3b",em("STEP 3b. Network Parameter Settings"),
                                      status = "info", icon = icon("forward"))
                )
              )
            ),
            column(
              width = 9,
              div(style = "margin-bottom: 20px;"),
              fluidRow(
                uiOutput("cre_choose_ui")
              ),
              conditionalPanel(
                condition = "input.cre_choose === 'plot'  && input.cre_run >0",
                dropdownButton(
                  sliderInput("sftplot_width", "Plot Width:", min = 500, max = 1200, value = 700, step = 10),
                  sliderInput("sftplot_height", "Plot Height:", min = 300, max = 1200, value = 400, step = 10),
                  circle = F, status = "danger", icon = icon("cog"), width = "300px", size = 'sm',
                  tooltip = tooltipOptions(title = "Adjust Plot Size")
                ),
                tags$div(
                  style = "overflow-x: auto; overflow-y: auto; width: 100%;", 
                  plotOutput("cre_sftplot", height = "auto", width = "auto")
                ),
                shinyjs::hidden(div(id = "cre_download", downloadButton("credownload", "Download Plot")))
              ),
              conditionalPanel(
                condition = "input.cre_choose === 'table'",
                reactableOutput("cre_sfttab"),
                shinyjs::hidden(div(id = "cre_download_2", downloadButton("credownload2", "Download Data")))
              )
            )
          )
        ),
        tabPanel(
          title = "3b. Network Parameter Settings",
          icon = icon("tree"),
          value = "module_detection_tab",
          fluidRow(
            column(
              width = 3,
              div(style = "margin-bottom: 20px;"),
              
              # select power 
              uiOutput("powerselect_ui"),
              
              # select minModuleSize
              sliderInput(
                inputId = 'minMsize',
                label = tags$span(style = "font-size: 18px;", "Minimum Gene Group Size: ",
                                  " ",
                                  actionButton("mms2",label = "",icon = icon("info"), class = "btn-info btn-xs")), 
                min = 0, max = 200, value = 30,
              ),
              
              # select mergeCutHeight  
              sliderInput(
                inputId = 'mch',
                label = tags$span(style = "font-size: 18px;", "Gene Group Merge Cut Height: ",
                                  " ",
                                  actionButton("mch2",label = "",icon = icon("info"), class = "btn-info btn-xs")), 
                min = 0.0, max = 1.0, value = 0.25,
              ),
              
              # prettyRadioButtons(
              #   inputId = "TOMtype",
              #   label = tags$span(
              #     style = "font-size: 18px;", "Topological Overlap Matrix (TOM) Type: ", " ", 
              #     actionButton("tom",label = "",icon = icon("info"), class = "btn-info btn-xs")
              #   ),
              #   choices = c("None" = "none", "Unsigned" = "unsigned", "Signed" = "signed"),
              #   selected = "signed",
              #   inline = TRUE, status = "info", fill = TRUE
              # ),
              hr(),
              bs4Dash::actionButton(
                "cre_run2",em("Construct Network"), status = "primary",
                style = "box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.4); border-radius: 8px;"
              ),
              div(style = "margin-bottom: 20px;"),
              fluidRow(
                column(
                  width = 12,
                  shinyjs::hidden(div(
                    id = "step4_ui", 
                    bs4Dash::actionButton("step4",em("STEP 4. Gene Group-Trait Association"),
                                          status = "info", icon = icon("forward"))
                  
                    )
                  ),
                  shinyjs::hidden(div(
                    id = "step6_ui", 
                    div(style = "margin-bottom: 10px;"),
                    bs4Dash::actionButton("step6",em("STEP 6 Functional Enrichment Analysis"),
                                          status = "info", icon = icon("forward"))
                    )
                  )
                )
              )
            ),
            column(
              width = 9,
              div(style = "margin-bottom: 20px;"),
              fluidRow(
                uiOutput("cre_choose2_ui")
              ),
              conditionalPanel(
                condition = "input.cre_choose2 === 'plot'  && input.cre_run2 >0",
                prettyRadioButtons(
                  inputId = "cre_plot_type", label = "", inline = TRUE, status = "info", fill = TRUE,
                  choices = c("Gene Clustering Dendrogram", "TOM Heatmap of Co-expression Gene Groups"), 
                  selected = "Gene Clustering Dendrogram",
                ),
                conditionalPanel(
                  condition = "input.cre_plot_type === 'Gene Clustering Dendrogram'",
                  dropdownButton(
                    sliderInput("creplot_width", "Plot Width:", min = 400, max = 1600, value = 800, step = 10),
                    sliderInput("creplot_height", "Plot Height:", min = 300, max = 1200, value = 400, step = 10),
                    circle = F, status = "danger", icon = icon("cog"), width = "300px", size = 'sm',
                    tooltip = tooltipOptions(title = "Adjust Plot Size")
                  ),
                  tags$div(
                    style = "overflow-x: auto; overflow-y: auto; width: 100%;", 
                    plotOutput("cre_clustPlot", height = "auto", width = "auto")
                    ),
                  shinyjs::hidden(div(id = "cre_download_4", downloadButton("credownload4", "Download Plot")))
                  ),
                conditionalPanel(
                  condition = "input.cre_plot_type === 'TOM Heatmap of Co-expression Gene Groups'",
                  dropdownButton(
                    sliderInput("cre_tom_width", "Plot Width:", min = 400, max = 1500, value = 500, step = 10),
                    sliderInput("cre_tom_height", "Plot Height:", min = 300, max = 1500, value = 500, step = 10),
                    sliderInput("cre_tom_num", "Gene number:", min = 100, max = 10000, value = 600, step = 50),
                    circle = F, status = "danger", icon = icon("cog"), width = "300px", size = 'sm',
                    tooltip = tooltipOptions(title = "Adjust Plot Size")
                    ),
                  tags$div(
                    style = "overflow-x: auto; overflow-y: auto; width: 100%;", 
                    plotOutput("cre_TOMPlot", height = "500px", width = "800px")
                    ),
                  shinyjs::hidden(div(id = "cre_download_5", uiOutput("credownload5")))
                  ),
                ),
              conditionalPanel(
                condition = "input.cre_choose2 === 'table'",
                prettyRadioButtons(
                  inputId = "cre_module", label = "", inline = TRUE, status = "info", fill = TRUE,
                  choices = c("Table" = "table", "Bar Plot" = "plot"), 
                  selected = "table",
                ),
                conditionalPanel(
                  condition = "input.cre_module === 'table'",
                  reactableOutput("module_gene"),
                  hr(),
                  shinyjs::hidden(
                    div(id = "cre_download_7", downloadButton("credownload7", "Download Data"))
                  )
                ),
                conditionalPanel(
                  condition = "input.cre_module === 'plot'",
                  dropdownButton(
                    sliderInput("cre_module_barplot_width", "Plot Width:", min = 500, max = 1000, value = 800, step = 10),
                    sliderInput("cre_module_barplot_height", "Plot Height:", min = 200, max = 1000, value = 500, step = 10),
                    circle = F, status = "danger", icon = icon("cog"), width = "300px", size = 'sm',
                    tooltip = tooltipOptions(title = "Adjust Plot Size")
                  ),
                  tags$div(
                    style = "overflow-x: auto; overflow-y: auto; width: 100%;", 
                    plotOutput("cre_module_barplot", width = "800px", height = "500px"),
                  ),
                  hr(),
                  shinyjs::hidden(
                    div(id = "cre_download_8", uiOutput("credownload8"))
                  )
                )
              )
            )
          )
        ),
        id = "create_tabs"
      )
    )
  )
)




















