# ui
tabItem(
  tabName = 'preprocessing',
  fluidRow(
    box(
      id = 'data_preprocessing',
      height = "auto",
      width = 12,
      elevation = 4,
      title = div(style = "font-size: 20px; font-weight: bold;", "2. Preprocessing"),
      collapsible = TRUE,
      maximizable = TRUE,
      status = "info",
      tabsetPanel(
        tabPanel(
          title = "2a. Quality Control",
          icon = icon("table"),
          fluidRow(
            column(
              width = 3,
              div(style = "margin-bottom: 20px;"),
              
              ## Choose normalized method
              uiOutput("method_ui"),
              uiOutput("cutoff_ui"),
              
              ## Input Sample percentage
              numericInput(
                inputId = "sample_percentage",
                label = tags$div(
                  tags$span(
                    style = "font-size: 18px;", "Fraction of Samples：",
                    "",
                    actionButton("sp", label = "", icon = icon("info"), class = "btn-info btn-xs")
                  ),
                  tags$br(),
                  tags$small(
                    style = "font-size: 13px; color: #666;",
                    "Genes with expression below the cutoff in more than this fraction of samples will be removed (e.g., 0.9 = 90%)."
                  )
                ), 
                min = 0.0, max = 1.0, value = 0.9, step = 0.1
                ),
              ## Choose filter gene number method
              pickerInput(
                inputId = "cutMethod",
                label = tags$div(
                  tags$span(
                    style = "font-size: 18px;", "Method for Selecting Most Variable Genes：",
                    "",
                    actionButton("cm", label = "", icon = icon("info"), class = "btn-info btn-xs"),
                  ),
                  tags$br(),
                  tags$small(
                    style = "font-size: 13px; color: #666;",
                    "Select the metric used to identify highly variable genes."
                  )
                ), 
                choices = list("Median Absolute Deviation (MAD)" = "MAD", "Variance (VAR)" = "Var"),
                options = list(style = "btn-info"),
                width = "100%"
                ),
              ## Input Reserved gene number
              numericInput(
                inputId = "reserved_gene_num",
                label = tags$div(
                  tags$span(
                    style = "font-size: 18px;", "Retain Top-N Variable Genes: ",
                    "",
                    actionButton("rgn",label = "",icon = icon("info"), class = "btn-info btn-xs")
                  ),
                  tags$br(),
                  tags$small(
                    style = "font-size: 13px; color: #666;",
                    "Specify the number of highly variable genes to keep for downstream analysis."
                  )
                ), 
                min = 0, max = 80000, value = 10000, width = '100%'
                ),
              hr(),
              bs4Dash::actionButton(
                "pre_filter",em("Confirm Settings"), status = "primary",
                style = "box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.4); border-radius: 8px;"
              ),
              div(style = "margin-bottom: 20px;"),
              shinyjs::hidden(
                div(id = "step3_1_ui",
                    bs4Dash::actionButton("step3_1",em("STEP 3. Network Construction"),
                                          status = "info", icon = icon("forward"))
                      )
                )
              ),
            column(
              width = 9,
              div(style = "margin-bottom: 20px;"),
              fluidRow(
                uiOutput("pre_choose_ui")
                ),
              fluidRow(
                infoBoxOutput("infobox1"),
                infoBoxOutput("infobox2"),
                infoBoxOutput("infobox3")
                ),
              conditionalPanel(
                condition = "input.pre_choose === 'plot'",
                dropdownButton(
                  sliderInput("preplot_width", "Plot Width:", min = 800, max = 5000, value = 1000, step = 100),
                  sliderInput("preplot_height", "Plot Height:", min = 300, max = 1200, value = 400, step = 10),
                  selectInput("preplot_format", "Image Format for plot:", choices = c("png", "jpeg", "svg")),
                  circle = F, status = "danger", icon = icon("cog"), width = "300px", size = 'sm',
                  tooltip = tooltipOptions(title = "Adjust Plot Size")
                  ),
                tags$div(
                  style = "overflow-x: auto; overflow-y: auto; width: 100%;", 
                  plotlyOutput("pre_clustPlot1", height = "auto", width = "auto")
                  )
                ),
              conditionalPanel(
                condition = "input.pre_choose === 'table'",
                DTOutput("pre_table"),
                shinyjs::hidden(div(id = "predownload2_ui", downloadButton("predownload2", "Download Data")))
                )
              )
            )
          ),
        tabPanel(
          title = "2b. Outlier Detection and Removal",
          icon = icon("trash"),
          value = "outlier_tab",
          fluidRow(
            column(
              width = 3,
              div(style = "margin-bottom: 20px;"),
              selectInput(
                inputId = "outlier",
                label = "Remove the Samples: ",
                choices = NULL,
                multiple = TRUE
              ),
              hr(),
              bs4Dash::actionButton(
                "remove",em("Apply Removal"), status = "primary",
                style = "box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.4); border-radius: 8px;"
              ),
              div(style = "margin-bottom: 20px;"),
              shinyjs::hidden(
                div(id = "step2_2_remove_ui",
                    bs4Dash::actionButton("step2_2_remove",em("STEP 3. Network Construction"),
                                          status = "info", icon = icon("forward"))
                    )
                )
              ),
            column(
              width = 9,
              div(style = "margin-bottom: 20px;"),
              fluidRow(
                uiOutput("pre_choose_remove_ui")
                ),
              conditionalPanel(
                condition = "input.rm_choose === 'plot'  && input.remove >0",
                dropdownButton(
                  sliderInput("rm_width", "Plot Width:", min = 800, max = 5000, value = 1000, step = 100),
                  sliderInput("rm_height", "Plot Height:", min = 300, max = 1200, value = 400, step = 10),
                  selectInput("rm_format", "Image Format for plot:", choices = c("png", "jpeg", "svg")),
                  circle = F, status = "danger", icon = icon("cog"), width = "300px", size = 'sm',
                  tooltip = tooltipOptions(title = "Adjust Plot Size")
                  ),
                tags$div(
                  style = "overflow-x: auto; overflow-y: auto; width: 100%;", 
                  plotlyOutput("pre_clustPlot2",height = "auto", width = "auto")
                )
              ),
              conditionalPanel(
                condition = "input.rm_choose === 'table' ",
                DTOutput("rm_table"),
                shinyjs::hidden(div(id = "pre_download4_ui", downloadButton("predownload4", "Download Data")))
                )
              )
            )
          ),
        id = "pre_tabs"
        )
      )
    )
  )
























