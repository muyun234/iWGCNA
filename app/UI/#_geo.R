# ui
tabItem(
  tabName = "geo",
  fluidRow(
    box(
      id = "GEO", height = "auto", width = 12, elevation = 4,
      title = div(style = "font-size: 20px; font-weight: bold;", "GEO Data Retrieval"),
      collapsible = FALSE, solidHeader = FALSE, status = "info",
      tabsetPanel(
        id = "geo_tabs",
        tabPanel(
          title = "Download Expression Data",
          value = "expr_tab",
          icon = icon("table"),
          fluidRow(
            column(
              width = 3,
              div(style = "margin-bottom: 20px;"),
              textInput("geo_num", div(style = "font-size: 18px;", "Please Enter the GEO Accession (ex: GSE137342): ")),
              uiOutput("geo_error", style = "white-space: pre-wrap; color: red;"),
              bs4Dash::actionButton(
                "geo_update",em("Run"), status = "primary",
                style = "box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.4); border-radius: 8px;"
              ),
              hr(),
              uiOutput("geo_platform_ui"),
              shinyjs::hidden(
                div(id = "geo_nor_ui", checkboxInput("geo_nor", span(style = "font-size: 16px;", "Perform Normalization"), value = FALSE)
                ),
                tags$small(style = "color: #555;","Only for single-channel data.")
              ),
              shinyjs::hidden(
                div(id = "geo_nor_setting_ui",
                    selectInput(
                      "geo_nor_method", 
                      label = tags$span(
                        "Please Select the Method: ", "",
                        actionButton("geo_nor_method_info",label = "",icon = icon("info"), class = "btn-info btn-xs")
                      ), 
                      choices = c("quantile" = "quantile", "scale" = "scale", "cyclicloess" = "cyclicloess"), 
                      selected = "quantile"
                    )
                )
              ),
              shinyjs::hidden(
                div(id = "geo_platform_update_ui", 
                    bs4Dash::actionButton(
                      "geo_platform_update",em("Run"), status = "primary",
                      style = "box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.4); border-radius: 8px;"
                    )
                )
              ),
              shinyjs::hidden(
                div(id = "geo_next_step_ui", hr(),
                    bs4Dash::actionButton("geo_next_step",em("Next Step: Download Trait Data"), 
                                          status = "info", icon = icon("forward"))
                )
              )
            ),
            column(
              width = 9,
              div(style = "margin-bottom: 20px;"),
              uiOutput("geo_choose_ui"),
              conditionalPanel(
                condition = 'input.geo_choose === "table"',
                DTOutput("geo_exp"),
                br(),
                downloadButton("geo_download", "Download Data")
              ),
              conditionalPanel(
                condition = 'input.geo_choose === "box_plot"',
                div(style = "margin-bottom: 20px;"),
                plotOutput("geo_boxplot"),
                br(),
                uiOutput("geo_download2")
              )
            )
          )
        ),
        tabPanel(
          title = "Download Trait Data",
          value = "trait_tab",
          icon = icon("tag"),
          fluidRow(
            column(
              width = 3,
              div(style = "margin-bottom: 20px;"),
              div("Original phenotypic data from GEO", style = "font-size:18px; font-weight:bold;"),
              div("Please select the columns you want to keep. 
                  It is strongly recommended to retain the column containing GEO sample IDs, 
                  usually named 'geo_accession'.",
                  style = "color:#555; font-size:15px;"),
              br(),
              bs4Dash::actionButton(
                "geo_generate_filtered_data",em("Generate Filtered Data"), status = "primary",
                style = "box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.4); border-radius: 8px;"
              )
            ),
            column(
              width = 9,
              div(style = "margin-bottom: 20px;"),
              DTOutput("geo_phen"),
              shinyjs::hidden(div(id = "geo_download3_ui", downloadButton("geo_download3", "Download Data")))
            )
          ),
          hr(),
          fluidRow(
            column(
              width = 3,
              div(style = "margin-bottom: 20px;"),
              div("Table Based on Selected Columns", style = "font-size:18px; font-weight:bold;"),
              br(),
              bs4Dash::actionButton(
                "geo_convert_traits",em("Convert into 0/1 Encoding"), status = "primary",
                style = "box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.4); border-radius: 8px;"
              )
            ),
            column(
              width = 9,
              div(style = "margin-bottom: 20px;"),
              DTOutput("geo_trait"),
              shinyjs::hidden(div(id = "geo_download4_ui", downloadButton("geo_download4", "Download Data")))
            )
          ),
          hr(),
          fluidRow(
            column(
              width = 3,
              div(style = "margin-bottom: 20px;"),
              div("Table Transformed into 0/1 Encoding", style = "font-size:18px; font-weight:bold;"),
              uiOutput("geo_rm_sample_ui"),
              bs4Dash::actionButton(
                "geo_rm_update",em("Update"), status = "primary",
                style = "box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.4); border-radius: 8px;"
              ),
            ),
            column(
              width = 9,
              div(style = "margin-bottom: 20px;"),
              DTOutput("geo_transformed_trait"),
              shinyjs::hidden(div(id = "geo_download5_ui", downloadButton("geo_download5", "Download Data")))
            )
          )
        )
      )
    )
  )
)