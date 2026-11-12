# ui
tabItem(
  tabName = "interested",
  box(
    id = "trait_upload", height = "auto", width = 12, elevation = 4, maximizable = TRUE,
    title = div(style = "font-size: 20px; font-weight: bold;", "4. Module-Trait Association"),
    collapsible = FALSE, solidHeader = FALSE, status = "info",
    tabsetPanel(
      tabPanel(
        title = "4a. Inspect Trait/Phenotype Data",
        value = "inter_module_trait_tab",
        icon = icon("upload"),
        fluidRow(
          column(
            width = 4,
            div(style = "margin-bottom: 20px;"),
            uiOutput("inter_upload_ui"),
            shinyjs::hidden(
              div(
                id = "inter_upload_btn_ui",
                bs4Dash::actionButton(
                  "inter_data_upload",em("Load Table"), status = "primary",
                  style = "box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.4); border-radius: 8px;"
                )
              )
            ),
            hr(),
            div(
              id = "inter_demo_ui", 
              tags$p(style = "font-size: 18px;", strong("Choose a Demo Trait Dataset")),
              awesomeCheckbox(inputId = "inter_demo", label = span(style = "font-size: 18px;", "Demo Trait 1"), 
                              value = FALSE, status = "danger"),
              downloadButton("inter_download_demo", "Download"),
              hr()
            ),
            shinyjs::hidden(div(
              id = "step4b_ui", 
              bs4Dash::actionButton("step4b", em("STEP 4b. Association Heatmap"), 
                                    status = "info", icon = icon("forward"))
            ))
          ),
          column(
            width = 8,
            div(style = "margin-bottom: 20px;"),
            tags$div(
              style = "overflow-x: auto; overflow-y: auto; width: 100%;", 
              DTOutput("inter_upload_check")
            ),
            shinyjs::hidden(div(
              id = "inter_download4_ui",
              downloadButton("inter_download4", "Download"),
            ))
          )
        ),
      ),
      # calaulate association
      tabPanel(
        title = "4b. Association Heatmap",
        value = "inter_association_tab",
        icon = icon("code-merge"),
        fluidRow(
          column(
            width = 12,
            div(style = "margin-bottom: 20px;"),
            prettyRadioButtons(
              inputId = "inter_choose", label = "", inline = TRUE, status = "info", fill = TRUE,
              choices = c("Association Heatmap" = "plot", "Module Details" = "table"),
              selected = "plot"
            ),
            conditionalPanel(
              condition = "input.inter_choose === 'plot'",
              dropdownButton(
                sliderInput("mtplot_width", "Plot Width:", min = 200, max = 1800, value = 800, step = 10),
                sliderInput("mtplot_height", "Plot Height:", min = 200, max = 1800, value = 700, step = 10),
                circle = F, status = "danger", icon = icon("cog"), width = "300px", size = 'sm',
                tooltip = tooltipOptions(title = "Adjust Plot Size")
              ),
              tags$div(
                style = "overflow-x: auto; overflow-y: auto; width: 100%;", 
                plotOutput("inter_mtplot", height = "800px")
              ),
              uiOutput("inter_download")
            ),
            conditionalPanel(
              condition = "input.inter_choose === 'table'",
              reactableOutput("KME"),
              downloadButton("inter_download2", "Download Data")
            )
          ),
          column(
            width = 12,
            hr(),
            shinyjs::hidden(div(
                id = "inter_step5_ui", 
                bs4Dash::actionButton("inter_step5",em("STEP 5. Hub Gene Identification"), 
                                      status = "info", icon = icon("forward"))
            )),
            shinyjs::hidden(div(
              id = "inter_step6_ui", 
              div(style = "margin-bottom: 10px;"),
              bs4Dash::actionButton("inter_step6",em("STEP 6 Functional Enrichment Analysis"),
                                    status = "info", icon = icon("forward"))
            ))
          )
        )
      ),
      id = "inter_tabs"
      )
    )
  )
