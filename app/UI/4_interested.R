# UI
tabItem(
  tabName = "interested",
  tags$head(tags$style(HTML(paste(
    readLines("www/analysis_style.txt", warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )))),
  div(
    id = "iw-interested", class = "analysis-page",
    div(
      class = "analysis-page-heading",
      tags$h2("4. Module-Trait Association"),
      tags$p("Load trait data, match samples, and explore associations with gene co-expression modules.")
    ),
    box(
      id = "trait_upload", height = "auto", width = 12,
      elevation = 1, maximizable = TRUE, title = "Module-Trait Association",
      collapsible = FALSE, solidHeader = TRUE, status = "info",
      tabsetPanel(
        id = "inter_tabs",
        tabPanel(
          title = "4a. Inspect Trait/Phenotype Data",
          value = "inter_module_trait_tab", icon = icon("upload"),
          fluidRow(
            column(
              width = 4,
              div(
                class = "analysis-settings",
                div(
                  class = "analysis-step",
                  div(
                    class = "analysis-step-heading",
                    span(class = "analysis-step-number", "1"),
                    div(
                      tags$h3("Upload trait data"),
                      tags$p("Load the traits or phenotypes associated with your expression samples.")
                    )
                  ),
                  uiOutput("inter_upload_ui"),
                  tags$p(
                    class = "analysis-hint",
                    "The first column must contain sample names matching the expression matrix. Subsequent columns contain traits. Accepted formats: .txt, .tsv and .csv; maximum file size: 30 MB."
                  ),
                  tags$p(
                    class = "analysis-hint",
                    "For a two-column file, the trait is converted into a 0/1 column for each category. For multiple trait columns, supply numeric values or encode categories as 0/1 columns."
                  )
                ),
                shinyjs::hidden(div(
                  id = "inter_upload_btn_ui",
                  bs4Dash::actionButton(
                    "inter_data_upload", "Load Trait Data", status = "info",
                    icon = icon("upload"), class = "analysis-button analysis-run"
                  )
                )),
                div(
                  id = "inter_demo_ui", class = "analysis-step analysis-next",
                  div(
                    class = "analysis-step-heading",
                    span(class = "analysis-step-number", icon("flask"), `aria-hidden` = "true"),
                    div(
                      tags$h3("Try demo traits"),
                      tags$p("Use the trait dataset for the GSE15459 expression demo.")
                    )
                  ),
                  shinyWidgets::awesomeCheckbox(
                    inputId = "inter_demo", label = "Use Demo Traits",
                    value = FALSE, status = "info"
                  ),
                  tags$p(class = "analysis-hint", "Use with the matching expression demo and its retained samples."),
                  div(
                    class = "analysis-next",
                    downloadButton("inter_download_demo", "Download Demo", class = "analysis-button")
                  )
                ),
                shinyjs::hidden(div(
                  id = "step4b_ui", class = "analysis-next",
                  bs4Dash::actionButton(
                    "step4b", "STEP 4b. Association Heatmap", status = "info",
                    icon = icon("arrow-right"), class = "analysis-button analysis-next-button"
                  )
                ))
              )
            ),
            column(
              width = 8,
              div(
                class = "analysis-results",
                tags$h3(tags$strong("Trait data preview")),
                tags$p(class = "analysis-hint", "Loaded traits are arranged in the same sample order as the expression data used for network construction."),
                div(class = "analysis-plot-scroll", DTOutput("inter_upload_check")),
                shinyjs::hidden(div(
                  id = "inter_download4_ui", class = "analysis-next",
                  downloadButton("inter_download4", "Download Processed Traits")
                )),
                shinyjs::hidden(div(
                  id = "inter_trait_summary_ui", class = "analysis-summary",
                  tags$h3("Trait data summary"), uiOutput("inter_trait_summary")
                ))
              )
            )
          )
        ),
        tabPanel(
          title = "4b. Association Heatmap",
          value = "inter_association_tab", icon = icon("code-merge"),
          fluidRow(
            column(
              width = 12,
              div(
                class = "analysis-step",
                tags$h3("Explore association results"),
                tags$p(
                  class = "analysis-hint",
                  "Switch between the module-trait heatmap and gene-level module membership details."
                ),
                div(
                  style = "margin-top: 16px;",
                  prettyRadioButtons(
                    inputId = "inter_choose", label = NULL,
                    inline = TRUE, status = "info", shape = "round", fill = TRUE,
                    choices = c(
                      "Association Heatmap" = "plot", "Module Details" = "table"
                    ),
                    selected = "plot", width = "100%"
                  )
                )
              ),
              div(
                class = "analysis-results",
                conditionalPanel(
                  condition = "input.inter_choose === 'plot'",
                  dropdownButton(
                    sliderInput("mtplot_width", "Plot Width:", min = 200, max = 1800, value = 800, step = 10),
                    sliderInput("mtplot_height", "Plot Height:", min = 200, max = 1800, value = 700, step = 10),
                    circle = FALSE, status = "info", icon = icon("cog"), width = "300px", size = "sm",
                    tooltip = tooltipOptions(title = "Adjust Plot Size")
                  ),
                  div(class = "analysis-plot-scroll", plotOutput("inter_mtplot", height = "800px")),
                  uiOutput("inter_download")
                ),
                conditionalPanel(
                  condition = "input.inter_choose === 'table'",
                  reactableOutput("KME"),
                  downloadButton("inter_download2", "Download Data")
                ),
                shinyjs::hidden(div(
                  id = "inter_association_summary_ui", class = "analysis-summary",
                  tags$h3("Module-trait association summary"), uiOutput("inter_association_summary")
                )),
                shinyjs::hidden(div(
                  id = "inter_step5_ui", class = "analysis-next",
                  bs4Dash::actionButton(
                    "inter_step5", "STEP 5. Hub Gene Identification", status = "info",
                    icon = icon("arrow-right"), class = "analysis-button analysis-next-button"
                  )
                )),
                shinyjs::hidden(div(
                  id = "inter_step6_ui", class = "analysis-next",
                  bs4Dash::actionButton(
                    "inter_step6", "STEP 6. Functional Enrichment Analysis", status = "info",
                    icon = icon("arrow-right"), class = "analysis-button analysis-next-button"
                  )
                ))
              )
            )
          )
        )
      )
    )
  )
)
