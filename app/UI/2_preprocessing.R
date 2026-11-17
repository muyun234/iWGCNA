# UI
bs4Dash::tabItem(
  tabName = "preprocessing",
  tags$head(tags$style(HTML(paste(
    readLines("www/analysis_style.txt", warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )))),
  div(
    id = "iw-preprocessing", class = "analysis-page",
    div(class = "analysis-page-heading",
        tags$h2("2. Preprocessing"),
        tags$p("Prepare expression data for network analysis.")),
    fluidRow(
      bs4Dash::box(
        id = "data_preprocessing", width = 12, height = "auto",
        status = "info", solidHeader = TRUE, elevation = 1,
        title = "Preprocessing", collapsible = TRUE, maximizable = TRUE,
        tabsetPanel(
          id = "pre_tabs",
          tabPanel(
            title = "2a. Quality Control", icon = icon("sliders-h"),
            fluidRow(
              column(
                width = 4,
                div(class = "analysis-settings",
                    div(class = "analysis-step",
                        div(class = "analysis-step-heading",
                            span(class = "analysis-step-number", "1"),
                            div(tags$h3("Expression transformation"),
                                tags$p("Choose how to transform the uploaded expression values."))),
                        div(class = "analysis-choice-cards", uiOutput("method_ui")),
                        uiOutput("transformation_hint")
                    ),
                    div(class = "analysis-step",
                        div(class = "analysis-step-heading analysis-heading-with-switch",
                            span(class = "analysis-step-number", "2"),
                            div(class = "analysis-heading-text",
                                tags$h3("Low-expression filtering"),
                                tags$p("Remove genes with low expression across samples.")),
                            div(class = "analysis-switch",
                                shinyWidgets::switchInput(
                                  "low_expression_filter", label = NULL,
                                  value = FALSE, onLabel = "On", offLabel = "Off",
                                  onStatus = "info", offStatus = "default", size = "small"
                                ))),
                        div(id = "pre-low-controls", class = "analysis-fields-stack",
                            fluidRow(
                              column(6, uiOutput("cutoff_ui")),
                              column(6, numericInput(
                                "sample_percentage", "Maximum low-expression samples (%)",
                                value = 90, min = 0, max = 100, step = 1, width = "100%"
                              )))),
                        uiOutput("low_filter_hint")
                    ),
                    div(class = "analysis-step",
                        div(class = "analysis-step-heading",
                            span(class = "analysis-step-number", "3"),
                            div(tags$h3("Variable gene selection"),
                                tags$p("Rank genes by variability and choose how many to retain."))),
                        shinyWidgets::pickerInput(
                          "cutMethod", "Ranking method",
                          choices = c("Median absolute deviation (MAD)" = "MAD",
                                      "Variance (VAR)" = "Var"),
                          selected = "MAD", width = "100%",
                          options = shinyWidgets::pickerOptions(style = "btn-outline-info")
                        ),
                        div(
                          class = "analysis-radio-list",
                          shinyWidgets::prettyRadioButtons(
                            "gene_selection", label = NULL,
                            choiceNames = list(
                              tags$span(tags$strong("Keep all available genes"),
                                        tags$small(textOutput("available_gene_count", inline = TRUE))),
                              tags$strong("Keep top N genes")
                            ),
                            choiceValues = c("all", "top"), selected = "all",
                            status = "info", shape = "round", fill = TRUE, width = "100%"
                          )
                        ),
                        conditionalPanel(
                          "input.gene_selection === 'top'",
                          div(class = "analysis-indented-input", numericInput(
                            "reserved_gene_num", "Number of genes to retain",
                            min = 1, max = 100000, value = 10000, step = 100,
                            width = "100%"
                          )))
                    ),
                    bs4Dash::actionButton(
                      "pre_filter", "Apply Preprocessing", status = "info",
                      class = "analysis-button", icon = icon("check")
                    ),
                    shinyjs::hidden(
                      div(
                        id = "step3_1_ui", class = "analysis-next",
                        bs4Dash::actionButton(
                          "step3_1", "STEP 3. Network Construction",
                          status = "info", icon = icon("arrow-right"), class = "analysis-button"
                        )
                      )
                    )                
                )
              ),
              column(
                width = 8,
                div(class = "analysis-results",
                    fluidRow(infoBoxOutput("infobox1"), infoBoxOutput("infobox2"), infoBoxOutput("infobox3")),
                    uiOutput("pre_choose_ui"),
                    conditionalPanel(
                      "input.pre_choose === 'plot'",
                      shinyWidgets::dropdownButton(
                        sliderInput("preplot_width", "Plot Width:", 800, 5000, 1000, step = 100),
                        sliderInput("preplot_height", "Plot Height:", 300, 1200, 400, step = 10),
                        selectInput("preplot_format", "Image Format:", c("png", "jpeg", "svg")),
                        circle = FALSE, status = "info", icon = icon("cog"),
                        width = "300px", size = "sm",
                        tooltip = shinyWidgets::tooltipOptions(title = "Adjust Plot Size")
                      ),
                      div(class = "analysis-plot-scroll",
                          plotlyOutput("pre_clustPlot1", height = "auto", width = "auto"))
                    ),
                    conditionalPanel(
                      "input.pre_choose === 'table'", DTOutput("pre_table"),
                      shinyjs::hidden(div(id = "predownload2_ui", downloadButton("predownload2", "Download Data")))
                    ),
                    shinyjs::hidden(
                      div(
                        id = "pre_summary_ui", class = "analysis-summary", style = "margin-top: 20px;",
                        tags$h3("Preprocessing summary"), uiOutput("preprocessing_summary")
                      )
                    )                
                )
              )
            )
          ),
          tabPanel(
            title = "2b. Outlier Detection and Removal", icon = icon("trash"),
            value = "outlier_tab",
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
                        tags$h3("Outlier sample removal"),
                        tags$p("Select outlier samples identified from the sample clustering plot.")
                      )
                    ),
                    selectInput(
                      inputId = "outlier", label = "Samples to remove",
                      choices = NULL, multiple = TRUE, width = "100%"
                    ),
                    tags$p(
                      class = "analysis-hint",
                      "Selected samples will be excluded from downstream analysis after applying removal."
                    )
                  ),
                  bs4Dash::actionButton(
                    inputId = "remove", label = "Apply Removal", status = "info",
                    icon = icon("check"), class = "analysis-button"
                  ),
                  shinyjs::hidden(
                    div(
                      id = "step2_2_remove_ui", class = "analysis-next",
                      bs4Dash::actionButton(
                        inputId = "step2_2_remove", label = "STEP 3. Network Construction",
                        status = "info", icon = icon("arrow-right"), class = "analysis-button"
                      )
                    )
                  )
                )
              ),
              column(
                width = 8,
                div(
                  class = "analysis-results",
                  uiOutput("pre_choose_remove_ui"),
                  conditionalPanel(
                    condition = "input.rm_choose === 'plot' && input.remove > 0",
                    shinyWidgets::dropdownButton(
                      sliderInput(
                        "rm_width", "Plot Width:", min = 800, max = 5000, value = 1000, step = 100
                      ),
                      sliderInput(
                        "rm_height", "Plot Height:", min = 300, max = 1200, value = 400, step = 10
                      ),
                      selectInput(
                        "rm_format", "Image Format:", choices = c("png", "jpeg", "svg")
                      ),
                      circle = FALSE, status = "info", icon = icon("cog"), width = "300px", size = "sm",
                      tooltip = shinyWidgets::tooltipOptions(title = "Adjust Plot Size")
                    ),
                    div(
                      class = "analysis-plot-scroll",
                      plotlyOutput("pre_clustPlot2", height = "auto", width = "auto")
                    )
                  ),
                  conditionalPanel(
                    condition = "input.rm_choose === 'table'",
                    DTOutput("rm_table"),
                    shinyjs::hidden(
                      div(
                        id = "pre_download4_ui",
                        downloadButton("predownload4", "Download Data")
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )
)
