# UI
tabItem(
  tabName = 'create',
  tags$head(tags$style(HTML(paste(
    readLines("www/analysis_style.txt", warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )))),
  div(
    id = "iw-create", class = "analysis-page",
    div(class = "analysis-page-heading",
        tags$h2("3. Network Construction and Gene Group (Module) Detection"),
        tags$p("Evaluate soft-threshold power, then construct the network and detect gene groups.")),
    fluidRow(
      box(
        id = 'create_network',
        height = "auto",
        width = 12,
        elevation = 1,
        title = "Network Construction and Gene Group Detection",
        collapsible = TRUE,
        solidHeader = TRUE,
        maximizable = TRUE,
        status = "info",
        tabsetPanel(
          tabPanel(
            title = "3a. Soft-Threshold Power Evaluation (β)",
            value = "power_tab",
            icon = icon("table"),
            fluidRow(
              column(
                width = 4,
                div(class = "analysis-settings",
                    div(
                      class = "analysis-step",
                      div(
                        class = "analysis-step-heading",
                        span(class = "analysis-step-number", "1"),
                        div(
                          tags$h3("Network type"), 
                          tags$p(class = "analysis-hint", "Choose how positive and negative correlations contribute to the network.", 
                                 actionButton("network", label = "", icon = icon("info"), class = "btn-info btn-xs",
                                              `aria-label` = "Network type information"))
                        )
                      ),
                      div(
                        class = "analysis-choice-cards",
                        shinyWidgets::prettyRadioButtons(
                          inputId = "Networktype", label = NULL,
                          choiceNames = list(
                            tags$span(tags$strong("Unsigned"), tags$small("Treat positive and negative correlations equally.")),
                            tags$span(tags$strong("Signed"), tags$small("Emphasize positive correlations.")),
                            tags$span(tags$strong("Signed hybrid"), tags$small("Set negative correlations to zero."))
                          ),
                          choiceValues = c("unsigned", "signed", "signed hybrid"),
                          selected = "signed", inline = FALSE, status = "info", shape = "round", fill = TRUE, width = "100%"
                        )
                      ),
                      tags$p("Network type information ")
                    ),
                    div(
                      class = "analysis-step",
                      div(
                        class = "analysis-step-heading",
                        span(class = "analysis-step-number", "2"),
                        div(
                          tags$h3("R² cutoff"), 
                          tags$p(class = "analysis-hint", "Set the target scale-free topology fit for power evaluation.",
                                 actionButton("Cutoff", label = "", icon = icon("info"), class = "btn-info btn-xs",
                                              `aria-label` = "R squared cutoff information"))
                        )
                      ),
                      sliderInput(
                        inputId = "cutoff",
                        label = NULL,
                        min = 0, max = 1, value = 0.9, step = 0.05, width = "100%"
                      ),
                      uiOutput("recommend")
                    ),
                    bs4Dash::actionButton("cre_run", "Evaluate Soft-Threshold Power", status = "info", 
                                          icon = icon("calculator"), class = "analysis-button analysis-run"),
                    shinyjs::hidden(div(
                      id = "step3b_ui", class = "analysis-next", 
                      bs4Dash::actionButton("step3b", "STEP 3b. Network Parameter Settings", 
                                            status = "info", icon = icon("arrow-right"), 
                                            class = "analysis-button analysis-next-button")
                    ))
                )
              ),
              column(
                width = 8,
                div(class = "analysis-results",
                    fluidRow(
                      uiOutput("cre_choose_ui")
                    ),
                    conditionalPanel(
                      condition = "input.cre_choose === 'plot'  && input.cre_run >0",
                      dropdownButton(
                        sliderInput("sftplot_width", "Plot Width:", min = 500, max = 1200, value = 700, step = 10),
                        sliderInput("sftplot_height", "Plot Height:", min = 300, max = 1200, value = 400, step = 10),
                        circle = F, status = "info", icon = icon("cog"), width = "300px", size = 'sm',
                        tooltip = tooltipOptions(title = "Adjust Plot Size")
                      ),
                      tags$div(
                        class = "analysis-plot-scroll",
                        plotOutput("cre_sftplot", height = "auto", width = "auto")
                      ),
                      shinyjs::hidden(div(id = "cre_download", downloadButton("credownload", "Download Plot")))
                    ),
                    conditionalPanel(
                      condition = "input.cre_choose === 'table'",
                      reactableOutput("cre_sfttab"),
                      shinyjs::hidden(div(id = "cre_download_2", downloadButton("credownload2", "Download Data")))
                    ),
                    shinyjs::hidden(div(
                      id = "cre_power_summary_ui", class = "analysis-summary",
                      tags$h3("Soft-threshold evaluation summary"), uiOutput("cre_power_summary")
                    ))
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
                width = 4,
                div(class = "analysis-settings",
                    div(
                      class = "analysis-step",
                      div(class = "analysis-step-heading",
                          span(class = "analysis-step-number", "1"),
                          div(
                            tags$h3("Soft-threshold power (β)"),
                            tags$p("Use the recommended power or select a value based on the evaluation results."))),
                      uiOutput("powerselect_ui"), uiOutput("powerselect_hint")
                    ),
                    div(
                      class = "analysis-step",
                      div(class = "analysis-step-heading",
                          span(class = "analysis-step-number", "2"),
                          div(
                            tags$h3("Minimum gene group size"), 
                            tags$p(class = "analysis-hint", "Set the minimum number of genes required to form a module.", 
                                   actionButton("mms2", label = "", icon = icon("info"), class = "btn-info btn-xs",
                                                `aria-label` = "Minimum gene group size information"))
                          )),
                      sliderInput(
                        inputId = "minMsize",
                        label = NULL,
                        min = 0, max = 200, value = 30, width = "100%"
                      )
                    ),
                    div(
                      class = "analysis-step",
                      div(class = "analysis-step-heading",
                          span(class = "analysis-step-number", "3"),
                          div(
                            tags$h3("Gene group merge cut height"), 
                            tags$p(class = "analysis-hint", "Set the threshold for merging modules with similar expression profiles.", 
                                   actionButton("mch2", label = "", icon = icon("info"), class = "btn-info btn-xs",
                                                `aria-label` = "Merge cut height information")))),
                      sliderInput(
                        inputId = "mch",
                        label = NULL,
                        min = 0, max = 1, value = 0.25, width = "100%"
                      )
                    ),
                    bs4Dash::actionButton(
                      "cre_run2", "Construct Network", status = "info", icon = icon("project-diagram"), 
                      class = "analysis-button analysis-run"
                    ),
                    shinyjs::hidden(div(
                      id = "step4_ui", class = "analysis-next", 
                      bs4Dash::actionButton("step4", "STEP 4. Gene Group-Trait Association", 
                                            status = "info", icon = icon("arrow-right"), 
                                            class = "analysis-button analysis-next-button")
                    )),
                    shinyjs::hidden(div(
                      id = "step6_ui", class = "analysis-next",
                      bs4Dash::actionButton("step6", "STEP 6. Functional Enrichment Analysis", 
                                            status = "info", icon = icon("arrow-right"), 
                                            class = "analysis-button analysis-next-button")
                    ))
                )
              ),
              column(
                width = 8,
                div(class = "analysis-results",
                    fluidRow(
                      uiOutput("cre_choose2_ui")
                    ),
                    conditionalPanel(
                      condition = "input.cre_choose2 === 'plot'  && input.cre_run2 >0",
                      prettyRadioButtons(
                        inputId = "cre_plot_type", label = NULL, inline = TRUE, status = "info", fill = TRUE,
                        choices = c("Gene Clustering Dendrogram", "TOM Heatmap of Co-expression Gene Groups"),
                        selected = "Gene Clustering Dendrogram"
                      ),
                      conditionalPanel(
                        condition = "input.cre_plot_type === 'Gene Clustering Dendrogram'",
                        dropdownButton(
                          sliderInput("creplot_width", "Plot Width:", min = 400, max = 1600, value = 800, step = 10),
                          sliderInput("creplot_height", "Plot Height:", min = 300, max = 1200, value = 400, step = 10),
                          circle = F, status = "info", icon = icon("cog"), width = "300px", size = 'sm',
                          tooltip = tooltipOptions(title = "Adjust Plot Size")
                        ),
                        tags$div(
                          class = "analysis-plot-scroll",
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
                          circle = F, status = "info", icon = icon("cog"), width = "300px", size = 'sm',
                          tooltip = tooltipOptions(title = "Adjust Plot Size")
                        ),
                        tags$div(
                          class = "analysis-plot-scroll",
                          plotOutput("cre_TOMPlot", height = "500px", width = "800px")
                        ),
                        shinyjs::hidden(div(id = "cre_download_5", uiOutput("credownload5")))
                      )
                    ),
                    conditionalPanel(
                      condition = "input.cre_choose2 === 'table'",
                      prettyRadioButtons(
                        inputId = "cre_module", label = "", inline = TRUE, status = "info", fill = TRUE,
                        choices = c("Table" = "table", "Bar Plot" = "plot"),
                        selected = "table"
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
                          circle = F, status = "info", icon = icon("cog"), width = "300px", size = 'sm',
                          tooltip = tooltipOptions(title = "Adjust Plot Size")
                        ),
                        tags$div(
                          class = "analysis-plot-scroll",
                          plotOutput("cre_module_barplot", width = "800px", height = "500px")
                        ),
                        hr(),
                        shinyjs::hidden(
                          div(id = "cre_download_8", uiOutput("credownload8"))
                        )
                      )
                    ),
                    shinyjs::hidden(div(
                      id = "cre_network_summary_ui", class = "analysis-summary",
                      tags$h3("Network construction summary"), uiOutput("cre_network_summary")
                    ))
                )
              )
            )
          ),
          id = "create_tabs"
        )
      )
    )
  )
)
