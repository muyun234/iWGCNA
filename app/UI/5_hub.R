# ui
tabItem(
  tabName = 'hub',
  tags$head(tags$style(HTML(paste(
    readLines("www/analysis_style.txt", warn = FALSE, encoding = "UTF-8"), collapse = "\n"
  )))),
  div(
    id = "iw-hub", class = "analysis-page",
    div(class = "analysis-page-heading",
        tags$h2("5. Hub Gene Identification"),
        tags$p("Explore module membership, trait significance and hub gene networks.")),
    fluidRow(
      box(
        id = 'identify_hub_gene',
        height = "auto",
        width = 12,
        elevation = 1,
        title = "Hub Gene Identification",
        collapsible = TRUE,
        solidHeader = TRUE,
        maximizable = TRUE,
        status = "info",
        fluidRow(
          column(
            width = 3,
            div(class = "analysis-settings",
                div(
                  class = "analysis-step",
                  div(class = "analysis-step-heading",
                      span(class = "analysis-step-number", "1"),
                      div(tags$h3("Module and trait"), tags$p("Choose the module and trait to explore."))),
                  selectInput(
                    inputId = "hub_module", label = "Module", choices = NULL,
                    multiple = FALSE, width = "100%"
                  ),
                  selectInput(
                    inputId = "hub_trait", label = "Trait", choices = NULL,
                    multiple = FALSE, width = "100%"
                  )
                ),
                uiOutput("hub_sidebar_summary"),
                shinyjs::hidden(div(
                  id = "hub_step6_ui", class = "analysis-next",
                  bs4Dash::actionButton(
                    "hub_step6", "STEP 6. Functional Enrichment Analysis",
                    status = "info", icon = icon("arrow-right"), class = "analysis-button analysis-next-button"
                  )
                ))
            )
          ),
          column(
            width = 9,
            div(style = "margin-bottom: 10px;"),
            div(
              id = "hub_choose_ui",
              tabsetPanel(
                id = "hub_choose",
                selected = "mm_gs",
                tabPanel(
                  title = "Module Membership vs Gene Significance",
                  value = "mm_gs",
                  div(style = "margin-bottom: 10px;"),
                  div(
                    class = "analysis-step",
                    div(class = "analysis-step-heading",
                        span(class = "analysis-step-number", "2"),
                        div(tags$h3("Gene selection thresholds"), tags$p("Highlight genes meeting both absolute module membership and gene significance thresholds."))),
                    fluidRow(
                      column(
                        width = 6,
                        sliderInput(
                          inputId = "hub_MMcut",
                          label = tags$span(
                            style = "font-size: 18px;", "|MM| threshold ≥ ","",
                            actionButton("kmecut", label = "", icon = icon("info"), class = "btn-info btn-xs")
                          ),
                          min = 0, max = 1, step = 0.01, value = 0.5
                        )
                      ),
                      column(
                        width = 6,
                        sliderInput(
                          inputId = "hub_GScut",
                          label = tags$span(
                            style = "font-size: 18px;", "|GS| threshold ≥ ","",
                            actionButton("gscut", label = "", icon = icon("info"), class = "btn-info btn-xs")
                          ),
                          min = 0, max = 1, step = 0.01, value = 0.2
                        )
                      )
                    )
                  ),
                  fluidRow(
                    column(
                      width = 12,
                      div(
                        id = "hub_gsmm_run_ui",
                        bs4Dash::actionButton(
                          "hub_gsmm_run",em("Generate MM–GS Plot"), status = "info",
                          class = "analysis-button analysis-run", icon = icon("check")
                        )
                      )
                    )
                  ),
                  fluidRow(
                    column(
                      width = 12,
                      uiOutput("hub_gsmm_text")
                    )
                  ),
                  fluidRow(
                    column(
                      width = 12,
                      hr(),
                      shinyjs::hidden(div(
                        id = "hub_gs2mm_plot_all_ui",
                        dropdownButton(
                          sliderInput("hub_gsmmplot_width", "Plot Width (px):", min = 200, max = 800, value = 580, step = 10),
                          sliderInput("hub_gsmmplot_height", "Plot Height (px):", min = 200, max = 800, value = 460, step = 10),
                          circle = F, status = "info", icon = icon("cog"), width = "300px", size = 'sm',
                          tooltip = tooltipOptions(title = "Adjust Plot Size")
                        ),
                        tags$div(
                          class = "analysis-plot-scroll",
                          uiOutput("hub_gs2mm_plot_ui")
                        ),
                        div(style = "margin-bottom: 10px;"),
                        uiOutput("hub_download")
                      ))
                    ),
                    column(
                      width = 12,
                      shinyjs::hidden(div(
                        id = "hub_gs2mmdt_ui",
                        hr(),
                        tags$div(
                          class = "analysis-plot-scroll",
                          DTOutput("hub_kMEhub")
                        ),
                        div(style = "margin-bottom: 10px;"),
                        downloadButton("hub_download2", "Download Data")
                      ))
                    )
                  )
                ),
                tabPanel(
                  title = "Hub Gene Networks",
                  value = "hub",
                  div(style = "margin-bottom: 10px;"),
                  fluidRow(
                    style = "display: flex; flex-wrap: wrap; align-items: flex-end;",
                    column(
                      width = 7,
                      div(
                        class = "analysis-step",
                        div(
                          class = "analysis-step-heading", span(class = "analysis-step-number", "2"),
                          div(
                            tags$h3("Connectivity-based selection"),
                            tags$p("Select the top fraction of module genes by intramodular connectivity.")
                          )
                        ),
                        sliderInput(
                          inputId = "hub_top_cut", label = NULL,
                          min = 1, max = 100, step = 0.5, value = 10.0, post = "%"
                        )
                      )
                    ),
                    
                    column(
                      width = 5,
                      div(
                        id = "hub_run_ui", style = "margin-bottom: 18px;",
                        bs4Dash::actionButton(
                          "hub_run", "Generate Hub Gene Network",
                          status = "info", class = "analysis-button analysis-run", icon = icon("check")
                        )
                      )
                    )
                  ),
                  fluidRow(
                    column(
                      width = 12,
                      hr(),
                      prettyRadioButtons(
                        inputId = "hub_connect_mm_choose", label = NULL, inline = TRUE,
                        status = "info", shape = "round", fill = TRUE,
                        choices = c("Connectivity vs MM" = "plot", "Gene Details" = "table"),
                        selected = "plot", width = "100%"
                      ),
                      conditionalPanel(
                        condition = "input.hub_connect_mm_choose === 'plot'",
                        shinyjs::hidden(div(
                          id = "hub_connect_mm_plot_all_ui",
                          dropdownButton(
                            sliderInput("hub_connect_mm_plot_width", "Plot Width (px):", min = 200, max = 1000, value = 850, step = 10),
                            sliderInput("hub_connect_mm_plot_height", "Plot Height (px):", min = 200, max = 800, value = 380, step = 10),
                            circle = FALSE, status = "info", icon = icon("cog"), width = "300px", size = "sm",
                            tooltip = tooltipOptions(title = "Adjust Plot Size")
                          ),
                          div(class = "analysis-plot-scroll", uiOutput("hub_connect_mm_plot_ui")),
                          div(style = "margin-bottom: 10px;"),
                          downloadButton("hub_download6", "Download PNG")
                        ))
                      ),
                      conditionalPanel(
                        condition = "input.hub_connect_mm_choose === 'table'",
                        shinyjs::hidden(div(
                          id = "hub_connect_mm_dt_ui",
                          div(class = "analysis-plot-scroll", DTOutput("hub_connect_mm_dt")),
                          div(style = "margin-bottom: 10px;"),
                          downloadButton("hub_download7", "Download Data")
                        ))
                      )
                    )
                  ),
                  fluidRow(
                    column(
                      width = 12,
                      hr(),
                      div(
                        class = "analysis-step",
                        div(class = "analysis-step-heading",
                            span(class = "analysis-step-number", "3"),
                            div(tags$h3("Network view"), tags$p("Choose the WGCNA network or compare connectivity with PPI degree."))),
                        prettyRadioButtons(
                          inputId = "hub_network_choose", label = NULL, inline = TRUE, status = "info", fill = TRUE,
                          choices = c("WGCNA network" = "wgcna_network", "Connectivity vs PPI Degree" = "connect_ppidegree"),
                          selected = "wgcna_network"
                        )
                      ),
                      conditionalPanel(
                        condition = "input.hub_network_choose === 'wgcna_network'",
                        fluidRow(
                          column(
                            width = 12,
                            shinyjs::hidden(div(
                              id = "hub_wgcna_plot_ui", style = "position: relative;",
                              dropdownButton(
                                sliderInput("hub_wgcna_width", "Plot Width (px):", min = 500, max = 2000, value = 600, step = 50),
                                sliderInput("hub_wgcna_height", "Plot Height (px):", min = 500, max = 2000, value = 600, step = 50),
                                sliderInput("hub_wgcna_tex_size", "Text size:", min = 1, max = 8, value = 4, step = 0.1),
                                sliderInput("hub_wgcna_hub_num", "Top hub genes to label:", min = 1, max = 50, value = 10, step = 1),
                                # sliderInput("hub_wgcna_edge_width", "Edge width:", min = 1, max = 5, value = 3, step = 0.5),
                                circle = F, status = "info", icon = icon("cog"), width = "300px", size = 'sm',
                                tooltip = tooltipOptions(title = "Adjust Plot Setting")
                              ),
                              div(
                                class = "analysis-plot-scroll",
                                uiOutput("hub_wgcna_network_ui")
                              ),
                              div(style = "margin-bottom: 10px;"),
                              downloadButton("hub_download3", "Download Data")
                            ))
                          )
                        )
                      ),
                      conditionalPanel(
                        condition = "input.hub_network_choose === 'connect_ppidegree'",
                        fluidRow(
                          column(
                            width = 7,
                            div(
                              class = "analysis-step",
                              div(
                                class = "analysis-step-heading",
                                span(class = "analysis-step-number", "4"),
                                div(
                                  tags$h3("PPI settings"),
                                  tags$p("Choose the STRING database version and interaction confidence threshold.")
                                )
                              ),
                              div(
                                id = "hub_db_type_ui",
                                selectInput(
                                  "hub_string_version", "STRING database version:",
                                  choices = c("11.5", "12.0", "12.5"),
                                  selected = "12.5", multiple = FALSE, width = "100%"
                                )
                              ),
                              selectInput(
                                inputId = "hub_ppi_level", label = "PPI network confidence level:",
                                choices = c(
                                  "Low (STRING score ≥ 400)" = "400",
                                  "Medium (STRING score ≥ 700)" = "700",
                                  "High (STRING score ≥ 900)" = "900"
                                ),
                                selected = "900", width = "100%"
                              )
                            )
                          ),
                          column(
                            width = 5,
                            div(
                              class = "analysis-step",
                              div(
                                class = "analysis-step-heading",
                                span(class = "analysis-step-number", "5"),
                                div(
                                  tags$h3("PPI degree selection"),
                                  tags$p("Choose the top fraction of genes by PPI degree.")
                                )
                              ),
                              sliderInput(
                                inputId = "hub_ppi_degree_percent", label = NULL,
                                min = 1, max = 100, value = 10, step = 0.5, post = "%"
                              )
                            )
                          )
                        ),
                        fluidRow(
                          column(
                            width = 12,
                            shinyjs::hidden(
                              div(
                                id = "hub_connect_ppi_plot_all_ui",
                                dropdownButton(
                                  sliderInput("hub_connect_ppi_width", "Plot Width (px):", min = 500, max = 1000, value = 850, step = 50),
                                  sliderInput("hub_connect_ppi_height", "Plot Height (px):", min = 300, max = 1000, value = 450, step = 50),
                                  sliderInput("hub_connect_ppi_text", "Label Size:", min = 1, max = 10, value = 4, step = 1),
                                  circle = FALSE, status = "info", icon = icon("cog"), width = "300px",
                                  size = "sm", tooltip = tooltipOptions(title = "Adjust Plot Setting")
                                ),
                                div(class = "analysis-plot-scroll", uiOutput("hub_connect_ppi_plot_ui")),
                                div(style = "margin-bottom: 10px;"), 
                                uiOutput("hub_download5")
                              )
                            )
                          ),
                          column(
                            width = 12,
                            shinyjs::hidden(
                              div(
                                id = "hub_connect_ppi_dt_ui",
                                hr(),
                                div(class = "analysis-plot-scroll", DTOutput("hub_connect_ppi_dt")),
                                div(style = "margin-bottom: 10px;"),
                                downloadButton("hub_download8", "Download Data")
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
      )
    )
  )
)




