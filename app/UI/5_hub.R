# ui
tabItem(
  tabName = 'hub',
  fluidRow(
    box(
      id = 'identify_hub_gene',
      height = "auto",
      width = 12,
      elevation = 4,
      title = div(style = "font-size: 20px; font-weight: bold;", "5. Hub Gene Identification"),
      collapsible = TRUE,
      solidHeader = FALSE,
      maximizable = TRUE,
      status = "info",
      fluidRow(
        column(
          width = 3,
          div(style = "margin-bottom: 10px;"),
          selectInput(
            inputId = "hub_module", label = div(style = "font-size: 18px;", "Select a Gene Group: "), 
            choices = NULL, multiple = F
          ),
          selectInput(
            inputId = "hub_trait", label = div(style = "font-size: 18px;", "Select a Trait: "), 
            choices = NULL, multiple = F
          ),
          hr(),
          bs4Dash::actionButton(
            "hub_step6",em("STEP 6. Functional Enrichment Analysis"), 
            status = "info", icon = icon("forward")
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
                fluidRow(
                  column(
                    width = 4,
                    sliderInput(
                      inputId = "hub_MMcut",
                      label = tags$span(
                        style = "font-size: 18px;", "MM Threshold ≧ ","",
                        actionButton("kmecut", label = "", icon = icon("info"), class = "btn-info btn-xs")
                      ), 
                      min = 0, max = 1, step = 0.01, value = 0.5
                    )
                  ),
                  column(
                    width = 4,
                    sliderInput(
                      inputId = "hub_GScut",
                      label = tags$span(
                        style = "font-size: 18px;", "GS Threshold ≧ ","",
                        actionButton("gscut", label = "", icon = icon("info"), class = "btn-info btn-xs")
                      ), 
                      min = 0, max = 1, step = 0.01, value = 0.2
                    )
                  ),
                  column(
                    width = 12,
                    div(
                      id = "hub_gsmm_run_ui",
                      bs4Dash::actionButton(
                        "hub_gsmm_run",em("Generate MM–GS Plot"), status = "primary",
                        style = "box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.4); border-radius: 8px;"
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
                        sliderInput("hub_gsmmplot_width", "Plot Width (px):", min = 200, max = 800, value = 450, step = 10),
                        sliderInput("hub_gsmmplot_height", "Plot Height (px):", min = 200, max = 800, value = 350, step = 10),
                        circle = F, status = "danger", icon = icon("cog"), width = "300px", size = 'sm',
                        tooltip = tooltipOptions(title = "Adjust Plot Size")
                      ),
                      tags$div(
                        style = "overflow-x: auto; overflow-y: auto; width: 100%;", 
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
                        style = "overflow-x: auto; overflow-y: auto", 
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
                  column(
                    width = 4,
                    sliderInput(
                      inputId = "hub_top_cut",
                      label = tags$span(
                        style = "font-size: 18px;", "Top Genes by Connectivity (%)"
                      ), 
                      min = 1, max = 100, step = 0.5, value = 10.0, post = "%"
                    )
                  ),
                  column(
                    width = 12,
                    div(
                      id = "hub_run_ui",
                      bs4Dash::actionButton(
                        "hub_run",em("Generate Hub Gene Network"), status = "primary",
                        style = "box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.4); border-radius: 8px;"
                      )
                    )
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    shinyjs::hidden(div(
                      id = "hub_connect_mm_plot_all_ui",
                      hr(),
                      dropdownButton(
                        sliderInput("hub_connect_mm_plot_width", "Plot Width (px):", min = 200, max = 800, value = 500, step = 10),
                        sliderInput("hub_connect_mm_plot_height", "Plot Height (px):", min = 200, max = 800, value = 350, step = 10),
                        circle = F, status = "danger", icon = icon("cog"), width = "300px", size = 'sm',
                        tooltip = tooltipOptions(title = "Adjust Plot Size")
                      ),
                      tags$div(
                        style = "overflow-x: auto; overflow-y: auto; width: 100%;", 
                        uiOutput("hub_connect_mm_plot_ui")
                      ),
                      div(style = "margin-bottom: 10px;"),
                      downloadButton(outputId = "hub_download6", label = "Download PNG")
                    ))
                  ),
                  column(
                    width = 6,
                    shinyjs::hidden(div(
                      id = "hub_connect_mm_dt_ui",
                      hr(),
                      tags$div(
                        style = "overflow-x: auto; overflow-y: auto", 
                        DTOutput("hub_connect_mm_dt")
                      ),
                      div(style = "margin-bottom: 10px;"),
                      downloadButton("hub_download7", "Download Data")
                    ))
                  )
                ),
                fluidRow(
                  column(
                    width = 12,
                    hr(),
                    prettyRadioButtons(
                      inputId = "hub_network_choose", label = "", inline = TRUE, status = "info", fill = TRUE,
                      choices = c("WGCNA network" = "wgcna_network", "Connectivity vs PPI Degree" = "connect_ppidegree"),
                      selected = "wgcna_network"
                    ),
                    conditionalPanel(
                      condition = "input.hub_network_choose === 'wgcna_network'",
                      fluidRow(
                        column(
                          width = 12,
                          div(
                            style = "display: flex; gap: 15px; overflow-x: auto; padding-bottom: 10px; margin-left: 10px",
                            div(style = "min-width: 250px;", sliderInput(
                              "hub_wgcna_hub_num", "Number of Top Hub Genes to Label: ",
                              min = 1, max = 50, value = 10, step = 1
                            ))
                          )
                        )
                      ),
                      fluidRow(
                        column(
                          width = 12,
                          shinyjs::hidden(div(
                            id = "hub_wgcna_plot_ui", style = "position: relative;",
                            dropdownButton(
                              sliderInput("hub_wgcna_width", "Plot Width (px):", min = 500, max = 2000, value = 600, step = 50),
                              sliderInput("hub_wgcna_height", "Plot Height (px):", min = 500, max = 2000, value = 600, step = 50),
                              sliderInput("hub_wgcna_tex_size", "Text size:", min = 1, max = 8, value = 4, step = 0.1),
                              # sliderInput("hub_wgcna_edge_width", "Edge width:", min = 1, max = 5, value = 3, step = 0.5),
                              circle = F, status = "danger", icon = icon("cog"), width = "300px", size = 'sm',
                              tooltip = tooltipOptions(title = "Adjust Plot Setting")
                            ),
                            div(
                              style = "overflow-x: auto; overflow-y: auto; width: 100%;", 
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
                          width = 12,
                          div(
                            prettyRadioButtons(
                              inputId = "hub_ppi_level",
                              label = tagList(tags$div("PPI network confidence level: ")),
                              choices = c(
                                "Low (STRING score ≥ 400)" = 400, 
                                "Medium (STRING score ≥ 700)" = 700,
                                "High (STRING score ≥ 900)" = 900
                              ), 
                              selected = 400, inline = TRUE, status = "info", fill = TRUE
                            )
                          ),
                          shinyjs::hidden(div(
                            id = "hub_ppi_update_ui",
                            bs4Dash::actionButton(
                              "hub_ppi_update",em("Update"), status = "primary",
                              style = "box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.4); border-radius: 8px;"
                            ),
                            div(style = "margin-bottom: 10px;")
                          ))
                        )
                      ),
                      fluidRow(
                        column(
                          width = 12,
                          shinyjs::hidden(div(
                            id = "hub_connect_ppi_plot_ui",
                            dropdownButton(
                              sliderInput("hub_connect_ppi_width", "Plot Width (px):", min = 500, max = 1000, value = 850, step = 50),
                              sliderInput("hub_connect_ppi_height", "Plot Height (px):", min = 300, max = 1000, value = 450, step = 50),
                              sliderInput("hub_connect_ppi_text", "Label Size:", min = 1, max = 10, value = 4, step = 1),
                              circle = F, status = "danger", icon = icon("cog"), width = "300px", size = 'sm',
                              tooltip = tooltipOptions(title = "Adjust Plot Setting")
                            ),
                            div(
                              style = "overflow-x: auto; overflow-y: auto; width: 100%;", 
                              uiOutput("hub_connect_ppi_plot_ui")
                            ),                        
                            div(style = "margin-bottom: 10px;"),
                            uiOutput("hub_download5")
                          ))
                        ),
                        column(
                          width = 12,
                          shinyjs::hidden(div(
                            id = "hub_connect_ppi_dt_ui",
                            hr(),
                            tags$div(
                              style = "overflow-x: auto; overflow-y: auto", 
                              DTOutput("hub_connect_ppi_dt")
                            ),
                            div(style = "margin-bottom: 10px;"),
                            downloadButton("hub_download8", "Download Data")
                          ))
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









