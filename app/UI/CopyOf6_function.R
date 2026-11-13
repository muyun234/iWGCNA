# ui
tabItem(
  tabName = 'function',
  fluidRow(
    box(
      id = 'functional_analysis',
      height = "auto",
      width = 12,
      elevation = 4,
      title = div(style = "font-size: 20px; font-weight: bold;", "6. Functional Enrichment Analysis"),
      collapsible = TRUE,
      solidHeader = FALSE,
      maximizable = TRUE,
      status = "info",
      fluidRow(
        column(
          width = 3,
          # div(style = "margin-bottom: 20px;"),
          selectInput("func_module_select", label = div(style = "font-size: 18px;", "Select one or more Gene Groups: "), 
                      choices = NULL, multiple = TRUE),
          hr(),
          div("Functional Enrichment Analysis Options: ", style = "font-size:20px; font-weight:bold;"),
          div(style = "margin-bottom: 10px;"),
          # div("Enrichment Databases", style = "font-size:19px; font-weight:bold;"),
          box(
            title = div("Enrichment Databases", style = "font-size:19px; font-weight:bold;"),
            height = "auto", width = 12, elevation = 1, collapsible = TRUE, solidHeader = TRUE, status = "lightblue", 
            div(
              id = "func_db_type_ui",
              prettyRadioButtons(
                inputId = "func_db_type", 
                label = tagList(
                  tags$div(style = "font-size: 18px; font-weight:bold;", "Enrichment Analysis Tool")
                ),
                choices = c("STRING" = "stringdb", "ClusterProfiler" = "clusterprofiler"), 
                selected = "clusterprofiler", inline = TRUE, status = "info", fill = TRUE
              ),
              conditionalPanel(
                condition = "input.func_db_type === 'stringdb'",
                selectInput(
                  "func_string_version", label = div(style = "font-size: 18px;", "Select STRING Database Version: "),
                  choices = c("11.5", "12.0", "12.5"), selected = "12.5",multiple = FALSE
                )
              )
            ),
            div("Gene Ontology (GO) Categories", style = "font-size:18px; font-weight:bold;"),
            div(
              style = "margin-left: 20px;",
              pickerInput(
                inputId = "func_go_category",
                label = tags$span(style = "font-size: 16px;", "Select GO Category: "),
                choices = c(
                  "None" = "",
                  "Biological Process" = "BP",
                  "Molecular Function" = "MF",
                  "Cellular Component" = "CC",
                  "All (multi-select) " = "all"
                ),
                selected = "BP", multiple = FALSE, options = pickerOptions(clearable = TRUE)
              ),
              shinyjs::hidden(
                div(id = "func_go_setting_ui",
                    # prettyRadioButtons(
                    #   inputId = "func_go_db",
                    #   label = tags$span(
                    #     style = "font-size: 16px;", "Select GO Database:"
                    #     ),
                    #   choices = (c("org.Hs.eg.db (Bioconductor)" = "orgdb", "DAVID" = "david")),
                    #   selected = "orgdb",
                    #   inline = TRUE, status = "info", fill = TRUE
                    #   ),
                    numericInput(
                      inputId = "func_go_padj", label = div(style = "font-size: 16px;", "P-Adjust Threshold: "),
                      value = 0.05, min = 0, max = 1, step = 0.01,
                    ),
                    numericInput(
                      inputId = "func_go_topn", label = div(style = "font-size: 16px;", "Number of GO Terms per Gene Group: "),
                      value = 3, min = 1, step = 1
                    )
                )
              )
            ),
            div(
              id = "func_kegg_ui",
              style = "margin-left: 20px;",
              hr(style = "border-top: 1px dashed #999;"),
              div("Pathway Databases", style = "font-size:18px; font-weight:bold;"),
              prettyCheckbox(
                inputId = "func_kegg", 
                label = tags$span(style = "font-size: 18px;", "KEGG Pathway Analysis"),
                value = FALSE, status = "info", shape = "curve", outline = TRUE
              ),
              # tags$small(
              #   style = "color: #555;",
              #   "⚠️ Only supported with org.Hs.eg.db (Bioconductor database)"
              #   ),
              shinyjs::hidden(
                div(id = "func_kegg_setting_ui",
                    numericInput(
                      inputId = "func_kegg_padj", 
                      label = div(style = "font-size: 16px;", "P-Adjust Threshold: "),
                      value = 0.05, min = 0, max = 1, step = 0.01,
                    ),
                    numericInput(
                      inputId = "func_kegg_topn", label = div(style = "font-size: 16px;", "Number of KEGG Pathways per Gene Group: "),
                      value = 3, min = 1, step = 1
                    )
                )
              )
            )
          ),
          hr(),
          selectInput("func_genelist", label = tags$div(style = "font-size: 18px;", "Specify genes to highlight on the heatmap (optional): "),
                      choices = NULL, selected = NULL, multiple = TRUE
          ),
          bs4Dash::actionButton(
            "func_update",em("Run"), status = "primary",
            style = "box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.4); border-radius: 8px;"
            )
          ),
        column(
          width = 9,
          div(style = "margin-bottom: 20px;"),
          fluidRow(
            uiOutput("func_choose_ui")
            ),
          conditionalPanel(
            condition = "input.func_choose === 'heatmap'  && input.func_update >0",
            box(
              height = "auto", width = 12, elevation = 2, title = "Combined Heatmap and Annotation",
              collapsible = FALSE, solidHeader = FALSE, maximizable = TRUE, status = "info",
              dropdownButton(
                sliderInput("func_heatmap_width", "Plot Width:", min = 300, max = 3000, value = 1200, step = 100),
                sliderInput("func_heatmap_height", "Plot Height:", min = 300, max = 3000, value = 1500, step = 100),
                sliderInput("func_heatmap_anno_size", "Text Size:", min = 1, max = 20, value = 12, step = 1),
                radioButtons("func_heatmap_box_style", "Text Box Style:", 
                             choices = c("Link" = "link", "Block" = "block"),
                             selected = "link", inline = TRUE
                ),
                # radioButtons("func_heatmap_type", "Plot type:", choices = c("Heatmap" = "heatmap", "Heatmap & Line" = "both"), selected = "both", inline = TRUE),
                div(id = "func_heatmap_show_mode_ui", 
                    radioButtons(
                      "func_heatmap_show_mode", "Annotation display mode:",
                      choices = c("Together" = "together", "Separate" = "separate"),
                      selected = "together", inline = TRUE
                    )
                ),
                # selectInput("func_heatmap_module_color", "Gene groups with black annotations: ", choices = NULL, multiple = TRUE),
                awesomeCheckbox("func_heatmap_module_color", "Chage annotation color with black", value = FALSE),
                circle = F, status = "danger", icon = icon("cog"), width = "300px", size = 'sm',
                tooltip = tooltipOptions(title = "Adjust Plot Size")
              ),
              tags$div(
                style = "overflow-x: auto; overflow-y: auto; width: 100%;", 
                uiOutput("func_heatmap_ui")
              ),
              br(),
              shinyjs::hidden(div(id = "func_download_ui", downloadButton("func_download", "Download Plot")))
            )
          ),
          conditionalPanel(
            condition = "input.func_choose === 'bar_plot'",
            dropdownButton(
              sliderInput("func_bar_width", "Plot Width:", min = 500, max = 1500, value = 800, step = 10),
              sliderInput("func_bar_height", "Plot Height:", min = 200, max = 1500, value = 500, step = 10),
              sliderInput("func_bar_x_expand", "Extend X-axis Range:", min = 0, max = 3, value = 0.2, step = 0.05),
              sliderInput("func_bar_text_size", "Text size on bars:", min = 1, max = 10, value = 3, step = 1),
              numericInput("func_bar_facet_ncol", "Plots per Row:", value = 3, min = 1, max = 6, step = 1),
              shinyjs::hidden(
                div(id = "func_bar_show_mode_ui", 
                    radioButtons(
                      "func_bar_show_mode", "Annotation display mode:",
                      choices = c("Together" = "together", "Separate" = "separate"),
                      selected = "separate", inline = TRUE
                      )
                    )
                ),
              circle = F, status = "danger", icon = icon("cog"), width = "300px", size = 'sm',
              tooltip = tooltipOptions(title = "Adjust Plot Size")
              ),
            tags$div(
              style = "overflow-x: auto; overflow-y: auto; width: 100%;", 
              uiOutput("func_bar_plot_ui")
              ),
            br(),
            shinyjs::hidden(div(
              id = "func_download2_ui", downloadButton("func_download2", "Download Plot")
            ))
          )
        )
      )
    )
  )
)











