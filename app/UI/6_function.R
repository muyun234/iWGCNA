# UI: Step 6. Keep the original 3 / 9 column layout.
tabItem(
  tabName = "function",
  div(
    id = "iw-function", class = "analysis-page",
    div(class = "analysis-page-heading",
        tags$h2("6. Functional Enrichment Analysis"),
        tags$p("Explore the biological functions and pathways of selected modules.")),
    fluidRow(
      box(
        id = "functional_analysis", height = "auto", width = 12,
        elevation = 1, title = "Functional enrichment",
        collapsible = TRUE, solidHeader = FALSE, maximizable = TRUE,
        status = "info",
        fluidRow(
          column(
            width = 3,
            div(
              class = "analysis-step",
              div(class = "analysis-step-heading",
                  span(class = "analysis-step-number", "1"),
                  div(tags$h3("Modules"),
                      tags$p("Choose the modules to analyze."))),
              selectInput("func_module_select", "Modules:",
                          choices = NULL, multiple = TRUE, width = "100%"),
              tags$p(class = "analysis-hint", "Select one or more modules, or All to include every module except grey.")
            ),
            div(
              class = "analysis-step",
              div(class = "analysis-step-heading",
                  span(class = "analysis-step-number", "2"),
                  div(tags$h3("Functional enrichment"),
                      tags$p("Choose GO categories and/or KEGG pathways."))),
              div(id = "func_db_type_ui",
                  selectInput("func_string_version", "STRING database version:",
                              choices = c("11.5", "12.0", "12.5"),
                              selected = "12.5", multiple = FALSE, width = "100%")),
              div(
                class = "analysis-enrichment-section",
                tags$h4("Gene Ontology (GO)"),
                pickerInput(
                  inputId = "func_go_category", label = "GO category:",
                  choices = c("None" = "", "Biological Process" = "BP",
                              "Molecular Function" = "MF", "Cellular Component" = "CC",
                              "All (BP, MF, CC)" = "all"),
                  selected = "BP", multiple = FALSE, width = "100%",
                  options = pickerOptions(clearable = TRUE)
                ),
                shinyjs::hidden(div(
                  id = "func_go_setting_ui",
                  numericInput("func_go_padj", "Adjusted P-value cutoff:",
                               value = 0.05, min = 0, max = 1, step = 0.01, width = "100%"),
                  numericInput("func_go_topn", "GO terms per module (per category):",
                               value = 10, min = 1, step = 1, width = "100%")
                ))
              ),
              div(
                id = "func_kegg_ui", class = "analysis-enrichment-section",
                tags$h4("KEGG pathways"),
                prettyCheckbox("func_kegg", "KEGG pathway analysis",
                               value = FALSE, status = "info", shape = "curve", outline = TRUE),
                shinyjs::hidden(div(
                  id = "func_kegg_setting_ui",
                  numericInput("func_kegg_padj", "Adjusted P-value cutoff:",
                               value = 0.05, min = 0, max = 1, step = 0.01, width = "100%"),
                  numericInput("func_kegg_topn", "KEGG pathways per module:",
                               value = 20, min = 1, step = 1, width = "100%")
                ))
              )
            ),
            div(
              class = "analysis-step",
              div(class = "analysis-step-heading",
                  span(class = "analysis-step-number", "3"),
                  div(tags$h3("Gene labels"),
                      tags$p("Optionally label genes on the heatmap."))),
              selectInput("func_genelist", "Genes to label:", choices = NULL,
                          selected = NULL, multiple = TRUE, width = "100%"),
              tags$p(class = "analysis-hint", "This selection only controls heatmap labels. Enrichment uses the selected modules.")
            ),
            bs4Dash::actionButton("func_update", "Run enrichment analysis",
                                  icon = icon("play"), status = "info", class = "analysis-button analysis-run"),
            uiOutput("func_summary")
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
                  sliderInput("func_heatmap_height", "Plot Height:", min = 300, max = 3000, value = 700, step = 100),
                  sliderInput("func_heatmap_anno_size", "Text Size:", min = 1, max = 20, value = 12, step = 1),
                  radioButtons("func_heatmap_box_style", "Text Box Style:", 
                               choices = c("Link" = "link", "Block" = "block"),
                               selected = "block", inline = TRUE
                  ),
                  # radioButtons("func_heatmap_type", "Plot type:", choices = c("Heatmap" = "heatmap", "Heatmap & Line" = "both"), selected = "both", inline = TRUE),
                  div(id = "func_heatmap_show_mode_ui", 
                      radioButtons(
                        "func_heatmap_show_mode", "Annotation display mode:",
                        choices = c("Together" = "together", "Separate" = "separate"),
                        selected = "separate", inline = TRUE
                      )
                  ),
                  # selectInput("func_heatmap_module_color", "Gene groups with black annotations: ", choices = NULL, multiple = TRUE),
                  awesomeCheckbox("func_heatmap_module_color", "Use black annotation text", value = FALSE),
                  circle = F, status = "info", icon = icon("cog"), width = "300px", size = 'sm',
                  tooltip = tooltipOptions(title = "Adjust Plot Size")
                ),
                tags$div(
                  class = "analysis-plot-scroll", 
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
                circle = F, status = "info", icon = icon("cog"), width = "300px", size = 'sm',
                tooltip = tooltipOptions(title = "Adjust Plot Size")
              ),
              tags$div(
                class = "analysis-plot-scroll", 
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
)
