# UI

bs4Dash::tabItem(
  tabName = "upload",
  tags$head(tags$style(HTML(paste(
    readLines("www/analysis_style.txt", warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )))),
  div(
    id = "iw-upload", class = "analysis-page",
    div(
      class = "analysis-page-heading",
      tags$h2("1. Data Upload"),
      tags$p("Choose your data format and organism, then upload an expression matrix or try the demo.")
    ),
    bs4Dash::box(
      width = 12, elevation = 1, title = "Data Upload",
      collapsible = FALSE, solidHeader = TRUE, status = "info",
      fluidRow(
        column(
          width = 8,
          div(
            class = "analysis-settings",
            div(
              class = "analysis-step",
              div(
                class = "analysis-step-heading",
                span(class = "analysis-step-number", "1"),
                div(
                  tags$h3("Expression data format"),
                  tags$p("Choose the type of expression values in your input file.")
                )
              ),
              div(
                class = "analysis-choice-cards",
                shinyWidgets::prettyRadioButtons(
                  inputId = "raw_format", label = NULL,
                  choiceNames = list(
                    tags$span(
                      tags$strong("Read count"),
                      tags$small("Raw read counts before normalization.")
                    ),
                    tags$span(
                      tags$strong("Normalized read count"),
                      tags$small("Normalized expression values, such as FPKM, RPKM, TPM or CPM.")
                    )
                  ),
                  choiceValues = c("count", "normalized_count"),
                  selected = "count", inline = FALSE,
                  status = "info", shape = "round", fill = TRUE, width = "100%"
                )
              ),
              tags$p(
                class = "analysis-hint",
                "More about normalized expression values ",
                actionButton("rf", label = "", icon = icon("info"),
                             class = "btn-info btn-xs", title = "Expression format information",
                             `aria-label` = "Expression format information")
              )
            ),
            div(
              class = "analysis-step",
              div(
                class = "analysis-step-heading",
                span(class = "analysis-step-number", "3"),
                div(
                  tags$h3("Upload expression data"),
                  tags$p("Select a file or drag it into the upload area.")
                )
              ),
              div(id = "upload_ui", uiOutput("exp_file_ui")),
              tags$p(class = "analysis-hint", "Accepted formats: .txt, .tsv and .csv. Maximum file size: 30 MB.")
            ),
            shinyjs::hidden(div(
              id = "run_ui",
              bs4Dash::actionButton(
                "data_upload", "Load Expression Data", status = "info",
                class = "analysis-button analysis-run", icon = icon("upload")
              )
            ))
          )
        ),
        column(
          width = 4,
          div(
            class = "analysis-step",
            div(
              class = "analysis-step-heading",
              span(class = "analysis-step-number", "2"),
              div(
                tags$h3("Organism"),
                tags$p("Select the organism represented in your expression matrix.")
              )
            ),
            uiOutput("select_organism_ui"),
            tags$p(
              class = "analysis-hint",
              "Both Gene Symbol and Ensembl Gene ID formats are supported.",
              tags$br(),
              "Note: For Ensembl Gene IDs, please remove version suffixes 
                (e.g., ENSG00000141510.18 → ENSG00000141510) before uploading to facilitate PPI analysis."
            )            
          ),
          div(
            class = "analysis-step",
            div(
              class = "analysis-step-heading",
              span(class = "analysis-step-number", icon("flask"), `aria-hidden` = "true"),
              div(
                tags$h3("Try demo data"),
                tags$p("Explore iWGCNA with an example dataset.")
              )
            ),
            shinyWidgets::awesomeCheckbox(
              inputId = "demo_on", label = "Use Demo", value = FALSE, status = "info"
            ),
            tags$details(
              tags$summary(style = "cursor: pointer;", "Show description"),
              tags$p(
                class = "analysis-hint",
                "The demonstration dataset originates from GSE15459 and contains normalized gene 
                expression profiles from 192 gastric cancer patient samples."
              ),
              tags$p(class = "analysis-hint", "Reference: J Cancer. 2024;15(11):3284–3296")
            ),
            div(
              class = "analysis-next",
              downloadButton("download", "Download Demo", class = "btn-outline-info analysis-button")
            )
          )
        )
      ),
      fluidRow(
        column(
          width = 12,
          shinyjs::hidden(div(
            id = "gene_type_ui",
            hr(),
            shinyWidgets::prettyRadioButtons(
              inputId = "gene_type", label = "Gene ID Conversion:",
              choices = c("Symbol", "Ensembl ID"), selected = character(0),
              inline = TRUE, status = "info", fill = TRUE
            )
          )),
          shinyjs::hidden(div(
            id = "step2_ui", class = "analysis-next",
            bs4Dash::actionButton(
              "step2", "STEP 2. Preprocessing", status = "info",
              icon = icon("arrow-right"), class = "analysis-button analysis-next-button"
            )
          ))
        )
      ),
      fluidRow(
        column(
          width = 12,
          hr(),
          div(
            class = "analysis-results",
            div(
              id = "information_ui",
              tags$h3(tags$strong("Input file requirements")),
              tags$p(
                class = "analysis-hint",
                "The input file must be a gene expression matrix in .txt, .tsv or .csv format. The first column 
                contains gene names, and all subsequent columns represent samples. Matrix entries should be 
                numeric expression values. Both raw read counts and normalized values (e.g., FPKM, RPKM, TPM 
                and CPM) are supported."
              ),
              img(src = "data_input.png", width = "50%",
                  style = "display: block; max-width: 100%; height: auto; margin-top: 14px;",
                  alt = "Example expression matrix with gene names in the first column and 
                  samples in subsequent columns")
            ),
            shinyjs::hidden(div(
              id = "preview_ui",
              tags$h3("Data Preview"),
              fluidRow(
                column(width = 4, uiOutput("gene_count_ui")),
                column(width = 4, uiOutput("sample_count_ui"))
              ),
              hr(),
              DT::dataTableOutput("raw_data")
            ))
          )
        )
      )
    )
  )
)
