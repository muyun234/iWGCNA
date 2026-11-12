#ui
tabItem(
  tabName = "upload",
  box(width = 12, elevation = 4, title = div(style = "font-size: 20px; font-weight: bold;", "1. Data Upload"),
      collapsible = FALSE, solidHeader = FALSE, status = "info",
      fluidRow(
        column(
          width = 9,
          div(
            prettyRadioButtons(
              inputId = "raw_format",
              label = tags$span(
                style = "font-size: 18px;", "Expression Data Format: ",
                actionButton("rf", label = "", icon = icon("info"), class = "btn-info btn-xs")
              ), 
              choices = c("Read Count" = "count", "Normalized Read Count" = "normalized_count"), selected = "count",
              inline = TRUE, status = "info", fill = TRUE
            )
          ),
          div(
            tagList(
              uiOutput("select_organism_ui"),
              tags$small(
                style = "color: #555; display: block; margin-bottom: 10px;",
                "⚠️ The platform supports both Gene Symbol and Ensembl Gene ID formats"
              )
            )
          ),
          div(
            id = "upload_ui",
            uiOutput("exp_file_ui")
          ),
          shinyjs::hidden(div(
            id = "run_ui",
            bs4Dash::actionButton(
              "data_upload",em("Run"), status = "primary",
              style = "box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.4); border-radius: 8px;"
            )
          )),
        ),
        column(
          width = 3,
          div(
            style = "margin-bottom: 10px;",
            awesomeCheckbox(
              inputId = "demo_on", 
              label = tagList(
                tags$div(style = "font-size: 18px;", "Use Demo"), 
                tags$details(
                  tags$summary(
                    style = "cursor: pointer; color: #555; font-weight: normal;",
                    "Show description"
                  ),
                  tags$small(
                    style = "color: #555;",
                    "The demonstration dataset originates from GSE15459, a study characterizing 
                    gene expression patterns in gastric tumors to identify oncogenic pathway activity 
                    and its impact on patient outcomes.. The dataset includes 
                    normalized gene expression profiles from 192 patient samples.",
                    tags$br(),
                    "Reference: ",
                    tags$br(),
                    "J Cancer. 2024;15(11):3284–3296"
                  )
                )
              ),
              value = FALSE, status = "danger"),
            downloadButton("download", "Download")
            )
          )
        ),
      fluidRow(
        column(
          width = 12,
          shinyjs::hidden(div(
            id = "gene_type_ui",
            hr(),
            prettyRadioButtons(
              inputId = "gene_type",
              label = tagList(
                tags$div(style = "font-size: 18px;", "Gene ID Conversion: ")
              ),
              choices = c("Symbol", "Ensembl ID"), selected = character(0),
              inline = TRUE, status = "info", fill = TRUE
            )
          )),
          shinyjs::hidden(div(
            id = "step2_ui", 
            div(style = "margin-bottom: 10px;"),
            bs4Dash::actionButton("step2",em("STEP 2. Preprocessing"),
                                  style = "box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.4); 
                                            border-radius: 8px;",
                                  status = "info", icon = icon("forward"))
          ))
        )
      ),
      fluidRow(
        column(
          width = 12,
          hr(),
          div(
            id = "information_ui",
            div(
              "The input file must be a gene expression matrix in .txt, .tsv or .csv format, 
              where the first column contains gene names and all subsequent columns represent samples. 
              Matrix entries should be numeric expression values. Both raw read counts and normalized 
              values (e.g., FPKM, RPKM, TPM, CPM) are supported.",
              style = "text-align: left; text-transform: none; padding: 10px; font-size: 18px;"
            ),
            img(src = "data_input.png", align = "left", width = "50%")
          ),
          shinyjs::hidden(div(
            id = "preview_ui",
            div(style = "font-size: 20px; font-weight: bold;", "Data Preview"),
            fluidRow(
              column(
                width = 4,
                uiOutput("gene_count_ui")
              ),
              column(
                width = 4,
                uiOutput("sample_count_ui")
              )
            ),
            hr(),
            DT::dataTableOutput("raw_data")
          ))
        )
      )
    )
  )









