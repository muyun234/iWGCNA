# UI
bs4Dash::tabItem(
  tabName = "home",
  tags$head(tags$style(HTML(paste(
    readLines("www/analysis_style.txt", warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )))),
  div(
    class = "analysis-page iw-home",
    div(
      class = "home-hero",
      fluidRow(
        bs4Dash::box(
          width = 12, title = NULL, status = NULL, solidHeader = FALSE,
          headerBorder = FALSE, collapsible = FALSE, maximizable = FALSE,
          elevation = 0,
          fluidRow(
            column(
              width = 6,
              tags$p(class = "home-eyebrow", "iWGCNA WEB SERVER"),
              tags$h1(class = "home-title", "Integrated gene", tags$br(),
                      "co-expression analysis"),
              tags$p(class = "home-lead",
                     "Build networks, prioritize hub genes, and interpret biological functions in one browser-based workflow."),
              div(
                class = "home-actions",
                shinyWidgets::actionBttn(
                  inputId = "home_start", label = "Start analysis",
                  style = "simple", color = "primary", size = "lg",
                  no_outline = FALSE, class = "home-button home-button-primary"
                ),
                shinyWidgets::actionBttn(
                  inputId = "home_demo", label = "Try gastric cancer demo",
                  style = "bordered", color = "primary", size = "lg",
                  no_outline = FALSE, class = "home-button home-button-secondary"
                )
              ),
              tags$p(class = "home-note", "No programming or local installation required.")
            ),
            column(
              width = 6
#               tags$img(src = "home_illustration.png", class = "home-illustration",
#                        alt = "Co-expression and protein interaction networks linked to functional enrichment analysis.")
            )
          )
        )
      )
    ),
    div(
      class = "home-features",
      fluidRow(
        bs4Dash::box(
          width = 4, title = NULL, status = NULL, solidHeader = FALSE,
          headerBorder = FALSE, collapsible = FALSE, maximizable = FALSE,
          elevation = 0,
          div(class = "home-feature",
              div(class = "home-feature-icon", `aria-hidden` = "true", icon("th")),
              div(tags$h3("Accelerated network construction"),
                  tags$p("OpenBLAS-optimized WGCNA workflow.")))
        ),
        bs4Dash::box(
          width = 4, title = NULL, status = NULL, solidHeader = FALSE,
          headerBorder = FALSE, collapsible = FALSE, maximizable = FALSE,
          elevation = 0,
          div(class = "home-feature",
              div(class = "home-feature-icon", `aria-hidden` = "true", icon("project-diagram")),
              div(tags$h3("Co-expression + PPI"),
                  tags$p("Compare connectivity to prioritize candidate hub genes.")))
        ),
        bs4Dash::box(
          width = 4, title = NULL, status = NULL, solidHeader = FALSE,
          headerBorder = FALSE, collapsible = FALSE, maximizable = FALSE,
          elevation = 0,
          div(class = "home-feature",
              div(class = "home-feature-icon", `aria-hidden` = "true", icon("globe")),
              div(tags$h3("14-species enrichment"),
                  tags$p("Gene mapping, GO and KEGG analysis, and integrated module heatmaps.")))
        )
      )
    ),
    div(
      class = "home-performance",
      fluidRow(
        bs4Dash::box(
          width = 12, title = NULL, status = NULL, solidHeader = FALSE,
          headerBorder = FALSE, collapsible = FALSE, maximizable = FALSE,
          elevation = 0,
          tags$h2("Performance under the tested configurations"),
          tags$p(class = "home-benchmark-context",
                 "30,000 genes \u00b7 100 RNA-seq samples \u00b7 full TOM"),
          fluidRow(
            column(
              width = 9,
              div(class = "home-chart", HTML('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 840 232" role="img" aria-labelledby="home-chart-title home-chart-desc">
            <title id="home-chart-title">Network analysis elapsed time</title>
            <desc id="home-chart-desc">R WGCNA with reference BLAS: 230.80 minutes. PyWGCNA: 34.49 minutes. iWGCNA with OpenBLAS: 6.99 minutes. Lower is faster.</desc>
            <g font-family="Arial, sans-serif" font-size="16" fill="#292d32">
              <text x="151" y="38" text-anchor="end"><tspan x="151">R WGCNA</tspan><tspan x="151" dy="19">(reference BLAS)</tspan></text>
              <text x="151" y="105" text-anchor="end">PyWGCNA</text>
              <text x="151" y="147" text-anchor="end"><tspan x="151" font-weight="700">iWGCNA</tspan><tspan x="151" dy="19">(OpenBLAS)</tspan></text>
              <rect x="170" y="28" width="553.92" height="29" rx="2" fill="#adb5bd"/>
              <rect x="170" y="83" width="82.776" height="29" rx="2" fill="#79d4da"/>
              <rect x="170" y="138" width="16.776" height="29" rx="2" fill="#17a2b8"/>
              <text x="735.92" y="48">230.80 min</text>
              <text x="264.776" y="103">34.49 min</text>
              <text x="198.776" y="158">6.99 min</text>
              <g stroke="#626970" stroke-width="1.5" fill="none">
                <path d="M170 15V182H770"/>
                <path d="M170 182v7 M290 182v7 M410 182v7 M530 182v7 M650 182v7 M770 182v7"/>
              </g>
              <g text-anchor="middle" font-size="14">
                <text x="170" y="207">0</text><text x="290" y="207">50</text>
                <text x="410" y="207">100</text><text x="530" y="207">150</text>
                <text x="650" y="207">200</text><text x="770" y="207">250</text>
                <text x="470" y="229">Elapsed time (min)</text>
              </g>
            </g>
          </svg>'))
            ),
            column(
              width = 3,
              div(
                class = "home-speedups",
                div(class = "home-speedup", tags$strong("33.02\u00d7"),
                    tags$p("speedup vs R WGCNA with reference BLAS")),
                div(class = "home-speedup", tags$strong("4.93\u00d7"),
                    tags$p("vs PyWGCNA"))
              )
            )
          ),
          tags$p(class = "home-benchmark-note",
                 "Timed network-analysis steps; three runs per configuration. Speedups reflect the tested BLAS configurations.")
        )
      )
    )
  )
)