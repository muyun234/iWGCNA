# server

bullet <- function(..., icon_name = "diamond", icon_color = "#2c7be5") {
  tags$div(
    style = "
      display: flex;
      align-items: flex-start;
      font-size: 18px;
      line-height: 1.6;
      margin-bottom: 6px;
    ",
    
    icon(icon_name, style = paste0(
      "margin-right: 8px; margin-top: 3px; color: ", icon_color, ";"
    )),
    
    tags$div(...)
  )
}

output$background <- renderUI({
  box(
    width = 12, elevation = 4, collapsible = TRUE, solidHeader = TRUE, status = "primary",
    title = span(style = "font-size: 20px; font-weight: bold;", "Background (Why iWGCNA?)"),
    icon = icon("diamond"),
    bullet(
      "WGCNA is a widely used method for identifying gene co-expression networks 
      and functional gene groups."
    ),
    bullet(
      "With the decreasing cost of RNA sequencing, researchers increasingly", 
      tags$strong("analyze large cohorts"), " and multiple datasets."
    ),
    bullet(
      "However, traditional WGCNA workflows are ", tags$strong("computationally intensive"), 
      ", fragmented, and difficult to interpret, often requiring long runtimes and 
      multiple external tools."
    ),
    bullet(
      "As a result, no dedicated ", tags$strong("end-to-end web-based WGCNA platform"), 
      " has been available."
    ),
    bullet(
      "iWGCNA bridges this gap by enabling ", 
      tags$strong("fast, interactive, and interpretable co-expression analysis"), 
      " directly through a web browser."
      )
  )
})

output$user <- renderUI({
  box(
    width = 12, elevation = 4, collapsible = TRUE, solidHeader = TRUE, status = "success",
    title = span(style = "font-size: 20px; font-weight: bold;", "Who Should Use iWGCNA?"),
    icon = icon("users"),
    bullet(
      "Biologists performing RNA-seq data analysis", 
      icon_name = "users", icon_color = "#28a745"
    ),
    bullet(
      "Researchers studying gene co-expression and pathways", 
      icon_name = "users", icon_color = "#28a745"
    ),
    bullet(
      "Users who want WGCNA without complex scripting", 
      icon_name = "users", icon_color = "#28a745"
    ),
    bullet(
      "Labs analyzing multiple datasets or large cohorts", 
      icon_name = "users", icon_color = "#28a745"
    ),
    bullet(
      "Anyone seeking interpretable network biology results.", 
      icon_name = "users", icon_color = "#28a745"
    )
  )
})

output$provide <- renderUI({
  box(
    width = 12, elevation = 4, collapsible = TRUE, solidHeader = TRUE, status = "orange",
    title = span(style = "font-size: 20px; font-weight: bold;", "What iWGCNA Provides?"),
    icon = icon("rocket"),
    div(
      "iWGCNA is a one-stop WGCNA web server that integrates: ",
      style = "font-size: 20px; line-height: 1.6; margin-bottom: 6px;"
      ),
    bullet(
      "Fast co-expression network construction", 
      icon_name = "rocket", icon_color = "#ff851b"
    ),
    bullet(
      "Functional enrichment analysis", 
      icon_name = "rocket", icon_color = "#ff851b"
    ),
    bullet(
      "Interactive network visualization", 
      icon_name = "rocket", icon_color = "#ff851b"
    ),
    bullet(
      "Integrative heatmaps combining expression and functional annotation", 
      icon_name = "rocket", icon_color = "#ff851b"
    ),
    bullet(
      "All within a single, unified, browser-based platform", 
      icon_name = "rocket", icon_color = "#ff851b"
    )
  )
})

output$highlight <- renderUI({
  box(
    width = 12, elevation = 4, collapsible = TRUE, solidHeader = TRUE, status = "warning",
    title = span(style = "font-size: 20px; font-weight: bold;", "Performance Highlight"),
    icon = icon("bolt"),
    div(
      strong("Optimized for large-scale transcriptomic data"),
      style = "font-size: 20px; line-height: 1.6; margin-bottom: 6px;"
    ),
    img(src = "speed.png", style = "display: block; margin: 0 auto;", width = "80%"),
    div(
      "Enables ", strong("real-time, interactive WGCNA analysis on the web."),
      style = "font-size: 18px; line-height: 1.6; margin-bottom: 6px;"
    )
  )
})


output$access <- renderUI({
  box(
    width = 12, elevation = 4, collapsible = TRUE, solidHeader = TRUE, status = "indigo",
    title = span(style = "font-size: 20px; font-weight: bold;", "Access"),
    icon = icon("globe"),
    bullet(
      "Web-based (no installation required)", 
      icon_name = "globe", icon_color = "#6610f2"
    ),
    bullet(
      "Free and open access", 
      icon_name = "globe", icon_color = "#6610f2"
    ),
    bullet(
      "Runs on all modern web browsers", 
      icon_name = "globe", icon_color = "#6610f2"
    ),
    div(style = "margin-bottom: 5px;"),
    tags$a(
      href = "https://shiny.cgu.edu.tw/app/iWGCNA",
      target = "_blank",
      icon(
        "hand-point-right", style = 
        "margin-right: 8px; margin-top: 3px; color: #6610f2;"
      ),
      "https://shiny.cgu.edu.tw/app/iWGCNA"
    )
  )
})



