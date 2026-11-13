# Server
shiny::hideTab(inputId = "create_tabs", target = "module_detection_tab")

##################################################################
# step 1 : Calculate Soft threshold
networktype = eventReactive(input$cre_run,{
  as.character(input$Networktype)
})

rscut = eventReactive(input$cre_run ,{
  as.numeric(input$cutoff)
})

mbs = eventReactive(input$cre_run,{
  as.numeric(input$reserved_gene_num)
})

observeEvent(input$cre_run, {
  
  if (is.null(exp.ds$data) || anyNA(exp.ds$data)){
    shinyalert("Please upload a file", type = "warning")
    updateTabItems(session, "steps", "upload")
    return()
    }
  
  if(is.null(exp.ds$table2)){
    shinyalert("Please preprocess the data", type = "warning")
    updateTabItems(session, "steps", "preprocessing")
    return()
    }
  
  exp.ds$sft = list()
  
  showModal(modalDialog(
    title = "Loading...",
    easyClose = FALSE,
    footer = NULL,
    div(id = "creat_step1", tags$p("Calculate the power......")),
    div(id = "creat_step2", shinyjs::hidden(tags$p("Ploting......")))
  ))
  
  print(networktype())
  power <- c(c(1:10), seq(from = 12, to = 30, by = 2))
  
  exp.ds$sft <- WGCNA::pickSoftThreshold(
    exp.ds$table2, powerVector = power, networkType = networktype(), 
    RsquaredCut = rscut(), verbose = 5, blockSize = mbs())
  
  sft.data <- exp.ds$sft$fitIndices
  sft.data$signed_R2 <- -sign(sft.data$slope) * sft.data$SFT.R.sq
  exp.ds$sft.data <- sft.data
  
  shinyjs::show("creat_step2")
  ## infobox4 for Recommend power
  eligible <- is.finite(sft.data$signed_R2) & sft.data$slope < 0 & sft.data$signed_R2 > rscut()
  
  exp.ds$recommended_power <- if (any(eligible)) {
    min(sft.data$Power[eligible])
  } else {
    NA_real_
  }
  
  output$recommend <- renderUI({
    power <- exp.ds$recommended_power
    summaryBox::summaryBox(
      "Recommended power",
      if (is.na(power)) "No power meets the cutoff" else as.character(power),
      width = 12, icon = "fas fa-clipboard-list", style = if (is.na(power)) "warning" else "success"
    )
  })
  
  output$cre_choose_ui <- renderUI({
    radioGroupButtons(
      inputId = "cre_choose",
      choices = c(`<i class='fa fa-pie-chart'></i>` = "plot",
                  `<i class='fa fa-table'></i>` = "table"),
      selected = "plot", justified = TRUE, width = "50", status = "info"
    )
  })
  
  ## show sft plot
  output$cre_sftplot <- renderPlot({
    if(is.null(exp.ds$table2)){return()}
    if(length(exp.ds$sft) == 0){return()}
    print("start plotting")
    
    p1 <- ggplot(exp.ds$sft.data, aes(Power, signed_R2, label = Power)) +
      geom_point(color = "gray10") +
      ggrepel::geom_text_repel(aes(label = Power)) +
      geom_hline(yintercept = rscut(), color = "red") +
      labs(
        x = "Soft threshold (power)",
        y = expression(paste("Scale-free topology fit - ", R^{2})),
        title = "Scale independence"
      ) + 
      theme_bw()
    
    p2 <- ggplot(exp.ds$sft.data, aes(Power, mean.k., label = Power)) +
      geom_point(color = "gray10") +
      ggrepel::geom_text_repel(aes(label = Power)) +
      labs(
        x = "Soft threshold (power)",
        y = "Mean connectivity (k)",
        title = "Mean connectivity"
      ) + 
      theme_bw()
    
    patchwork::wrap_plots(p1, p2)
    
  }, width = reactive(input$sftplot_width), height = reactive(input$sftplot_height))
  
  ## Show sft table 
  output$cre_sfttab = renderReactable({
    if(is.null(exp.ds$table2)){return()}
    if(length(exp.ds$sft) == 0){return()}
    
    df <- exp.ds$sft.data
    power_estimate <- exp.ds$recommended_power
    
    # 確保df的所有row都是num
    numeric_columns <- sapply(df, is.numeric)
    df[numeric_columns] <- lapply(df[numeric_columns], round, digits = 2)    
    
    reactable(
      df, bordered = TRUE, striped = TRUE, paginateSubRows = TRUE, showPageSizeOptions = TRUE, 
      pageSizeOptions = c(5, 10, 20), defaultPageSize = 10, searchable = TRUE, 
      rowClass = function(index) {
        
        if (!is.na(power_estimate) && df[index, "Power"] == power_estimate) {
          "highlight"
        } else if (!is.na(power_estimate) && df[index, "Power"] >= power_estimate) {
          "bold"
        }
      },
      columns = list(
        Power = colDef(
          style = function(value) {
            if (!is.na(power_estimate) && value >= power_estimate) {
              color <- "#17a2b8"
            } else {
              color <- "#000"
            }
            list(color = color)
            }
          )
        )
      )
     })  
  
  output$powerselect_ui <- renderUI({
    sliderInput(
      inputId = "powerselect",
      label = tags$span(style = "font-size: 18px;", "Soft-Threshold Power (β): "), 
      min = 1, max = 33, value = paste0(exp.ds$recommended_power)
    )
  })
  
  shinyjs::show("cre_download")
  shinyjs::show("cre_download_2")
  shinyjs::show("step3b_ui")
  shiny::showTab(inputId = "create_tabs", target = "module_detection_tab")
  removeModal()
  })

## call next step btn
observeEvent(input$step3b,{
  updateTabsetPanel(session, "create_tabs", "module_detection_tab")
  print("Step 3-2 button clicked, proceed to 3-2_Module detection.")
})

##################################################################
## step 2 : Module identification
power = eventReactive(input$cre_run2,{
  as.numeric(input$powerselect)
})
mms = eventReactive(input$cre_run2,{
  as.numeric(input$minMsize)
})
mch = eventReactive(input$cre_run2,{
  as.numeric(input$mch)
})
tomtype = eventReactive(input$cre_run2,{
  # as.character(input$TOMtype)
  if(networktype() == "unsigned") {
    return("unsigned")  
  } else if (networktype() == "signed" || networktype() == "signed hybrid") {
    return ("signed")
  }
})

observe({
  req(input$reserved_gene_num)
  Max <- as.numeric(input$reserved_gene_num)
  Min <- if (Max < 100) Max else 100
  Value <- if (Max < 600) Max else 600
  updateSliderInput(session, "cre_tom_num", max = Max, min = Min, value = Value)
})

observeEvent(input$cre_run2, {
  if(is.null(exp.ds$table2)){
    shinyalert("Please preprocess the data", type = "warning")
    updateTabItems(session, "steps", "preprocessing")
    return()
  }
  
  if(is.null(power())){
    shinyalert("Please complete the calculate Power step before proceeding", type = "warning")
    updateTabItems(session, "steps", "create")
    return()
  }
  
  print("Creat_module_btn is clicked!!")

  showModal(modalDialog(
    title = "Loading...",
    easyClose = FALSE,
    footer = NULL,
    div(id = "creat_step1", tags$p("Building the network......")),
    div(id = "creat_step2", shinyjs::hidden(tags$p("Analyzing results & plotting......")))
  ))
  
  exp.ds$plot_tomhm <- FALSE
  
  # 確保所有 datExpr 的列都是數值類型
  datExpr <- as.data.frame(lapply(exp.ds$table2, as.numeric))
  
  print(power())
  print(mms())
  print(mch())
  print(networktype())
  print(tomtype())
  print(mbs())
  
  tom_path <- file.path(data_folder$folder_path, "TOM")
  # mbs <- nrow(table2)
  cor <- WGCNA::cor
  exp.ds$net <- WGCNA::blockwiseModules(
    exp.ds$table2, maxBlockSize = mbs(), networkType = networktype(), TOMType = tomtype(), 
    minModuleSize = mms(), power = power(), mergeCutHeight = mch(), 
    numericLabels = TRUE, verbose = 3, saveTOMs = TRUE, saveTOMFileBase = tom_path
    )
  cor <- stats::cor
  
  shinyjs::show("creat_step2")
  
  exp.ds$nSamples = nrow(exp.ds$table2)
  print(table(exp.ds$net$colors))
  
  exp.ds$moduleLabels = exp.ds$net$colors
  exp.ds$moduleColors <- WGCNA::labels2colors(exp.ds$net$colors)
  exp.ds$moduleUnColors <- WGCNA::labels2colors(exp.ds$net$unmergedColors)
  
  exp.ds$MEs = exp.ds$net$MEs
  MEs_num <- as.numeric(sub("ME", "", colnames(exp.ds$MEs)))
  num_to_colcor <- paste0("ME", WGCNA::labels2colors(MEs_num))
  exp.ds$MEs_col <- exp.ds$MEs
  colnames(exp.ds$MEs_col) <- num_to_colcor

  exp.ds$Gene2module <- data.frame(
    GID = names(exp.ds$net$colors),
    Module = exp.ds$moduleColors
  )
  rownames(exp.ds$Gene2module) <- NULL
  # print("Gene2module: ")
  # print(head(exp.ds$Gene2module$GID[exp.ds$Gene2module$Module %in% "red"]))
  
  # Save RDS files here
  module_numbers <- as.data.frame(table(exp.ds$moduleColors))
  colnames(module_numbers) <- c("Gene_Group_name","Gene_count")  
  print(head(module_numbers))
  
  module_numbers$Gene_count <- sapply(module_numbers$Gene_Group_name, function(name) {
    paste(exp.ds$Gene2module[exp.ds$Gene2module$Module == name, 1], collapse = ", ")
  })
  exp.ds$module_numbers <- module_numbers
  
  file2_path <- file.path(data_folder$folder_path, "gene_mapping.rds")
  saveRDS(module_numbers, file2_path)
  
  # kWithin
  df <- WGCNA::intramodularConnectivity.fromExpr(datExpr = exp.ds$table2, colors = exp.ds$moduleColors, 
                                                 power = power(), networkType = networktype())
  connectivity_df <- data.frame(
    Gene = colnames(exp.ds$table2),
    Module = exp.ds$moduleColors,
    df, row.names = NULL
  )
  
  exp.ds$connectivity_df <- connectivity_df
  
  # exp.ds$TOM <- TOMsimilarityFromExpr(exp.ds$table2, power = power(), TOMType = tomtype(), networkType = networktype())
  e <- new.env()
  load(exp.ds$net$TOMFiles[1], envir = e)
  exp.ds$TOM <- as.matrix(e$TOM)
  probes = names(exp.ds$table2)
  rownames(exp.ds$TOM) <- probes
  colnames(exp.ds$TOM) <- probes
  
  exp.ds$geneTree = exp.ds$net$dendrograms[[1]]
  
  output$cre_choose2_ui <- renderUI({
    radioGroupButtons(
      inputId = "cre_choose2",
      choices = c(`<i class='fa fa-pie-chart'></i>` = "plot",
                  `<i class='fa fa-table'></i>` = "table"),
      selected = "plot", justified = TRUE, width = "50", status = "info"
    )
  })
  
  ## show cluster plot
  output$cre_clustPlot <- renderPlot({
    if(is.null(exp.ds$net)){return()}
    
    cluster_color <- exp.ds$moduleColors
    names(cluster_color) <- names(exp.ds$net$colors)
    
    cluster_uncolor <- exp.ds$moduleUnColors
    names(cluster_uncolor) <- names(exp.ds$net$unmergedColors)
    
    WGCNA::plotDendroAndColors(exp.ds$net$dendrograms[[1]], cbind(cluster_uncolor, cluster_color),
                               c("unmerged", "merged"), dendroLabels = FALSE, addGuide = TRUE, hang = 0.03,
                               guideHang = 0.05)
  }, width = reactive(input$creplot_width), height = reactive(input$creplot_height))
  
  # plot TOM heatmap
  exp.ds$plot_tomhm <- TRUE
  
  ## show module to gene
  tbl <- table(exp.ds$Gene2module$Module)
  exp.ds$module_tb <- data.frame(Module = paste0(names(tbl), " (", as.vector(tbl), " genes", ")"), module = paste0(names(tbl)))
  
  genes_per_module <- tapply(exp.ds$Gene2module$GID, exp.ds$Gene2module$Module, paste, collapse = ",")
  exp.ds$module_tb$Genes <- genes_per_module[exp.ds$module_tb$module]
  
  output$module_gene <- renderReactable({
    csvDownloadButton <- function(tableId, label = "Download CSV", filename = "data.csv") {
      htmltools::tags$button(
        label,
        onclick = sprintf("Reactable.downloadDataCSV('%s', '%s')", tableId, filename)
      )
    }
    hex_to_rgba <- function(hex, alpha = 0.5) {
      rgb <- grDevices::col2rgb(hex) / 255
      sprintf("rgba(%d,%d,%d,%.2f)", rgb[1]*255, rgb[2]*255, rgb[3]*255, alpha)
    }
    
    reactable(
      exp.ds$module_tb, columns = list(
        Module = colDef(name = "Gene Group (Number of Genes)", minWidth = 50),
        module = colDef(show = FALSE),
        Genes =  colDef(name = "Gene IDs", minWidth = 100)
      ), 
      searchable = TRUE, wrap = FALSE,
      details = function(index) {
        plant_data <- exp.ds$Gene2module[exp.ds$Gene2module$Module == exp.ds$module_tb$module[index], ]
        plant_data <- plant_data %>%
          dplyr::rename(Gene = GID, Module = Module) %>%
          dplyr::select(Module, Gene)
        
        tableId <- paste0("plant_data_", index)
        fileName <- paste0("Gene_Group_", exp.ds$module_tb$module[index], ".csv")
        
        htmltools::div(
          style = "padding: 1rem;",
          htmltools::browsable(
            tagList(
              reactable(
                plant_data, outlined = TRUE, searchable = TRUE, defaultPageSize = 10,
                showPageSizeOptions = TRUE, pageSizeOptions = c(5, 10, 20, 50),
                theme = reactableTheme(
                  headerStyle = list(backgroundColor = hex_to_rgba(exp.ds$module_tb$module[index], 0.5))
                ),
                elementId = tableId
              ),
              csvDownloadButton(tableId, "Download CSV", filename = fileName)
            )
          )
        )
      }
    )
  })
  
  # module bar plot
  module_size <- exp.ds$Gene2module %>%
    dplyr::filter(Module != "grey") %>%
    count(Module, name = "Gene_count") %>%
    arrange(desc(Gene_count))
  
  output$cre_module_barplot <- renderPlot({
    ggplot(module_size, aes(x = Module, y = Gene_count, fill = Module)) +
      geom_col(width = 0.8) +
      geom_text(aes(label = Gene_count), vjust = -0.4, size = 3.5) +
      scale_fill_manual(
        values = setNames(
          as.character(module_size$Module),
          as.character(module_size$Module)
        )
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      labs(x = "Module", y = "Number of genes") +
      theme_classic(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
  }, width = input$cre_module_barplot_width, input$cre_module_barplot_height)
  
  exp.ds$KME <- NULL
  shinyjs::show("cre_download_4")
  shinyjs::show("cre_download_5")
  shinyjs::show("cre_download_7")
  shinyjs::show("cre_download_8")
  shinyjs::show("step4_ui")
  shinyjs::show("step6_ui")
  removeModal()
})

output$cre_TOMPlot <- renderPlot({
  req(exp.ds$plot_tomhm)
  
  set.seed(123)
  gene = names(exp.ds$table2)
  partgene <- sample(gene, size = as.integer(input$cre_tom_num))
  
  dissTOM <- 1 - exp.ds$TOM
  partTOM <- dissTOM[partgene, partgene]
  
  selectTree = hclust(as.dist(partTOM), method = "average")
  selectColors <- exp.ds$Gene2module$Module[match(partgene, exp.ds$Gene2module$GID)]
  
  plotTOM = partTOM ^ 7
  diag(plotTOM) = NA
  
  colors <- grDevices::heat.colors(50)
  
  name <- paste0("Topological Overlap Matrix heatmap (",as.integer(input$cre_tom_num) , " genes)")
  
  WGCNA::TOMplot(plotTOM, selectTree, selectColors, main = name, col = colors)
}, width = reactive(input$cre_tom_width), height = reactive(input$cre_tom_height))


## call next step btn
observeEvent(input$step4,{
  updateTabItems(session, "steps", "interested")
  print("Step 4 button clicked, proceed to 4_Gene Group-Trait Association.")
})
observeEvent(input$step6,{
  updateTabItems(session, "steps", "function")
  print("Step 6 button clicked, proceed to 6_Functional Enrichment Analysis.")
})

##################################################################
## download

output$credownload <- downloadHandler(
  filename = function() {
    paste("03_SftPlot.png")
  },
  content = function(file) {
    if(is.null(exp.ds$table2)){return()}
    if(length(exp.ds$sft) == 0){return()}
    
    plot_width <- input$sftplot_width / 96 # 將像素轉換為英寸（DPI 假設為 96）
    plot_height <- input$sftplot_height / 96 # 將像素轉換為英寸（DPI 假設為 96）
    
    sft.data <- exp.ds$sft$fitIndices
    p1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
      geom_point(color = "gray10") +
      ggrepel::geom_text_repel(aes(label = Power)) +
      geom_hline(yintercept = 0.8, color = "red") +
      labs(
        x = "Soft threshold (power)",
        y = expression(paste("Scale-free topology fit - ", R^{2})),
        title = "Scale independence"
      ) + 
      theme_bw()
    
    p2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
      geom_point(color = "gray10") +
      ggrepel::geom_text_repel(aes(label = Power)) +
      labs(
        x = "Soft threshold (power)",
        y = "Mean connectivity (k)",
        title = "Mean connectivity"
      ) + 
      theme_bw()
    
    sft_plot <- p1 + p2
    ggsave(filename = file, plot = sft_plot, 
           width = plot_width, height = plot_height, dpi = 300)
  }
)

output$credownload2 <- downloadHandler(
  filename = function() {
    paste("03_SftTable.csv")
  },
  content = function(file) {
    df <- round(exp.ds$sft$fitIndices, digits = 2)
    write.table(x = df, file = file, sep = ",", quote = FALSE)
  }
)

output$credownload4 <- downloadHandler(
  filename = function() {
    paste("03_ClusterDendrogram.png")
  },
  content = function(file) {
    if (is.null(exp.ds$net)) return()
    
    # 因為Layout，所以直接設定好長寬才不會跑掉
    png(filename = file, width = 3000, height = 2000, res = 300)
    par(mar = c(2, 2, 1, 1)) # 增加上下左右的邊距
    # par(mar = c(bottom, left, top, right))
    
    cluster_color <- exp.ds$moduleColors
    names(cluster_color) <- names(exp.ds$net$colors)
    
    cluster_uncolor <- exp.ds$moduleUnColors
    names(cluster_uncolor) <- names(exp.ds$net$unmergedColors)
    
    WGCNA::plotDendroAndColors(
      exp.ds$net$dendrograms[[1]], cbind(cluster_uncolor, cluster_color), c("unmerged", "merged"), 
      dendroLabels = FALSE, addGuide = TRUE, hang = 0.03, guideHang = 0.05
    )
    dev.off()
  }
)

output$credownload5 <- renderUI({
  screenshotButton(
    label = "Download Plot",
    download = TRUE, 
    server_dir = tempdir(), 
    id = "cre_TOMPlot", 
    filename = "03_TOMPlot")
})

output$credownload7 <- downloadHandler(
  filename = function() {
    paste("03_Gene2GeneGroup.csv")
  },
  content = function(file) {
    write.table(x = exp.ds$module_numbers, file = file, sep = ",", quote = FALSE, row.names = FALSE)
  }
) 

output$credownload8 <- renderUI({
  screenshotButton(
    label = "Download Plot",
    download = TRUE, 
    server_dir = tempdir(), 
    id = "cre_module_barplot", 
    filename = "03_Module_bar_plot"
  )
})

##################################################################
## show icon info.
# choose R2 cutoff
addPopover(
  id = "Cutoff",
  options = list(
    content = paste0("Set the cut-off value for the  scale free topology model fit (default : 0.8)",
                     br(),
                     "Range : 0.0-1.0"),
    title = "More information",
    placement = "right",
    trigger = "click",
    html = TRUE  # Ensure HTML content is parsed
  )
)

# select power
# addPopover(
#   id = "ps1",
#   options = list(
#     content = paste0("Choose an appropriate or recommended power."),
#     title = "More information",
#     placement = "right",
#     trigger = "click",
#     html = TRUE  # Ensure HTML content is parsed
#   )
# )

# select minModuleSize
addPopover(
  id = "mms2",
  options = list(
    content = paste0("The gene group minimum size for gene group detection is set to 30 by default",
                     br(),
                     "If there are not enough genes, they cannot be divided into a gene group",
                     br(),
                     a("",href = "https://rdrr.io/cran/WGCNA/man/blockwiseModules.html", target="_blank")
    ),
    title = "More information",
    placement = "right",
    trigger = "click",
    html = TRUE  # Ensure HTML content is parsed
  )
)
# select mergeCutHeight  
addPopover(
  id = "mch2",
  options = list(
    content = paste0("Dendrogram cut height for gene group merging.",br(),
                     "By default, the value is 0.25."),
    title = "More information",
    placement = "right",
    trigger = "click",
    html = TRUE  # Ensure HTML content is parsed
  )
)

# Network type
addPopover(
  id = "network",
  options = list(
    content = paste0(
      "<b>Unsigned</b>: uses the absolute value of correlations (|cor|^β). Both positive and negative correlations are treated equally.",br(),
      "<b>Signed</b>: keeps the sign of correlation ((1 + cor)/2)^β. Only positive correlations are emphasized.", br(),
      "<b>Signed Hybrid</b>: negative correlations are set to zero; positive ones use ((1 + cor)/2)^β."
      ),
    title = "More information",
    placement = "right",
    trigger = "click",
    html = TRUE
  )
)

# TOM type
addPopover(
  id = "tom",
  options = list(
    content = paste0(
      "<b>None</b>: uses adjacency only", br(),
      "<b>Unsigned</b>: ignores the sign of correlation", br(),
      "<b>Signed</b>: considers the correlation sign — positive links reinforce, negative links weaken."
    ),
    title = "More information",
    placement = "right",
    trigger = "click",
    html = TRUE
  )
)










