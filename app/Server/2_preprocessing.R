# server

## step 1 : Set Normalization set.
shiny::hideTab(inputId = "pre_tabs", target = "outlier_tab")

output$method_ui <- renderUI({
  req(input$raw_format)
  shinyWidgets::prettyRadioButtons(
    "normalized_method", label = NULL,
    choiceNames = list(
      tags$span(tags$strong(if (input$raw_format == "count") "VST (DESeq2)" else "Use input values"),
                tags$small(if (input$raw_format == "count") "Variance-stabilizing transformation" else "No additional transformation")),
      tags$span(tags$strong("Log2 transformation"), tags$small("Apply log2(x + 1)"))
    ),
    choiceValues = c(if (input$raw_format == "count") "vst" else "raw", "logarithm"),
    selected = if (input$raw_format == "count") "vst" else "raw",
    inline = FALSE, status = "info", shape = "round", fill = TRUE, width = "100%"
  )
})

output$transformation_hint <- renderUI({
  req(input$raw_format)
  tags$p(class = "analysis-hint", if (input$raw_format == "normalized_count")
    "Already log-transformed? Choose Use input values." else
      "For non-negative raw counts. VST uses DESeq2; log2 applies log2(x + 1).")
})

output$cutoff_ui <- renderUI({
  req(input$raw_format)
  numericInput("expression_cutoff", "Expression cutoff", min = 0,
               value = if (input$raw_format == "normalized_count") 0 else 10,
               step = 0.1, width = "100%")
})

# when cutoff_ui is rendered again after a format change.
observe({
  input$raw_format
  input$expression_cutoff
  shinyjs::toggleState("expression_cutoff", condition = isTRUE(input$low_expression_filter))
  shinyjs::toggleState("sample_percentage", condition = isTRUE(input$low_expression_filter))
})
observe({
  shinyjs::toggleState("cutMethod", condition = identical(input$gene_selection, "top"))
})

output$low_filter_hint <- renderUI({
  if (!isTRUE(input$low_expression_filter)) {
    return(tags$p(class = "analysis-hint", "Filtering is off. No genes will be removed by this step."))
  }
  req(input$expression_cutoff, !is.null(input$sample_percentage))
  tags$p(class = "analysis-hint", paste0(
    "Remove genes with expression < ", input$expression_cutoff,
    " in more than ", input$sample_percentage,
    "% of samples. Filtering uses input values before transformation."
  ))
})

pre_input_matrix <- reactive({
  req(exp.ds$data, ncol(exp.ds$data) >= 2L)
  values <- exp.ds$data[, -1, drop = FALSE]
  validate(need(all(vapply(values, is.numeric, logical(1))), "Expression values must be numeric."))
  x <- as.matrix(values)
  validate(need(nrow(x) > 0L && ncol(x) > 0L && all(is.finite(x)),
                "Expression data must contain finite numeric values."))
  x
})

pre_available_genes <- reactive({
  x <- pre_input_matrix()
  if (!isTRUE(input$low_expression_filter)) return(nrow(x))
  req(!is.null(input$expression_cutoff), !is.null(input$sample_percentage))
  sum(rowMeans(x < input$expression_cutoff) <= input$sample_percentage / 100)
})

output$available_gene_count <- renderText({
  if (is.null(exp.ds$data)) return("Upload data to see available genes")
  paste(format(pre_available_genes(), big.mark = ",", scientific = FALSE, trim = TRUE), "genes available")
})

preprocessing_summary_content <- reactive({
  method_names <- c(raw = "None", logarithm = "Log2(x + 1)", vst = "VST (DESeq2)")
  method <- if (is.null(input$normalized_method)) "Awaiting input" else unname(method_names[input$normalized_method])
  filtering <- if (isTRUE(input$low_expression_filter)) paste0(
    "On; expression < ", input$expression_cutoff, " in > ", input$sample_percentage, "% of samples"
  ) else "Off"
  retained <- if (is.null(exp.ds$data)) "Upload expression data first" else {
    available <- pre_available_genes()
    fmt <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)
    if (identical(input$gene_selection, "all")) paste0("All ", fmt(available), " available genes") else {
      n <- input$reserved_gene_num
      if (is.null(n) || !is.finite(n) || n < 1 || n != floor(n)) "Enter a positive integer for N" else
        paste0("Top ", fmt(min(n, available)), " / ", fmt(available), " available genes")
    }
  }
  tagList(
    tags$p(tags$strong("Transformation: "), method),
    tags$p(tags$strong("Low-expression filtering: "), filtering),
    tags$p(tags$strong("Genes to retain: "), retained),
    tags$p(class = "analysis-summary-note",
           "Settings used for this analysis. WGCNA quality checks may remove additional genes or samples.")
  )
})

preprocessing_summary_saved <- reactiveVal(NULL)
output$preprocessing_summary <- renderUI({
  preprocessing_summary_saved()
})


observe({
  req(input$demo_on)
  if(input$demo_on){
    updatePrettyRadioButtons(session, inputId = "normalized_method", selected = "raw")
    updateSwitchInput(session, "low_expression_filter", value = FALSE)
    updatePrettyRadioButtons(session, "gene_selection", selected = "all")
    updateNumericInput(session, inputId = "expression_cutoff", value = 0)
    updateNumericInput(session, inputId = "sample_percentage", value = 90)
    updatePickerInput(session, inputId = "cutMethod", selected = "MAD")
    updateNumericInput(session, inputId = "reserved_gene_num", value = 21653)
    updatePrettyRadioButtons(session, inputId = "Networktype", selected = "unsigned")
    updateSliderInput(session, inputId = "cutoff", value = 0.9)
    updateSliderInput(session, inputId = "minMsize", value = 30)
    updateSliderInput(session, inputId = "mch", value = 0.25)
    updateSliderInput(session, inputId = "hub_top_cut", value = 15.0)
    updatePrettyRadioButtons(session, inputId = "hub_ppi_level", selected = 900)
    updateSliderInput(session, inputId = "hub_ppi_degree_percent", value = 15.0)
  }
})

observeEvent(input$pre_filter, {
  if (is.null(exp.ds$data) || anyNA(exp.ds$data)){
    shinyalert("Please upload a file", type = "warning")
    updateTabItems(session, "steps", "upload")
    return()
  }
  
  # creat folder
  if(input$demo_on){
    current_time <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
    data_folder$folder_path <- file.path(paste0("Data/", "demo_", current_time))
    dir.create(data_folder$folder_path)
    
    # Save the demo file to the SAVE folder
    file.copy("Demo/GSE15459.tsv", file.path(data_folder$folder_path))
    
    # reset the preplot width
    updateSliderInput(session = session, inputId = "preplot_width", value = 2800)
    updateSliderInput(session = session, inputId = "rm_width", value = 2800)
    
  } else {
    current_time <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
    file_name <- sub("\\.(txt|tsv)$", "", input$exp_file$name)
    data_folder$folder_path <- file.path(paste0("Data/", file_name, "_", current_time))
    dir.create(data_folder$folder_path)
    
    file_path <- file.path(data_folder$folder_path, input$exp_file$name)
    file.copy(input$exp_file$datapath, file_path)
  }
  
  showModal(modalDialog(
    title = "Loading...",
    easyClose = FALSE,
    footer = NULL,
    div(id = "pre_step1", tags$p("Loading data......")),
    div(id = "pre_step2", shinyjs::hidden(tags$p("Filtering & clustering data......"))),
    div(id = "pre_step3", shinyjs::hidden(tags$p("Ploting......")))
  ))
  shiny::hideTab(inputId = "pre_tabs", target = "outlier_tab")
  
  # 清空remove outlier的圖
  shinyjs::hide("step2_2_remove_ui")
  shinyjs::hide("pre_download4_ui")
  shinyjs::hide("step3_1_ui")
  shinyjs::hide("predownload2_ui")
  shinyjs::hide("pre_summary_ui")
  output$pre_choose_remove_ui <- renderUI(NULL)
  output$pre_clustPlot2 = renderPlotly({
    NULL
  })
  output$pre_table = DT::renderDataTable({
    NULL
  })
  
  print("Start preprocessing")
  
  method = input$normalized_method
  
  datatype <- input$raw_format
  filter_enabled <- isTRUE(input$low_expression_filter)
  RcCutoff <- if (filter_enabled) as.numeric(input$expression_cutoff) else 0
  samplePerc <- if (filter_enabled) as.numeric(input$sample_percentage) / 100 else 0
  keep_all <- identical(input$gene_selection, "all")
  requested_n <- input$reserved_gene_num
  cutmethod <- input$cutMethod
  if (filter_enabled && (length(RcCutoff) != 1L || !is.finite(RcCutoff) || RcCutoff < 0 ||
                         length(samplePerc) != 1L || !is.finite(samplePerc) || samplePerc < 0 || samplePerc > 1)) {
    removeModal()
    shinyalert("Invalid filtering settings", type = "warning")
    return()
  }
  if (!keep_all && (length(requested_n) != 1L || !is.finite(requested_n) ||
                    requested_n < 1 || requested_n > 80000 || requested_n != floor(requested_n))) {
    removeModal()
    shinyalert("Please enter an integer N between 1 and 80,000", type = "warning")
    return()
  }
  GNC <- if (keep_all) nrow(exp.ds$data) else as.integer(requested_n)
  
  exp.ds$table = data.frame()
  exp.ds$table2 = data.frame()
  exp.ds$gnccheck = list()
  exp.ds$GNC_check = list()
  exp.ds$GNC = list()
  
  # 在這裡加入對參數的檢查
  print(paste("RcCutoff:", RcCutoff))
  print(paste("Sample Percentage:", samplePerc))
  print(paste("Datatype:", datatype))
  print(paste("Method:", method))
  print(paste("Cut Method:", cutmethod))
  
  table_result <- tryCatch(
    {
      x <- exp.ds$data %>%
        setNames(c("id",colnames(.)[-1])) %>%
        mutate(id = as.character(id)) %>%
        tibble::column_to_rownames("id") 
      
      if (!all(vapply(x, is.numeric, logical(1))) || !all(is.finite(as.matrix(x)))) {
        stop("Expression values must be finite numbers.")
      }
      if (filter_enabled) {
        # Maximum LOW-expression fraction, not minimum expressed fraction.
        x <- x[rowMeans(x < RcCutoff) <= samplePerc, , drop = FALSE]
      }
      if (!nrow(x)) stop("No genes were retained under the current filtering settings.")
      if (method %in% c("vst", "logarithm") && any(x < 0)) {
        stop("VST and log2 transformation require non-negative input values.")
      }
      
      if (datatype == "expected_count") {
        x <- x %>%
          mutate_if(is.numeric,ceiling)
      }
      ## Normalization
      if(method == "vst") {
        dx <- x %>%
          as.matrix() %>%
          DESeq2::varianceStabilizingTransformation(., blind = TRUE)
      } else if (method == "logarithm") {
        dx <- log2(x + 1)
      } else if (method == "raw") {
        dx <- x
      } else {
        stop("Unknown normalization method: ", method)
      }
      
      dx
    },
    error = function(e) {
      shinyalert(
        title = "Error",
        text = paste("Data processing failed:", e$message),
        type = "error"
      )
      return(NULL)
    }
  )
  
  if (is.null(table_result)) {
    removeModal() 
    return()
  }
  
  exp.ds$table <- table_result
  
  if (nrow(exp.ds$table) == 0) {
    shinyalert("No genes were retained under the current conditions", type = "warning")
    removeModal()
  } else {
    shinyjs::show("pre_step2")
    exp.ds$GNC_check = GNC - nrow(exp.ds$table)
    if(exp.ds$GNC_check > 0 ) {
      exp.ds$GNC = nrow(exp.ds$table)
      exp.ds$gnccheck = "The number of Genes you want to retain is greater than the total number of genes after the first filter. The number of genes retained here is equal to the total number after the first filter."
    } else {
      exp.ds$GNC = GNC
      exp.ds$gnccheck = "All good"
    }
    
    print(exp.ds$gnccheck)
    
    if (keep_all || exp.ds$GNC == nrow(exp.ds$table)) {
      datExprVar <- exp.ds$table
    } else {
      score <- if (cutmethod == "MAD") {
        apply(exp.ds$table, 1, mad)
      } else {
        apply(exp.ds$table, 1, var)
      }
      
      # Stable row order breaks ties; retain exactly N before WGCNA QC.
      selected <- head(order(score, decreasing = TRUE, na.last = TRUE), exp.ds$GNC)
      datExprVar <- exp.ds$table[selected, , drop = FALSE]
    }
    
    dim(datExprVar)
    datExpr <- as.data.frame(t(datExprVar))
    gsg = WGCNA::goodSamplesGenes(datExpr, verbose = 3)
    if (!gsg$allOK){
      if (sum(!gsg$goodGenes) > 0)
        printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ",")))
      if (sum(!gsg$goodSamples) > 0)
        printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ",")))
      datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
    }
    exp.ds$table2 <- datExpr
    
    exp.ds$sample_tree <- hclust(dist(exp.ds$table2), method = "average")
    
    shinyjs::show("pre_step3")
    # infobox1 for Expression cutoff
    output$infobox1 <- renderInfoBox({
      infoBox(
        "Total Genes", paste0(nrow(exp.ds$data)), icon = icon("list"),
        color = "purple"
      )
    })
    
    # infobox2 for CutMethod
    output$infobox2 <- renderInfoBox({
      infoBox(
        "Keep Top-N Variable Genes", paste0(ncol(exp.ds$table2)," / ",nrow(exp.ds$table)), icon = icon("scissors", lib = "glyphicon"),
        color = "lightblue"
      )
    })
    
    output$pre_choose_ui <- renderUI({
      radioGroupButtons(
        inputId = "pre_choose",
        #label = "Choose a graph :", 
        choices = c(`<i class='fa fa-pie-chart'></i>` = 'plot',
                    `<i class='fa fa-table'></i>` = 'table'),
        selected = "plot", justified = TRUE, width = "50", status = "info"
      )
    })
    
    output$infobox3 <- renderInfoBox({
      infoBox("Retained Samples", nrow(exp.ds$table2), icon = icon("vial"), color = "info")
    })
    
    # sample tree
    output$pre_clustPlot1 <- renderPlotly({
      if (is.null(exp.ds$data) || any(is.na(exp.ds$data)) || is.null(exp.ds$table2)) {
        return(NULL)
      }
      # Set a temporary graphics device to avoid opening a PDF file
      pdf(file = NULL)
      
      # Create ggplot object using ggdendrogram
      p <- ggdendrogram(exp.ds$sample_tree, rotate = FALSE, size= 4) +
        labs(title = "Sample clustering for outlier detection", x = "", y = "") +
        theme_minimal() +
        theme(          
          axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 12)
        )
      # Convert ggplot object to plotly object and set configurations
      plotly_obj <- ggplotly(p, width = input$preplot_width, height = input$preplot_height) %>%
        plotly::layout(
          title = list(
            text = "Sample clustering to detect outlier",  
            font = list(
              family = "Courier New, monospace",  
              size = 20,  
              color = "#7f7f7f"  
            )
          ),
          xaxis = list(
            tickangle = -75,  
            automargin = TRUE  
          )
        ) %>%
        config(toImageButtonOptions = list(format = input$preplot_format))
      # Close the temporary graphics device
      dev.off()
      
      return(plotly_obj)
    })
    
    # summary num
    output$pre_table = DT::renderDataTable({
      if(is.null(exp.ds$data)){return()}
      if(length(which(is.na(exp.ds$data))) != 0) {return()}
      df <- round(exp.ds$table2, digits = 2)
      datatable(t(df),class = 'compact stripe',
                # selection = list(mode = 'single', selected = c(1)),
                options = list(pageLength = 10, scrollX = TRUE, autoWidth = TRUE),
                colnames = c('GeneID' = 1))
    })
    
    preprocessing_summary_saved(
      isolate(preprocessing_summary_content())
    )
    
    shinyjs::show("step3_1_ui")
    shinyjs::show("predownload2_ui")
    shinyjs::show("pre_summary_ui")
    shiny::showTab(inputId = "pre_tabs", target = "outlier_tab")
    removeModal() 
  }
})

observeEvent(input$step3_1, {
  updateTabItems(session, "steps", "create")
  updateTabsetPanel(session, "create_tabs", "power_tab")
  print("Step 3-1 button clicked, proceed to STEP 3-1_Create network.")
})

###########################################################################
## step 2 : remove the outlier.

s_outlier = reactive({
  req(exp.ds$table2)
  rownames(exp.ds$table2)
})

observe({
  req(s_outlier)
  if(input$demo_on){
    updateSelectInput(session, "outlier", choices = s_outlier(), selected = "GSM387791")
  } else {
    updateSelectInput(session, "outlier", choices = s_outlier())
  }
})

observeEvent(input$remove, {
  
  if (length(input$outlier) == 0) {
    shinyalert("Error", "Please select the sample to remove", type = "error")
    return()
  }
  
  showModal(modalDialog(
    title = "Loading...",
    easyClose = FALSE,
    footer = NULL,
    div(id = "pre_remove_step1", tags$p("Loading data......")),
    div(id = "pre_remove_step2", shinyjs::hidden(tags$p("Ploting......")))
  ))
  
  exp.ds$outliersamples <- as.character(input$outlier)
  print(exp.ds$outliersamples)
  
  # 確保輸入的欄位存在於資料中
  if (!all(exp.ds$outliersamples %in% s_outlier())) {
    return()
  }
  
  # 使用 slicet() 函數移除選定的Row
  # 因為exp.ds$table2經過t()處理，所以要移除row而不是column
  exp.ds$table2 <<- exp.ds$table2 %>%
    dplyr::slice(-which(rownames(exp.ds$table2) %in% exp.ds$outliersamples))
  
  print(nrow(exp.ds$table2))
  
  output$pre_choose_remove_ui <- renderUI({
    radioGroupButtons(
      inputId = "rm_choose",
      #label = "Choose a graph :", 
      choices = c(`<i class='fa fa-pie-chart'></i>` = "plot",
                  `<i class='fa fa-bar-chart'></i>` = "table"),
      selected = "plot", justified = TRUE, width = "50", status = "info"
    )
  })
  
  # 重新畫圖
  shinyjs::show("pre_remove_step2")
  
  exp.ds$rm_sample_tree <- hclust(dist(exp.ds$table2), method = "average")
  
  output$pre_clustPlot2 = renderPlotly({
    if(is.null(exp.ds$table2)){return()}
    
    pdf(file = NULL)
    
    # 使用 ggdendrogram 生成 ggplot 物件
    gg_plot <- ggdendrogram(exp.ds$rm_sample_tree, rotate = FALSE, size = 2) +
      labs(title = "Sample clustering to detect outlier", x = "", y = "") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 12)
      )
    
    # Convert ggplot to plotly and configure labels
    rmplotly_obj <- ggplotly(gg_plot, width = input$rm_width, height = input$rm_height) %>%
      plotly::layout(
        title = list(
          text = "Sample clustering to detect outlier",  # Title text
          font = list(
            family = "Courier New, monospace",  # Font family
            size = 20,  # Font size
            color = "#7f7f7f"  # Font color
          )
        ),
        xaxis = list(
          tickangle = -75,  # Rotate tick labels
          automargin = TRUE  # Adjust margins to fit labels
        )
      ) %>%
      config(toImageButtonOptions = list(format = input$rm_format))
    
    dev.off()
    return(rmplotly_obj)
  })
  
  output$rm_table = DT::renderDataTable({
    if(is.null(exp.ds$table2)){return()}
    
    df2 <- round(exp.ds$table2, digits = 2)
    datatable(t(df2),class = 'compact stripe', 
              # selection = list(mode = 'single', selected = c(1)),
              options = list(pageLength = 10, scrollX = TRUE, autoWidth = TRUE),
              colnames = c('GeneID' = 1))
  })
  
  shinyjs::show("step2_2_remove_ui")
  shinyjs::show("pre_download4_ui")
  removeModal()
})

observeEvent(input$step2_2_remove, {
  updateTabItems(session, "steps", "create")
  print("Step 2-2 button clicked, proceed to STEP 3-1_Create network.")
})

###########################################################################
## download

output$predownload2 <- downloadHandler(
  filename = function() {
    paste("02_Normalizeddata.csv")
  },
  content = function(file) {
    df <- round(exp.ds$table2, digits = 2)
    
    df_transposed <- as.data.frame(t(df))    
    
    # 將"Gene_ID"作為第一列標題
    colnames(df_transposed)[1] <- "Gene_ID"   
    
    write.table(x = df_transposed , file = file, sep = ",",row.names = T,quote = F)
  }
)

output$predownload4 <- downloadHandler(
  filename = function() {
    paste("02_RemoveOutlierdata.csv")
  },
  content = function(file) {
    df2 <- round(exp.ds$table2, digits = 2)
    
    write.table(x = t(df2), file = file, sep = ",", quote = FALSE, row.names = T)
  }
)

