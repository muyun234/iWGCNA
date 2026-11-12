# server
## step 1 : Set Normalization set.
shiny::hideTab(inputId = "pre_tabs", target = "outlier_tab")

output$method_ui <- renderUI({
  if (input$raw_format == "normalized_count"){
    pickerInput(
      inputId = "normalized_method",
      label = tags$div(
        tags$span(
          style = "font-size: 18px;", "Normalization Method："
        ),
        tags$br(),
        tags$small(
          style = "font-size: 13px; color: #666;",
          "Choose the normalization method for the input expression data."
        )
      ),          
      choices = c("original value (FPKM, RPKM, TPM, or CPM)" = "raw", "logarithm, Log2(x+1)" = "logarithm"),
      inline = TRUE,
      options = list(style = "btn-info")
    )
  } else if (input$raw_format == "count"){
    pickerInput(
      inputId = "normalized_method",
      label = tags$div(
        tags$span(
          style = "font-size: 18px;", "Normalization Method："
        ),
        tags$br(),
        tags$small(
          style = "font-size: 13px; color: #666;",
          "Choose the normalization method for the input expression data."
        )
      ),
      choices = c("VST (DESeq2)" = "vst", "logarithm, Log2(x+1)" = "logarithm"),
      inline = TRUE,
      options = list(style = "btn-info")
    )
  }
})

output$cutoff_ui <- renderUI({
  if(input$raw_format == "normalized_count") {
    # numericInput(
    #   inputId = "expression_cutoff",
    #   label = tags$span(
    #     style = "font-size: 18px;", "Expression Cutoff ：", "",
    #     actionButton("ec_alt",label = "",icon = icon("info"), class = "btn-info btn-xs")
    #   ), 
    #   min = 0,  max = 100, value = 1, width = '100%'
    # )
  } else if (input$raw_format == "count") {
    numericInput(
      inputId = "expression_cutoff",
      label = tags$div(
        tags$span(
          style = "font-size: 18px;", "Low-Expression Gene Cutoff ："
        ),
        tags$br(),
        tags$small(
          style = "font-size: 13px; color: #666;",
          "Minimum expression value to be considered expressed (e.g., counts < 10)."
        )
      ), 
      min = 0, max = 100, value = 10, width = '100%'
    )
  }
})

observe({
  req(input$demo_on)
  if(input$demo_on){
    updatePickerInput(session, inputId = "normalized_method", selected = "raw")
    updateNumericInput(session, inputId = "sample_percentage", value = 0)
    updatePickerInput(session, inputId = "cutMethod", selected = "MAD")
    updateNumericInput(session, inputId = "reserved_gene_num", value = 21653)
    updatePrettyRadioButtons(session, inputId = "Networktype", selected = "unsigned")
    updateSliderInput(session, inputId = "cutoff", value = 0.9)
    updateSliderInput(session, inputId = "minMsize", value = 30)
    updateSliderInput(session, inputId = "mch", value = 0.25)
    
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
  shinyjs::hide()
  output$pre_clustPlot2 = renderPlotly({
    NULL
  })
  output$pre_table = DT::renderDataTable({
    NULL
  })

  print("Start preprocessing")
  
  method = input$normalized_method
  
  if (input$raw_format == "normalized_count") {
    RcCutoff <- 0
  } else if (input$raw_format == "count") {
    RcCutoff <- as.numeric(input$expression_cutoff)
  }  
  
  # 其他輸入
  samplePerc <- as.numeric(input$sample_percentage)
  # GNC <- as.numeric(input$reserved_gene_num)
  GNC <- as.integer(input$reserved_gene_num)
  cutmethod <- input$cutMethod
  datatype <- input$raw_format
  
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
        tibble::column_to_rownames("id") %>%
        dplyr::filter(
          rowSums(. > RcCutoff) > (samplePerc*ncol(.))
        )
      
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
    
    GeneNumCut <- 1-exp.ds$GNC/nrow(exp.ds$table)
    
    if(GeneNumCut == 0){
      datExprVar <- exp.ds$table
    } else if(cutmethod == "MAD"){
      m.mad <- apply(exp.ds$table, 1, mad)
      datExprVar <- exp.ds$table[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, GeneNumCut))[2], 0.01)),]
    } else if(cutmethod == "Var"){
      m.var <- apply(exp.ds$table, 1, var)
      datExprVar <- exp.ds$table[which(m.var > max(quantile(m.var, probs=seq(0, 1, GeneNumCut))[2], 0.01)),]
    } else {
      print("Error: wrong filter method")
    }
    
    dim(datExprVar)
    datExpr <- as.data.frame(t(datExprVar))
    gsg = WGCNA::goodSamplesGenes(datExpr, verbose = 3)
    if (!gsg$allOK){
      # Optionally, print the gene and sample names that were removed:
      if (sum(!gsg$goodGenes) > 0)
        printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ",")))
      if (sum(!gsg$goodSamples) > 0)
        printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ",")))
      # Remove the offending genes and samples from the data:
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
    
    shinyjs::show("step3_1_ui")
    shinyjs::show("predownload2_ui")
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

###################################################################
## show icon info.

# Choose normalized method when format === count 
# addPopover(
#   id = "nor",
#   options = list(
#     content = paste0("We recommend a variance-stabilizing transformation(vst) by DESeq2 ,
#                      but one could also start with normalized counts (e,g: RPKM/FPKM data) and 
#                      log-transform them using log10(x+1)."
#                      # br(),
#                      # "Any info is here",
#                      # a("WGCNA Q&A",href = "https://horvath.genetics.ucla.edu/cgi-sys/suspendedpage.cgi", target="_blank")
#     ),
#     title = "More information",
#     placement = "right",
#     trigger = "click",
#     html = TRUE  # Ensure HTML content is parsed
#   )
# )

# Input Expression cutoff when format === count
addPopover(
  id = "ec",
  options = list(
    content = paste0(
      "Recommend : ",br(),
      "For count : Remove genes with counts < 10 in more than 90% of samples. ",br(),
      "For normalized count : Remove genes with counts < 1 in more than 90% of samples."
      # br(),
      # "Any info is here",
      # a("WGCNA Q&A", href = "https://horvath.genetics.ucla.edu/cgi-sys/suspendedpage.cgi", target = "_blank")
    ),
    title = "More information",
    placement = "right",
    trigger = "click",
    html = TRUE  # Ensure HTML content is parsed
  )
)

# Choose normalized method when format !== count 
addPopover(
  id = "nor",
  options = list(
    content = paste0(
      "We recommend a variance-stabilizing transformation(vst) by DESeq2,
      and log-transform them using log10(x+1)."
      ),
    title = "More information",
    placement = "right",
    trigger = "click",
    html = TRUE  # Ensure HTML content is parsed
  )
)

# Input Expression cutoff foramt !=== count
addPopover(
  id = "ec_alt",
  options = list(
    content = paste0(
      "Recommend : ",br(),
      "For count : Remove genes with counts < 'Cutoff' in more than 'Fraction of Sample'. ",br(),
      "For normalized count : Remove genes with counts < 0 in more than 'Fraction of Sample'."
      # br(),
      # "Any info is here",
      # a("WGCNA Q&A", href = "https://horvath.genetics.ucla.edu/cgi-sys/suspendedpage.cgi", target = "_blank")
    ),
    title = "More information",
    placement = "right",
    trigger = "click",
    html = TRUE  # Ensure HTML content is parsed
  )
)

# Input Sample percentage
addPopover(
  id = "sp",
  options = list(
    content = paste0(
      "Recommend : ",br(),
      "For count : Remove genes with low expression across samples. ",br(),
      "For normalized count : Remove genes by the fraction of samples."
      # br(),
      # "Any info is here",
      # a("WGCNA Q&A", href = "https://horvath.genetics.ucla.edu/cgi-sys/suspendedpage.cgi", target = "_blank")
    ),
    title = "More information",
    placement = "right",
    trigger = "click",
    html = TRUE  # Ensure HTML content is parsed
  )
)

# Choose filter gene number method
addPopover(
  id = "cm",
  options = list(
    content = paste0(
      "Filter the top N genes with the greatest variation through median absolute deviation (MAD) ",
      "or based on variance (VAR) for subsequent analysis." 
      # br(),
      # "Any info is here",
      # a("WGCNA Q&A", href = "https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html", target = "_blank")
    ),
    title = "More information",
    placement = "right",
    trigger = "click",
    html = TRUE  # Ensure HTML content is parsed
  )
)

# Input Reserved gene number
addPopover(
  id = "rgn",
  options = list(
    content = paste0(
      "According to experimental needs and data characteristics to adjust.",
      br(),
      "If the number exceeds the count of genes that remain after the filtering process, 
      all the filtered genes will be retained.”"
    ),
    title = "More information",
    placement = "right",
    trigger = "click",
    html = TRUE  # Ensure HTML content is parsed
  )
)






