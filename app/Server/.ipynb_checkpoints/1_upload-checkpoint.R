#server
exp.ds <- reactiveValues(
  data = NULL,
  genename_type = NULL,
  ori_gene_number = NULL,
  filtering_gene_number = NULL,
  sft = NULL,
  netout = NULL
  )

# Helper function to handle file operations
handleFileUpload <- function(file_path, destination_folder, file_name) {
  if (!dir.exists(destination_folder)) {
    dir.create(destination_folder)
  }
  file.copy(file_path, paste0(destination_folder, "/", file_name))
}

# Demo
observeEvent(input$demo_on, {
  if (input$demo_on) {
    shinyjs::hide("upload_ui")
    shinyjs::hide("data_upload_ui")
    shinyjs::hide("gene_filter_number_ui")
    shinyjs::hide("step2_ui")
    shinyalert("Success", "Demo gene expression file has been successfully loaded.", type = "success")
    
    # Read the demo file
    file1 <- read.delim(file = "Demo/normalized_GSE51588.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    exp.ds$data <- file1
    exp.ds$ori_gene_number <- nrow(file1)

    shinyjs::show("gene_type_ui")
    shinyjs::show("raw_data_ui")
    shinyjs::show("step2_ui")
    shinyjs::show("dataset_properties_ui")
    updatePrettyRadioButtons(session, inputId = "raw_format", selected = "normalized_count")
    updatePrettyRadioButtons(session, inputId = "gene_type", selected = character(0))
  } else {
    shinyjs::show("upload_ui")
    exp.ds$data <- NULL
    shinyjs::hide("gene_type_ui")
    shinyjs::hide("raw_data_ui")
    shinyjs::hide("step2_ui")
    shinyjs::hide("dataset_properties_ui")
    shinyjs::hide("gene_filter_number_ui")
  }
})

# data
observeEvent(input$exp_file, {
  if (!is.null(input$exp_file)) {
    shinyjs::show("data_upload_ui")
  }
})

observeEvent(input$data_upload, {
  if (is.null(input$exp_file)) {
    shinyalert("Error", "No file selected.", type = "error")
    return()
  }
  
  shinyjs::hide("gene_filter_number_ui")
  
  ext <- tools::file_ext(input$exp_file$name)
  if (ext != "txt" && ext != "tsv") {
    shinyalert("Error", "Please upload a .txt or .tsv file.", type = "error")
  } else {
    shinyalert("Success", paste(input$exp_file$name, "has been successfully uploaded."), type = "success")
  }
  
  file1 <- input$exp_file
  df <- read.delim(file = file1$datapath, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  df <- as.data.frame(lapply(df, function(x) {
    if (is.integer(x)) {
      as.numeric(x)
    } else {
      x
    }
  }))
  if (any(is.na(df))) {
    df[is.na(df)] <- 0
  }
  exp.ds$data <- df
  exp.ds$ori_gene_number <- nrow(df)
  
  shinyjs::show("dataset_properties_ui")
  shinyjs::show("gene_type_ui")
  shinyjs::show("raw_data_ui")
  shinyjs::show("step2_ui")
  updatePrettyRadioButtons(session, inputId = "raw_format", selected = "count")
  updatePrettyRadioButtons(session, inputId = "gene_type", selected = character(0))
})


output$raw_data = DT::renderDataTable({
  req(exp.ds$data)
  
  df <- exp.ds$data      
  
  # 只對數值類型進行四捨五入
  numeric_cols <- sapply(df, is.numeric)
  df[numeric_cols] <- round(df[numeric_cols], digits = 2)
  
  print(str(df))
  
  datatable(df,class = 'compact stripe', 
            selection = "none",
            options = list(pageLength = 5, scrollX = TRUE, autoWidth = F)
            ) 
})

observeEvent(input$gene_type, {
  req(exp.ds$data, cancelOutput = TRUE)
  
  showModal(modalDialog(
    title = "Loading...", easyClose = FALSE, footer = NULL
  ))
  
  exp.ds$filtering_gene_number <- NULL
  shinyjs::show("gene_filter_number_ui")
  
  if(input$gene_type == "Ensembl ID"){
    df <- exp.ds$data
    names(df)[1] <- "OriginalID"
    trans_id <- df$OriginalID[startsWith(df$OriginalID, "ENST")]
    symbol <- df$OriginalID[!startsWith(df$OriginalID, "ENST") & !startsWith(df$OriginalID, "ENSG")]
    
    if(length(trans_id) > 0){
      trans2ensembl <- mapIds(org.Hs.eg.db, keys = trans_id,  keytype = "ENSEMBLTRANS", column="ENSEMBL")
      trans_data <- df[startsWith(df$OriginalID, "ENST"),]
      trans_data$OriginalID <- trans2ensembl
      trans_data <- trans_data[!is.na(trans_data$OriginalID), ]
    } else {
      trans_data <- NULL
    }
    ensembl_data <- df[startsWith(df$OriginalID, "ENSG"),]
    if(length(symbol) > 0){
      symbol2ensembl <- mapIds(org.Hs.eg.db, keys = symbol, keytype = "SYMBOL", column = "ENSEMBL")
      symbol_data <- df[!startsWith(df$OriginalID, "ENST") & !startsWith(df$OriginalID, "ENSG"),]
      symbol_data$OriginalID <- symbol2ensembl
      symbol_data <- symbol_data[!is.na(symbol_data$OriginalID), ]
    } else {
      symbol_data <- NULL
    }
    final_data <- rbind(trans_data, ensembl_data, symbol_data)
    final_data <- final_data %>%
      group_by(OriginalID) %>%
      summarise(across(where(is.numeric), ~mean(., na.rm = TRUE))) %>%
      ungroup()
  } else if(input$gene_type == "Symbol"){
    
    df <- exp.ds$data
    names(df)[1] <- "OriginalID"
    trans_id <- df$OriginalID[startsWith(df$OriginalID, "ENST")]
    ensembl_ids <- df$OriginalID[startsWith(df$OriginalID, "ENSG")]

    if(length(trans_id) > 0){
      trans2symbol <- mapIds(org.Hs.eg.db, keys = trans_id,  keytype = "ENSEMBLTRANS", column="SYMBOL")
      trans_data <- df[startsWith(df$OriginalID, "ENST"),]
      trans_data$OriginalID <- trans2symbol
      trans_data <- trans_data[!is.na(trans_data$OriginalID), ]
    } else {
      trans_data <- NULL
    }
    if(length(ensembl_ids) > 0){
      ensembl2symbol <- mapIds(org.Hs.eg.db, keys = ensembl_ids,  keytype = "ENSEMBL", column="SYMBOL")
      ensembl_data <- df[startsWith(df$OriginalID, "ENSG"),]
      ensembl_data$OriginalID <- ensembl2symbol
      ensembl_data <- ensembl_data[!is.na(ensembl_data$OriginalID), ]
    } else {
      ensembl_data <- NULL
    }
    symbol_data <- df[!startsWith(df$OriginalID, "ENST") & !startsWith(df$OriginalID, "ENSG"),]
    
    final_data <- rbind(trans_data, ensembl_data, symbol_data)
    final_data <- final_data %>%
      group_by(OriginalID) %>%
      summarise(across(where(is.numeric), ~mean(., na.rm = TRUE))) %>%
      ungroup()
  }
  
  names(final_data)[1] <- "GeneName"
  exp.ds$data <- final_data
  exp.ds$filtering_gene_number <- nrow(final_data)
  
  output$gene_num_ui <- renderUI({
    infoBox(
      tabName = "gene_number",
      title = "Total gene number",
      width = NULL ,
      value = paste0(exp.ds$filtering_gene_number, " / ", exp.ds$ori_gene_number),
      color = "info",
      icon = icon("laptop-code")
    )
  })
  
  shinyalert("Success", paste("Gene IDs have been converted to", input$raw_format), type = "success")
  removeModal()
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

observeEvent(input$step2, {
  req(exp.ds$data, cancelOutput = TRUE)
  
  # 檢查基因名稱的格式是否統一
  df <- exp.ds$data
  names(df)[1] <- "GeneName"
  gene_names <- exp.ds$data$GeneName
  
  is_ensembl_format <- startsWith(gene_names, "ENSG")
  is_trans_format <- startsWith(gene_names, "ENST")
  is_symbol_format <- !is_ensembl_format & !is_trans_format
  
  if (all(is_ensembl_format)){
    showNotification("Gene name format check: All gene names are in Ensembl format.", type = "message")
    exp.ds$genename_type = "ENSEMBL"
  } else if (all(is_symbol_format)) {
    showNotification("Gene name format check: All gene names are in Symbol format.", type = "message")
    exp.ds$genename_type = "SYMBOL"
  } else {
    shinyalert("Gene names inconsistent! Please ensure all Ensembl IDs or all symbols",type = "error")
    return(NULL)
  }

  # Get the path of the temporary directory
  temp_directory <- tempdir()
  
  current_time <- format(Sys.time(), "%Y-%m-%d-%H-%M-%S")
  unzip_folder <<- file.path(temp_directory,paste0("Save_", current_time))
  
  if(input$demo_on){
    # Save the demo file to the SAVE folder
    file.copy("Demo/normalized_GSE51588.txt", file.path(unzip_folder))
    
    print("Step 2 button clicked, proceed to STEP 2_Data Preprocessing (demo).")
  } else {
    # process exp_file
    exp_file_path <- input$exp_file$datapath
    # copy file_ext
    ext <- tools::file_ext(input$exp_file$name)
    if (ext == "tsv" ||ext == "txt") {
      handleFileUpload(exp_file_path, unzip_folder, input$exp_file$name)
    } else {
      showNotification("Unsupported file format for exp_file", type = "error")
    }
    
    print("Step 2 button clicked, proceed to STEP 2_Data Preprocessing.")
    
  }
  
  updateTabItems(session, "steps", "preprocessing")
  shiny::hideTab(inputId = "pre_tabs", target = "outlier_tab")
  })

########################################
# show icon info.
addPopover(
  id = "rf",
  options = list(
    content = "Count : Read Count<br>Normalized count : RPKM,FPKM,TPM,CPM",
    title = "More information",
    placement = "right",
    trigger = "click",
    html = TRUE  # Ensure HTML content is parsed
  )
)





