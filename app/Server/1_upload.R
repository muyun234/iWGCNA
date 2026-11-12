#server
exp.ds <- reactiveValues(
  data = NULL,
  genename_type = NULL,
  ori_gene_number = NULL,
  switch_genename = FALSE,
  filtering_gene_number = NULL,
  sft = NULL,
  netout = NULL,
  trait = NULL,
  network_list = list(),
  organism = NULL,
  KME = NULL,
  hubgene_df = NULL
)
data_folder <- reactiveValues(folder_path = NULL)

organism_list <- list(
  # 哺乳類
  "Homo sapiens" = list(code = "hsa", orgnum = 9606),
  "Mus musculus" = list(code = "mmu", orgnum = 10090),
  "Rattus norvegicus" = list(code = "rno", orgnum = 10116),
  "Sus scrofa" = list(code = "ssc", orgnum = 9823),
  "Bos taurus" = list(code = "bta", orgnum = 9913),
  
  # 其他脊椎動物
  "Danio rerio" = list(code = "dre", orgnum = 7955),
  
  # 無脊椎動物
  "Drosophila melanogaster" = list(code = "dme", orgnum = 7227),
  "Caenorhabditis elegans" = list(code = "cel", orgnum = 6239),
  
  # 植物
  "Arabidopsis thaliana" = list(code = "ath", orgnum = 3702),
  "Glycine max" = list(code = "gmx", orgnum = 3847),
  "Oryza sativa Japonica" = list(code = "osa", orgnum = 39947),
  "Zea mays" = list(code = "zma", orgnum = 4577),
  
  # 真菌
  "Saccharomyces cerevisiae" = list(code = "sce", orgnum = 4932),
  
  # 原生生物
  "Plasmodium falciparum" = list(code = "pfa", orgnum = 36329),
  # "T. vaginalis G3 2022 (ref)" = list(code = "tva", orgdb = "org.Tvaginalis.eg.db", orgnum = 5722),
  "T. vaginalis G3" = list(code = "tva_g3", orgnum = "5722_g3")
)

output$select_organism_ui <- renderUI({
  pickerInput(
    inputId = "select_organism",
    label = tags$span(style = "font-size: 18px;", "Specify the organism:"),
    choices = sapply(organism_list, `[[`, "code"), selected = "hsa",
  )
})

observe({
  if (input$demo_on){
    exp.ds$organism <- "hsa"
  }else {
    exp.ds$organism <- input$select_organism
  }
})

observe({
  req(exp.ds$organism)
  
  org_name <- names(organism_list)[
    sapply(organism_list, `[[`, "code") == exp.ds$organism
  ]
  
  org_info <- organism_list[[org_name]]
  
  # exp.ds$orgdb_name <- org_info$orgdb
  # requireNamespace(exp.ds$orgdb_name, quietly = TRUE)
  # exp.ds$orgdb <- get(exp.ds$orgdb_name, envir = asNamespace(exp.ds$orgdb_name))
  
  exp.ds$orgnum <- if (!is.null(org_info$orgnum)) {
    org_info$orgnum
  } else {
    NA
  }
})

output$exp_file_ui <- renderUI({
  fileInput(
    "exp_file", label = div(style = "font-size: 18px;", "Please select a file or drag it into the upload area"),
    buttonLabel = tags$i(class = "fa-solid fa-upload", style = "font-size: 12px"),
    placeholder = "accept format .txt, .tsv, .csv; max: 30 Mb",
    accept = c(".txt", ".tsv", ".csv")
  )
})

# Demo
observeEvent(input$demo_on, {
  if (input$demo_on) {
    showModal(modalDialog(title = "Loading...", easyClose = FALSE, footer = NULL))
    shinyalert("Success", "Demo gene expression file has been successfully loaded.", type = "success")
    
    exp.ds$switch_genename <- FALSE
    
    output$exp_file_ui <- renderUI({
      fileInput(
        "exp_file", 
        label = div(style = "font-size: 18px;", "Please select a file or drag it into the upload area"),
        buttonLabel = tags$i(class = "fa-solid fa-upload", style = "font-size: 12px"),
        placeholder = "accept format .txt, .tsv, .csv; max: 30 Mb",
        accept = c(".txt", ".tsv", ".csv")
      )
    })
    
    # Read the demo file
    file1 <- read.delim(file = "Demo/GSE15459.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    
    exp.ds$data <- file1
    exp.ds$ori_gene_number <- nrow(file1)
    exp.ds$preview <- exp.ds$data

    # gene count and sample count
    exp.ds$gene_sample_count <- TRUE
    exp.ds$switch_genename <- TRUE
    exp.ds$genename_type = "SYMBOL"
    
    shinyjs::show("inter_demo_ui")
    shinyjs::hide("run_ui")
    shinyjs::hide("information_ui")
    shinyjs::show("preview_ui")
    shinyjs::show("step2_ui")
    
    updatePrettyRadioButtons(session, inputId = "raw_format", selected = "normalized_count")
    updatePickerInput(session, inputId = "select_organism", selected = "hsa")
    
    removeModal()
  } else if (!input$demo_on){
    exp.ds$data <- NULL
    exp.ds$gene_sample_count <- FALSE
    exp.ds$switch_genename <- FALSE
    shinyjs::hide("preview_ui")
    shinyjs::show("information_ui")
    shinyjs::hide("step2_ui")
    exp.ds$show_trait_df <- FALSE
  }
})

# data
observe({
  if (is.null(input$exp_file)) {
    shinyjs::hide("run_ui")
  } else {
    updateAwesomeCheckbox(session, inputId = "demo_on", value = FALSE)
    shinyjs::show("run_ui")
  }
})

observeEvent(input$data_upload, {
  if (is.null(input$exp_file)) {
    shinyalert("Error", "No file selected.", type = "error")
    return()
  }
  
  showModal(modalDialog(title = "Loading...", easyClose = FALSE, footer = NULL))
  
  exp.ds$gene_sample_count <- FALSE
  exp.ds$switch_genename <- FALSE

  ext <- tools::file_ext(input$exp_file$name)
  if (ext != "txt" && ext != "tsv" && ext != "csv") {
    shinyalert("Error", "Please upload a .txt, .tsv or csv file.", type = "error")
  } else {
    shinyalert("Success", paste(input$exp_file$name, "has been successfully uploaded."), type = "success")
  }
  
  file1 <- input$exp_file
  if(ext %in% c("txt", "tsv")) {
    df <- read.delim(file = file1$datapath, sep = "\t", header = TRUE, 
                     stringsAsFactors = FALSE, check.names = FALSE)
  } else if(ext == "csv") {
    df <- read.csv(input$exp_file$datapath, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  if (all(!is.na(as.numeric(df[[1]]))) && is.null(rownames(df[[1]]))) {
    df <- df[, -1, drop = FALSE] 
  }
  
  if (!is.character(df[[1]]) && !is.null(rownames(df))) {
    print("Change rowname into GeneName col")
    df <- cbind(GeneName = rownames(df), df)
    rownames(df) <- NULL
  }
  
  df <- as.data.frame(
    lapply(df, function(x) if (is.integer(x)) as.numeric(x) else x),
    check.names = FALSE
  )
  if (any(is.na(df))) {
    df[is.na(df)] <- 0
  }
  
  exp.ds$ori_gene_number <- nrow(df)
  exp.ds$data <- df
  
  exp.ds$gene_sample_count <- TRUE
  
  shinyjs::show("step2_ui")
  # updatePrettyRadioButtons(session, inputId = "raw_format", selected = "count")
  
  shinyjs::hide("information_ui")
  shinyjs::show("preview_ui")
  
  # change gene id into ensemble
  exp.ds$switch_genename <- TRUE
  exp.ds$preview <- exp.ds$data
  
  shinyjs::hide("inter_demo_ui")
  removeModal()
})

output$raw_data <- DT::renderDataTable({
  req(exp.ds$preview)
  
  df <- exp.ds$preview     
  # 只對數值類型進行四捨五入
  numeric_cols <- sapply(df, is.numeric)
  df[numeric_cols] <- round(df[numeric_cols], digits = 2)
  print(str(df))
  
  datatable(
    df,class = 'compact stripe', selection = "none",
    options = list(pageLength = 5, scrollX = TRUE, autoWidth = F)
  ) 
})

output$gene_count_ui <- renderUI({
  req(exp.ds$gene_sample_count, exp.ds$data)
  
  data_nrow <- nrow(exp.ds$data)
  infoBox(
    tabName = "gene_count",
    title = "Gene Count",
    width = NULL,
    value = data_nrow,
    color = "info",
    icon = icon("hashtag")
  )
})

output$sample_count_ui <- renderUI({
  req(exp.ds$gene_sample_count, exp.ds$data)
  
  sample_ncol <- ncol(exp.ds$data)
  infoBox(
    tabName = "sample_count",
    title = "Sample Count",
    width = NULL,
    value = sample_ncol - 1,
    color = "lightblue",
    icon = icon("file")
  )
})

observeEvent(input$step2, {
  req(exp.ds$data, cancelOutput = TRUE)
  
  print("Step 2 button clicked, proceed to STEP 2_Data Preprocessing.")
  updateTabItems(session, "steps", "preprocessing")
  shiny::hideTab(inputId = "pre_tabs", target = "outlier_tab")
})

########################################
# Download
output$download <- downloadHandler(
  filename = function() {
    paste("Demo_GSE15459.tsv")
  },
  content = function(file) {
    file.copy("Demo/GSE15459.tsv", file)
  }
)
########################################
# show icon info.
addPopover(
  id = "rf",
  options = list(
    content = "Normalized Read Count: RPKM, FPKM, TPM, CPM",
    title = "More information",
    placement = "right",
    trigger = "click",
    html = TRUE  # Ensure HTML content is parsed
  )
)





