# server
# library(GEOquery)
# library(Biobase)
# library(limma)
library(tibble)
library(readr)
# devtools::install_github("jacobkap/fastDummies")
library(fastDummies)

shiny::hideTab(inputId = "geo_tabs", target = "trait_tab")

##################################
geonum = eventReactive(input$geo_num,{ 
  as.character(trimws(input$geo_num))
})
geo.data <- reactiveValues(
  gset = NULL
)
observeEvent(input$geo_nor, {
  shinyjs::toggle(id = "geo_nor_setting_ui", condition = input$geo_nor)
})

observeEvent(input$geo_update, {
  if(geonum() == "") {
    shinyalert("Please enter the GEO accession", type = "warning")
    return()
  }
  
  output$geo_error <- renderText({ "" })
  output$geo_platform_ui <- renderUI({""})
  shinyjs::hide("geo_platform_update_ui")
  shinyjs::hide("geo_nor_ui")
  shinyjs::hide("geo_next_step_ui")
  
  showModal(modalDialog(
    title = "Loading...", easyClose = FALSE, footer = NULL,
    div(id = "step1", tags$p("Downloading the database: ", geonum()))
  ))
  
  tryCatch({
    geo.data$gset <- GEOquery::getGEO(geonum(), GSEMatrix = TRUE, getGPL = TRUE)
    
    output$geo_platform_ui <- renderUI({
      print(paste0("number of platforms: ", length(geo.data$gset)))
      platforms <- sapply(geo.data$gset, annotation)
      geo.data$platforms <- setNames(names(platforms), platforms)
      
      selectInput("geo_platform", div(style = "font-size: 18px;", "Please Select the Platform: "), 
                  choices = geo.data$platforms, selected = geo.data$platforms[1])
    })
    
    shinyjs::show("geo_platform_update_ui")
    removeModal()
  },
  error = function(e) {
    removeModal()
    output$geo_error <- renderText({
      paste0("⚠️ error: ", e$message, "\nPlease check the GEO accession")
    })
  })
})

observeEvent(input$geo_platform_update, {
  showModal(modalDialog(
    title = "Loading...", easyClose = FALSE, footer = NULL,
    div(id = "step1", tags$p("Processing the data......"))
  ))
  
  print(input$geo_platform)
  selected_gset <- geo.data$gset[[input$geo_platform]]
  
  geo.data$expr <- Biobase::exprs(selected_gset)
  geo.data$feature <- Biobase::fData(selected_gset)
  geo.data$phen <- pData(selected_gset)
  
  if(nrow(geo.data$expr) != 0) {
    if(input$geo_nor){
      print("Normalize the expr data")
      geo.data$expr <- limma::normalizeBetweenArrays(geo.data$expr, method = input$geo_nor_method)
    }
    
    # 探針名轉為基因名
    gene_col <- grep("symbol", colnames(geo.data$feature), ignore.case = TRUE, value = TRUE)
    ID <- data.frame(ori_id = geo.data$feature$ID, gene_symbol = geo.data$feature[[gene_col[1]]])
    # 如果有重組基因，取第一個
    ID$gene_symbol <- data.frame(sapply(ID$gene_symbol, function(x)unlist(strsplit(x, "///"))[1]), stringsAsFactors = F)[, 1]
    exp <- as.data.frame(geo.data$expr)
    exp$ori_id <- rownames(exp)
    exp <- merge(exp, ID, by = "ori_id")
    # 去除基因名稱的前後空白字元
    exp[, grep("gene_symbol", colnames(exp))] <- trimws(exp[, grep("gene_symbol", colnames(exp))]) 
    exp[exp == ""] <- NA
    exp <- exp[!is.na(exp$gene_symbol), ]
    
    # 檢查有多少基因名稱重複
    table(duplicated(exp$gene_symbol))
    # 將重複名稱的表達量取平均值
    exp <- limma::avereps(exp, ID = exp$gene_symbol)
    exp <- as.data.frame(exp)
    cols_to_convert <- colnames(exp)[!colnames(exp) %in% c("ori_id", "gene_symbol")]
    exp <- exp %>%
      mutate(across(all_of(cols_to_convert), as.numeric))
    rownames(exp) <- exp$gene_symbol
    geo.data$filtered_exp <- exp %>%
      dplyr::select(-ori_id, -gene_symbol)
    geo.data$final_exp <- geo.data$filtered_exp %>%
      tibble::rownames_to_column(var = "gene_symbol")
  } else {
    print("No data available in expr")
    shinyalert(
      title = "Cannot retrieve expression data for this GEO series",
      text = "It is recommended to download the expression data from NCBI’s Supplementary Files.",
      type = "error"
    )
    removeModal()
    return()
  }
  
  output$geo_choose_ui <- renderUI({
    radioGroupButtons(
      inputId = "geo_choose",
      choices = c(`<i class='fa fa-table'></i>` = "table",
                  `<i class='fa fa-bar-chart'></i>` = "box_plot"),
      selected = "table", justified = TRUE, width = "80", status = "info"
    )
  })
  
  output$geo_boxplot <- renderPlot({
    exp_matrix <- as.matrix(geo.data$filtered_exp)
    par(mar=c(7,4,2,1))
    title <- paste (geonum(), "/", input$geo_platform, sep = "")
    boxplot(exp_matrix, boxwex = 0.7, notch = T, main = title, outline = FALSE, las = 2)
  })
  
  output$geo_exp <- renderDT({
    datatable(
      geo.data$final_exp, selection = "none",
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })
  
  shinyjs::show("geo_nor_ui")
  shinyjs::show("geo_next_step_ui")
  shiny::showTab(inputId = "geo_tabs", target = "trait_tab")
  
  # trait
  rownames(geo.data$phen) <- NULL
  names(geo.data$phen) <- gsub("[:_]ch1", "", names(geo.data$phen)) 
  names(geo.data$phen) <- gsub("\\s", "_", names(geo.data$phen)) 
  
  output$geo_phen <- renderDT({
    datatable(
      geo.data$phen, selection = list(mode = "multiple", target = "column"), 
      extensions = c('Buttons', "ColReorder"),  escape = FALSE,
      options = list(dom = 'Blfrtip', buttons = c('colvis'), 
                     colReorder = TRUE, pageLength = 10, scrollX = TRUE, scrollY = "400px",
                     columnDefs = list(
                       list(
                         targets = "_all",
                         render = JS(
                           "function(data, type, row, meta) {",
                           "return type === 'display' && data.length > 15 ?",
                           "'<span title=\"' + data + '\">' + data.substr(0, 15) + '...</span>' : data;",
                           "}") 
                       )
                     )
      )
    )
  })
  
  shinyjs::show("geo_download3_ui")
  removeModal()
})

observeEvent(input$geo_next_step, {
  updateTabsetPanel(session, "geo_tabs", "trait_tab")
})

observeEvent(input$geo_generate_filtered_data, {
  if(is.null(input$geo_phen_columns_selected)) {
    shinyalert("Please select aleast one column", type = "warning")
    return()
  }
  
  if ("geo_accession" %in% names(geo.data$phen)) {
    geo.data$sample_colname <- "geo_accession"
  } else {
    # 偵測標註樣本名稱的欄位名
    gsm_pattern <- "^GSM[0-9]+$"
    candidate_cols <- sapply(geo.data$phen, function(col) any(grepl(gsm_pattern, as.character(col))))
    
    if (any(candidate_cols)) {
      geo.data$sample_colname <- names(geo.data$phen)[candidate_cols]
    } else {
      geo.data$sample_colname <- NULL
    }
  }
  
  selected_cols <- input$geo_phen_columns_selected
  geo.data$trait <- geo.data$phen[, selected_cols, drop = FALSE]
  
  if (!is.null(geo.data$sample_colname) && geo.data$sample_colname %in% colnames(geo.data$trait)) {
    print("Rearrange")
    geo.data$trait <- geo.data$trait %>%
      dplyr::select(geo.data$sample_colname, everything())    
  }
  
  output$geo_trait <- renderDT({
    datatable(
      geo.data$trait, selection = "none", extensions = c('Buttons'),
      options = list(
        colReorder = TRUE, pageLength = 10, scrollX = TRUE, scrollY = "400px"
      )
    )
  })
  shinyjs::show("geo_download4_ui")
})

observeEvent(input$geo_convert_traits, {
  
  choices <- base::setdiff(colnames(geo.data$final_exp), c("gene_symbol"))
  output$geo_rm_sample_ui <- renderUI({
    selectInput("geo_rm_sample", "Remove the Sample: ", multiple = TRUE,
                choices = choices, selected = NULL)
  })
  
  cols_to_dummy <- base::setdiff(names(geo.data$trait), geo.data$sample_colname)
  
  geo.data$trans_trait <- fastDummies::dummy_cols(geo.data$trait, select_columns = cols_to_dummy, 
                                                  remove_selected_columns = TRUE)
  output$geo_transformed_trait <- renderDT({
    datatable(
      geo.data$trans_trait, selection = "none", extensions = c('Buttons'),
      options = list(
        colReorder = TRUE, pageLength = 10, scrollX = TRUE, scrollY = "400px"
      )
    )
  })
  
  shinyjs::show("geo_download5_ui")
})

observeEvent(input$geo_rm_update, {
  if(is.null(input$geo_rm_sample)){
    shinyalert("Please select the sample names you want to remove", type = "warning")
    return()
  }
  
  if (is.null(geo.data$sample_colname) || !geo.data$sample_colname %in% names(geo.data$trans_trait)) {
    shinyalert("The selected column does not exist the column containing GEO sample IDs", type = "warning")
    return()
  }
  
  geo.data$trans_trait <- geo.data$trans_trait[
    !geo.data$trans_trait[[geo.data$sample_colname]] %in% input$geo_rm_sample, 
  ]
  
  output$geo_transformed_trait <- renderDT({
    datatable(
      geo.data$trans_trait, selection = "none", extensions = c('Buttons'),
      options = list(
        colReorder = TRUE, pageLength = 10, scrollX = TRUE, scrollY = "400px"
      )
    )
  })
})

##########################
# download
observeEvent(input$geo_platform_update, {
  req(geo.data$final_exp)
  req(geo.data$filtered_exp)
  req(geo.data$phen)
  
  output$geo_download <- downloadHandler(
    filename = function() {
      paste0(geonum(), ".txt")
    },
    content = function(file) {
      write_tsv(geo.data$final_exp, file)
    }
  )
  
  output$geo_download2 <- renderUI({
    screenshotButton(
      label = "Download Plot",
      download = TRUE, 
      server_dir = tempdir(), 
      id = "geo_boxplot", 
      filename = paste0(geonum(), "_boxplot"))
  })
  
  output$geo_download3 <- downloadHandler(
    filename = function() {
      paste0(geonum(), "_oriphen.txt")
    },
    content = function(file) {
      write_tsv(geo.data$phen, file)
    }
  )
})

observeEvent(input$geo_next_step, {
  output$geo_download4 <- downloadHandler(
    filename = function() {
      paste0(geonum(), "_trait.txt")
    },
    content = function(file) {
      write_tsv(geo.data$trait, file)
    }
  )
})

observeEvent(input$geo_convert_traits, {
  req(geo.data$trans_trait)
  
  output$geo_download5 <- downloadHandler(
    filename = function() {
      paste0(geonum(), "_transformed_trait.txt")
    },
    content = function(file) {
      write_tsv(geo.data$trans_trait, file)
    }
  )
})
##########################
# info
addPopover(
  id = "nor_method_info",
  options = list(
    content = paste0(
      "Recommend: ",br(),
      "For quantile: Common, makes sample distributions identical",br(),
      "For scale: Quick, use when sample differences are small", br(),
      "For cyclicloess: Fine correction for large sample differences", br(),
      "⚠️ Only for single-channel data"
    ),
    title = "More information",
    placement = "right",
    trigger = "click",
    html = TRUE 
  )
)
