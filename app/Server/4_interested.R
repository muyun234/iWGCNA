# server
shiny::hideTab(inputId = "inter_tabs", target = "inter_association_tab")

##################################################################
# step 4a : upload trait ####
phen <- reactiveVal(NULL) 

processTraitData <- function() {
  req(phen())
  
  print("Processing trait data...")
  # Process uploaded trait data
  if (ncol(phen()) == 2) {
    x <- phen()
    Tcol <- as.character(unique(x[, 2]))
    b <- list()
    for (i in 1:length(Tcol)) {
      b[[i]] <- data.frame(row.names = x[, 1],
                           levels = ifelse(x[, 2] == Tcol[i], 1, 0))
    }
    c <- dplyr::bind_cols(b)
    c <- data.frame(row.names = x[, 1], c)
    colnames(c) <- Tcol
    rownames(c) <- phen()[, 1]
    exp.ds$phen <- c
  } else {
    exp.ds$phen <- data.frame(row.names = phen()[, 1], phen()[, -1])
  }
  exp.ds$phen <- exp.ds$phen[match(rownames(exp.ds$table2), rownames(exp.ds$phen)), ]
}

output$inter_upload_ui <- renderUI({
  fileInput(
    "trait_file",
    buttonLabel = tags$i(class = "fa-solid fa-upload", style = "font-size: 12px"),
    label = tags$p(style = "font-size: 18px;", strong("Upload Your Trait/Phenotype Data")),
    placeholder = "Only accepted .txt, .tsv or .csv file format.;max: 30MB",
    accept = c(".txt",".tsv", ".csv")
  )
})

observeEvent(input$inter_demo, {
  if (input$inter_demo) {
    shinyjs::hide("inter_upload_btn_ui")
    
    if(is.null(exp.ds$net)){
      shinyalert("Error", "Please preprocess the gene group detection", type = "warning")
      updateAwesomeCheckbox(session, "inter_demo", value = FALSE)
      updateTabItems(session, "steps", "create")
      return()
    }

    output$inter_upload_ui <- renderUI({
      fileInput(
        "trait_file",
        buttonLabel = tags$i(class = "fa-solid fa-upload", style = "font-size: 12px"),
        label = div(style = "font-size: 18px;", "Choose a file or drag it here"),
        placeholder = "Only accepted .txt, .tsv or .csv file format.;max: 30MB",
        accept = c(".txt",".tsv", ".csv")
      )
    })
    
    df <- read.csv("Demo/GSE15459_trait.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    phen(df)
    processTraitData()
    shinyjs::show("inter_download4_ui")
    
    shinyalert("Success", "Demo trait file has been successfully loaded.", type = "success")
    shinyjs::show("step4b_ui")
    exp.ds$show_trait_df <- TRUE
    
  }else if (!input$inter_demo){
    if(is.null(input$trait_file)){
      shinyjs::hide("step4b_ui")
      phen(NULL)
      exp.ds$show_trait_df <- FALSE
    } else {
      shinyjs::show("step4b_ui")
    }
  }
})


observeEvent(input$trait_file, {
  if (!is.null(input$trait_file)) {
    shinyjs::show("inter_upload_btn_ui")
  }
})

observeEvent(input$inter_data_upload, {
  if(is.null(exp.ds$net)){
    shinyalert("Error", "Please preprocess the gene group detection", type = "warning")
    updateTabItems(session, "steps", "create")
    return()
  }
  if (is.null(input$trait_file)) {
    shinyalert("Error", "No file selected", type = "error")
    return()
  }

  ext <- tools::file_ext(input$trait_file$name)
  if (ext != "txt" && ext != "tsv" && ext != "csv") {
    shinyalert("Error", "Please upload a .txt, .tsv or .csv file.", type = "error")
    return()
  } else {
    shinyalert("Success", paste(input$trait_file$name, "has been successfully uploaded."), type = "success")
  }
  
  updateAwesomeCheckbox(session = session, inputId = "inter_demo", value = FALSE)

  file_path <- input$trait_file$datapath
  
  if(ext %in% c("txt", "tsv")) {
    df <- read.delim(file = file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  } else if(ext == "csv") {
    df <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }
  phen(df)
  processTraitData()
  shinyjs::show("inter_download4_ui")
  
  shinyjs::show("step4b_ui")
  exp.ds$show_trait_df <- TRUE
  })

observeEvent(input$step4b, {
  updateTabsetPanel(session = session, inputId = "inter_tabs", selected = "inter_association_tab")
  print("Step 4b button clicked, proceed to 4b_Compute Associations.")
})

output$inter_upload_check <- renderDT({
  req(exp.ds$phen, exp.ds$show_trait_df)
  datatable(
    exp.ds$phen, selection = "none",
    options = list(
      dom = 'lfrtip', paging = TRUE, ordering = FALSE, 
      scrollY = "210px", scrollCollapse = TRUE
      )
  )
})

##################################################################
# 4b correlations
observeEvent(input$step4b, {
  if (is.null(exp.ds$phen)) {
    shinyalert("Error", "Please upload a trait file", type = "error")
    return()
  }
  
  showModal(modalDialog(
    title = "Loading...", easyClose = FALSE, footer = NULL,
    div(id = "inter_step1", tags$p("Calculate the correlation.......")),
    div(id = "inter_step2", shinyjs::hidden(tags$p("Ploting......")))
  ))
  shinyjs::show("inter_cor_ui")
  
  # DEBUG
  print(head(exp.ds$phen))
  
  # copy file to folder
  if(input$inter_demo){
    file.copy("Demo/GSE326939_trait.csv", file.path(data_folder$folder_path))
  } else {
    file4_path <- file.path(data_folder$folder_path, input$trait_file$name)
    file.copy(input$trait_file$datapath, file4_path)
  }
  
  MEs_order <- WGCNA::orderMEs(exp.ds$MEs_col)
  moduleTraitCor <- cor(MEs_order, exp.ds$phen, use = "p")
  moduleTraitPvalue <- WGCNA::corPvalueStudent(moduleTraitCor, exp.ds$nSamples)
  
  textMatrix <-  paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "");
  
  stars <- matrix(
    "", nrow = nrow(moduleTraitPvalue),
    ncol = ncol(moduleTraitPvalue)
  )
  stars[!is.na(moduleTraitPvalue) & moduleTraitPvalue < 0.05] <- "*"
  stars[!is.na(moduleTraitPvalue) & moduleTraitPvalue < 0.01] <- "**"
  stars[!is.na(moduleTraitPvalue) & moduleTraitPvalue < 0.001] <- "***"
  
  textMatrix <- matrix(
    paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", stars),
    nrow = nrow(moduleTraitCor),
    ncol = ncol(moduleTraitCor),
    dimnames = dimnames(moduleTraitCor)
  )
  
  exp.ds$moduleTraitCor <- moduleTraitCor
  exp.ds$moduleTraitPvalue <- moduleTraitPvalue
  exp.ds$textMatrix <- textMatrix
  
  if (is.null(exp.ds$KME)) {
    print("calculate KME")
    exp.ds$KME <- WGCNA::signedKME(exp.ds$table2, exp.ds$MEs_col, outputColumnName="")
  }
  
  exp.ds$mod_color = gsub(pattern = "^..", replacement = "", rownames(exp.ds$moduleTraitCor))
  exp.ds$mod_color_anno = setNames(exp.ds$mod_color, exp.ds$mod_color)
  exp.ds$Left_anno = rowAnnotation(
    Module = exp.ds$mod_color,
    col = list(Module = exp.ds$mod_color_anno),
    show_legend = FALSE,
    show_annotation_name = FALSE
  )
  
  # 修改 moduleTraitCor 的行名
  rownames(exp.ds$moduleTraitCor) <- exp.ds$mod_color
  
  # 修改 moduleTraitPvalue 的行名
  rownames(exp.ds$moduleTraitPvalue) <- exp.ds$mod_color
  
  # 修改 textMatrix 的行名
  rownames(exp.ds$textMatrix) <- exp.ds$mod_color
  
  file3_path <- file.path(data_folder$folder_path, "trait_textMatrix.rds")
  saveRDS(exp.ds$textMatrix, file3_path)
  
  print("check mod_color:")
  print(head(exp.ds$mod_color))
  
  print("check Module")
  print(head(exp.ds$mod_color_anno))
  
  print("moduleTraitCor:")
  print(exp.ds$moduleTraitCor)
  
  print("moduleTraitPvalue:")
  print(exp.ds$moduleTraitPvalue)
  
  shinyjs::show("inter_step2")
  output$inter_mtplot <- renderPlot({
    exp.ds$heatmap = list()
    
    tryCatch({
      col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#1F77B4FF", "white", "#D62728FF"))
      
      exp.ds$heatmap = ComplexHeatmap::Heatmap(
        matrix = exp.ds$moduleTraitCor,
        cluster_rows = T,
        cluster_columns = T,
        left_annotation = exp.ds$Left_anno,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid::grid.text(
            textMatrix[i, j],
            x, y,
            gp = grid::gpar(fontsize = 10)
          )
        },
        row_names_side = "left",
        column_names_rot = 45,
        heatmap_legend_param = list(
          at = seq(-1, 1, 0.5),
          labels = as.character(seq(-1, 1, 0.5)),
          title = "Correlation",
          legend_height = unit(9, "cm"),
          title_position = "lefttop-rot"
        ),
        rect_gp = gpar(col = "black", lwd = 1.2),
        column_title = "Gene-group-trait correlations",
        column_title_gp = gpar(fontsize = 14, fontface = "bold"),
        col = col_fun,
        show_heatmap_legend = TRUE,
        row_split = NULL,
        column_split = colnames(exp.ds$moduleTraitCor),
        row_names_gp = gpar(fontsize = 14),
        column_names_gp = gpar(fontsize = 12)
      )
      
      draw(exp.ds$heatmap, heatmap_legend_side = "right")
      
    }, error = function(e) {
      showNotification(paste("Error in rendering plot:", e$message), type = "error")
    })
  }, width = reactive(input$mtplot_width), height = reactive(input$mtplot_height))
  
  output$KME <- renderReactable({
    df <- round(exp.ds$KME, digits = 2)
    df$Gene <- rownames(df)
    df <- merge(df, exp.ds$Gene2module, by.x = "Gene", by.y = "GID", all.x = TRUE)
    
    csvDownloadButton <- function(tableId, label = "Download CSV", filename = "data.csv") {
      htmltools::tags$button(
        label,
        onclick = sprintf("Reactable.downloadDataCSV('%s', '%s')", tableId, filename)
      )
    }
    hex_to_rgba <- function(hex, alpha = 0.05) {
      rgb <- grDevices::col2rgb(hex) / 255
      sprintf("rgba(%d,%d,%d,%.2f)", rgb[1]*255, rgb[2]*255, rgb[3]*255, alpha)
    }
    
    reactable(
      exp.ds$module_tb, columns = list(
        Module = colDef(name = "Module (Number of Genes)", minWidth = 70),
        module = colDef(show = FALSE),
        Genes =  colDef(name = "Genes", minWidth = 100)
      ), 
      searchable = TRUE, wrap = FALSE,
      details = function(index) {
        plant_data <- df[df$Module == exp.ds$module_tb$module[index], ]
        mm_col <- paste0(exp.ds$module_tb$module[index])
        plant_data <- plant_data %>%
          dplyr::select(Gene, Module, all_of(mm_col)) %>%
          dplyr::rename(`Module Membership (kME)` = all_of(mm_col))
        
        tableId <- paste0("plant_data_", index)
        fileName <- paste0("Gene_Group_", exp.ds$module_tb$module[index], ".csv")
        
        htmltools::div(
          style = "padding: 1rem;",
          htmltools::browsable(
            tagList(
              reactable(
                plant_data, outlined = TRUE, searchable = TRUE, defaultPageSize = 10,
                showPageSizeOptions = TRUE, pageSizeOptions = c(5, 10, 20, 50),
                defaultColDef = colDef(
                  align = "left", minWidth = 150
                ),
                columns = list(
                  Gene = colDef(minWidth = 200)
                ),
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
  
  shiny::showTab(inputId = "inter_tabs", target = "inter_association_tab")
  updateTabsetPanel(session = session, inputId = "inter_tabs", selected = "inter_association_tab")
  shinyjs::show("inter_step5_ui")
  shinyjs::show("inter_step6_ui")
  removeModal()
})

##############################
# 4b Eigengene ####
observeEvent(input$inter_step5,{
  updateTabItems(session, "steps", "hub")
  print("Step 5 button clicked, proceed to 5_Hub Gene Identification.")
})
observeEvent(input$inter_step6,{
  updateTabItems(session, "steps", "function")
  print("Step 6 button clicked, proceed to 6_Functional Enrichment Analysis.")
})
##############################
# download ####
output$inter_download_demo <- downloadHandler(
  filename = function() {
    paste("Demo_GSE15459_trait.csv")
  },
  content = function(file) {
    file.copy("Demo/GSE15459_trait.csv", file)
  }
)

observe({
  req(exp.ds$phen, exp.ds$show_trait_df)
  
  output$inter_download4 <- downloadHandler(
    filename = function() {
      "04_trait.tsv"
    },
    content = function(file) {
      df<- exp.ds$phen 
      df$SampleName <- rownames(df)
      df <- df %>%
        dplyr::select(SampleName, everything())
      readr::write_tsv(x = df, file = file)
    }
  )
})

observeEvent(exp.ds$phen, {
  output$inter_download <- renderUI({
    screenshotButton(
      label = "Download Plot",
      download = TRUE, 
      server_dir = tempdir(), 
      id = "inter_mtplot", 
      filename = "04_GeneGroup2Trait")
  })
  
  output$inter_download2 <- downloadHandler(
    filename = function() {
      paste("04_kMETable.csv")
    },
    content = function(file) {
      df<- round(exp.ds$KME, digits = 2)  
      colnames(df)[1] <- "Gene_ID" 
      write.table(x = df, file = file, sep = ",", quote = FALSE, row.names = TRUE)
    }
  )
})
###################################################################
