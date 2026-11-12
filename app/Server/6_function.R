# server
source("Server/clustergvis_function.R")

# step 6 : Functional Enrichment Analysis
selected_mod <- reactive({
  req(exp.ds$net)
  sorted_mod <- sort(unique(exp.ds$moduleColors))
  sorted_mod <- sorted_mod[sorted_mod != "grey"]
  sorted_mod <- c("All", sorted_mod)
  return(sorted_mod)
})
observe({
  updateSelectInput(session, "func_module_select", choices = selected_mod())
})
observe({
  req(input$func_module_select)
  genelist <- exp.ds$Gene2module[exp.ds$Gene2module$Module %in% input$func_module_select, 1]
  updateSelectInput(session, inputId = "func_genelist", choices = genelist)
})
observeEvent(input$func_go_category, {
  condition_to_show <- input$func_go_category != ""
  shinyjs::toggle(id = "func_go_setting_ui", condition = condition_to_show)
})

observeEvent(input$func_kegg, {
  shinyjs::toggle(id = "func_kegg_setting_ui", condition = input$func_kegg)
})
# david
# observe({
#   req(input$func_go_category)
#   
#   if(input$func_go_db == "david") {
#     category <- input$func_go_category
#     if(exp.ds$organism == "hsa") {
#       base_path <- "david_database/homo"
#       tail <- ".rds"
#     } else if(exp.ds$organism == "tva") {
#       base_path <- "david_database/tva"
#       tail <- "_tva.rds"
#     }
#     term2gene_path <- file.path(base_path, paste0(category, "_TERM2GENE", tail))
#     term2name_path <- file.path(base_path, paste0(category, "_TERM2NAME", tail))
#     print(paste0(term2gene_path))
#     exp.ds$TERM2GENE <- readRDS(term2gene_path)
#     exp.ds$TERM2NAME <- readRDS(term2name_path)
#   } else {
#     exp.ds$TERM2GENE <- NULL
#     exp.ds$TERM2NAME <- NULL
#   }
# })

go_topn = eventReactive(input$func_update,{ 
  as.integer(input$func_go_topn)
})

kegg_topn = eventReactive(input$func_update,{ 
  as.integer(input$func_kegg_topn)
})

# functional enrichment analysis ####
observeEvent(input$func_update, {
  if(is.null(exp.ds$table2)){
    shinyalert("Please preprocess the data", type = "warning")
    updateTabItems(session, "steps", "preprocessing")
    return()
  }
  if(is.null(exp.ds$net)){
    shinyalert("Please preprocess the gene group detection", type = "warning")
    updateTabsetPanel(session = session, inputId = "create_tabs", selected = "module_detection_tab")
    return()
  }
  if(input$func_go_category == "" && input$func_kegg == FALSE) {
    shinyalert("Please select at least one analysis option", type = "warning")
    return()
  }
  if(is.null(input$func_module_select)) {
    shinyalert("Please select at least one gene group", type = "warning")
    return()
  }
  print("Func_functional_analysis_btn is clicked!!")
  
  showModal(modalDialog(
    title = "Loading...",
    easyClose = FALSE,
    footer = NULL,
    div(id = "func_step1", tags$p("Perform enrichment analysis......"))
  ))
  
  if (!is.null(input$func_genelist)){
    print("There has a genelist")
    exp.ds$genelist <- input$func_genelist
  } else {
    print("There has not a genelist")
    exp.ds$genelist <- NULL
  }
  
  #debug
  if (!is.null(exp.ds$genelist)) print("ok") else print("no")
  
  print(paste0("Functional enrichment analysis method: ", "Stringdb"))
  
  stringdb <- TRUE
  string_org <- exp.ds$orgnum
  string_version <- input$func_string_version
  
  if(exp.ds$organism == "tva_g3"){
    string_org <- 5722
  }
  
  if (input$func_go_category != "" && input$func_kegg){
    exp.ds$enrich_mode <- "GO&KEGG"
  } else if (input$func_go_category != "" && input$func_kegg == FALSE) {
    exp.ds$enrich_mode <- "GO"
  } else if (input$func_go_category == "" && input$func_kegg) {
    exp.ds$enrich_mode <- "KEGG"
  }
  print(paste0("enrich mode: ", exp.ds$enrich_mode))
  
  # exp.ds$OrgDb <- exp.ds$orgdb
  print(paste0("Organism: ", exp.ds$organism))

  exp.ds$enrich_anno <- FALSE
  shinyjs::hide("func_heatmap_show_mode_ui")
  shinyjs::hide("func_bar_show_mode_ui")
  
  exp.ds$cw <- ClusterGVis::clusterData(obj = exp.ds$table2, cluster.method = "wgcna", object = exp.ds$net)
  
  if ("All" %in% input$func_module_select) {
    print("Select all modules")
    
    # remove grey module
    grey_cluster <- unique(exp.ds$cw$wide.res[exp.ds$cw$wide.re$modulecol == "grey", "cluster"])
    grey_name <- paste0("C", grey_cluster)

    exp.ds$select_modules_cw <- list(type = "wgcna", geneMode = "none", geneType = "none")
    exp.ds$select_modules_cw$wide.res <- exp.ds$cw$wide.res[exp.ds$cw$wide.re$modulecol != "grey", ]
    colnames(exp.ds$select_modules_cw$wide.res) <-
      sub("^X(?=[0-9])", "", colnames(exp.ds$select_modules_cw$wide.res), perl = TRUE)
    exp.ds$select_modules_cw$long.res <- exp.ds$cw$long.res[exp.ds$cw$long.res$cluster != grey_cluster, ]
    exp.ds$select_modules_cw$cluster.list[[grey_name]] <- NULL
  } else {
    select_modules <- input$func_module_select
    print("Select modules:")
    print(select_modules)
    
    select_modules_cluster <- unique(exp.ds$cw$wide.res[exp.ds$cw$wide.re$modulecol %in% select_modules, "cluster"])
    select_modules_name <- paste0("C", select_modules_cluster)
    
    exp.ds$select_modules_cw <- list(type = "wgcna", geneMode = "none", geneType = "none")
    exp.ds$select_modules_cw$wide.res <- exp.ds$cw$wide.res[exp.ds$cw$wide.re$modulecol %in% select_modules, ]
    exp.ds$select_modules_cw$long.res <- exp.ds$cw$long.res[exp.ds$cw$long.res$cluster %in% select_modules_cluster, ]
    exp.ds$select_modules_cw$cluster.list <- exp.ds$cw$cluster.list[names(exp.ds$cw$cluster.list) %in% select_modules_name]
  }
  
  exp.ds$moduel_anno <- data.frame(
    group = character(), Description = character(), p.adjust = numeric(), FoldEnrichment = numeric(), Count = numeric(), ratio = numeric()
  )
  exp.ds$enrich_go <- data.frame(
    group = character(), Description = character(), p.adjust = numeric(), FoldEnrichment = numeric(), Count = numeric(), ratio = numeric()
  )
  exp.ds$enrich_kegg <- data.frame(
    group = character(), Description = character(), p.adjust = numeric(), FoldEnrichment = numeric(), Count = numeric(), ratio = numeric()
  )
  exp.ds$enrich_bp <- data.frame(
    group = character(), Description = character(), p.adjust = numeric(), FoldEnrichment = numeric(), Count = numeric(), ratio = numeric()
  )
  exp.ds$enrich_mf <- data.frame(
    group = character(), Description = character(), p.adjust = numeric(), FoldEnrichment = numeric(), Count = numeric(), ratio = numeric()
  )
  exp.ds$enrich_cc <- data.frame(
    group = character(), Description = character(), p.adjust = numeric(), FoldEnrichment = numeric(), Count = numeric(), ratio = numeric()
  )
  
  if (exp.ds$enrich_mode == "GO&KEGG"){
    # if(input$func_go_db == "david") {
    #   print("use david database")
    #   exp.ds$enrich_go <- enrichCluster(
    #     object = exp.ds$select_modules_cw, type = "ownSet", id.trans = TRUE, 
    #      
    #     TERM2GENE = exp.ds$TERM2GENE, TERM2NAME = exp.ds$TERM2NAME, 
    #     pvalueCutoff = input$func_go_pvalue, topn = go_topn()
    #   )
    # } else {
    print("use orgdb database")
    
    if(input$func_go_category == "all"){
      exp.ds$enrich_bp <- myenrichCluster(
        object = exp.ds$select_modules_cw, type = "BP", 
        padjCutoff = input$func_go_padj, 
        topn = go_topn(), all_genes = colnames(exp.ds$table2), seed = 2025825,
        stringdb = stringdb, string_org = string_org, string_version = string_version
      )
      exp.ds$enrich_bp <- exp.ds$enrich_bp %>%
        mutate(Description = paste0("BP: ", Description))
      
      exp.ds$enrich_mf <- myenrichCluster(
        object = exp.ds$select_modules_cw, type = "MF", 
        padjCutoff = input$func_go_padj, 
        topn = go_topn(), all_genes = colnames(exp.ds$table2), seed = 2025825,
        stringdb = stringdb, string_org = string_org, string_version = string_version
      )
      exp.ds$enrich_mf <- exp.ds$enrich_mf %>%
        mutate(Description = paste0("MF: ", Description))
      
      exp.ds$enrich_cc <- myenrichCluster(
        object = exp.ds$select_modules_cw, type = "CC", 
        padjCutoff = input$func_go_padj, 
        topn = go_topn(), all_genes = colnames(exp.ds$table2), seed = 2025825,
        stringdb = stringdb, string_org = string_org, string_version = string_version
      )
      exp.ds$enrich_cc <- exp.ds$enrich_cc %>%
        mutate(Description = paste0("CC: ", Description))
      
    } else {
      exp.ds$enrich_go <- myenrichCluster(
        object = exp.ds$select_modules_cw, type = input$func_go_category, 
        padjCutoff = input$func_go_padj, 
        topn = go_topn(), all_genes = colnames(exp.ds$table2), seed = 2025825,
        stringdb = stringdb, string_org = string_org, string_version = string_version
      )
      
      exp.ds$enrich_go <- exp.ds$enrich_go %>%
        mutate(Description = paste0(input$func_go_category, ": ", Description))
    }
    # }
    
    exp.ds$enrich_kegg <- myenrichCluster(
      object = exp.ds$select_modules_cw, type = "KEGG",
      padjCutoff = input$func_kegg_padj, topn = kegg_topn(),
      all_genes = colnames(exp.ds$table2), seed = 2025825,
      stringdb = stringdb, string_org = string_org, string_version = string_version
    )
    exp.ds$enrich_kegg <- exp.ds$enrich_kegg %>%
      mutate(Description = paste0("KEGG: ", Description))
    
    shinyjs::show("func_heatmap_show_mode_ui")
    shinyjs::show("func_bar_show_mode_ui")
  } else if (exp.ds$enrich_mode == "GO") {
    print(paste0("Run go, category is: ", input$func_go_category))
    
    # if(input$func_go_db == "david") {
    #   print("use david database")
    #   exp.ds$enrich_go <- myenrichCluster(
    #     object = exp.ds$select_modules_cw, type = "ownSet", id.trans = TRUE, 
    #      
    #     TERM2GENE = exp.ds$TERM2GENE, TERM2NAME = exp.ds$TERM2NAME, 
    #     pvalueCutoff = input$func_go_padj, topn = go_topn()
    #   )
    # } else {
    print("use orgdb database")
    if(input$func_go_category == "all"){
      exp.ds$enrich_bp <- myenrichCluster(
        object = exp.ds$select_modules_cw, type = "BP", 
        padjCutoff = input$func_go_padj, 
        topn = go_topn(), all_genes = colnames(exp.ds$table2), seed = 2025825,
        stringdb = stringdb, string_org = string_org, string_version = string_version
      )
      exp.ds$enrich_bp <- exp.ds$enrich_bp %>%
        mutate(Description = paste0("BP: ", Description))
      
      exp.ds$enrich_mf <- myenrichCluster(
        object = exp.ds$select_modules_cw, type = "MF", 
        padjCutoff = input$func_go_padj, 
        topn = go_topn(), all_genes = colnames(exp.ds$table2), seed = 2025825,
        stringdb = stringdb, string_org = string_org, string_version = string_version
      )
      exp.ds$enrich_mf <- exp.ds$enrich_mf %>%
        mutate(Description = paste0("MF: ", Description))
      
      exp.ds$enrich_cc <- myenrichCluster(
        object = exp.ds$select_modules_cw, type = "CC", 
        padjCutoff = input$func_go_padj, 
        topn = go_topn(), all_genes = colnames(exp.ds$table2), seed = 2025825,
        stringdb = stringdb, string_org = string_org, string_version = string_version
      )
      exp.ds$enrich_cc <- exp.ds$enrich_cc %>%
        mutate(Description = paste0("CC: ", Description))
      
    } else {
      exp.ds$enrich_go <- myenrichCluster(
        object = exp.ds$select_modules_cw, type = input$func_go_category, 
        padjCutoff = input$func_go_padj, 
        topn = go_topn(), all_genes = colnames(exp.ds$table2), seed = 2025825,
        stringdb = stringdb, string_org = string_org, string_version = string_version
      )
      
      exp.ds$enrich_go <- exp.ds$enrich_go %>%
        mutate(Description = paste0(input$func_go_category, ": ", Description))
    }
    # }
    
    
  } else if (exp.ds$enrich_mode == "KEGG") {
    exp.ds$enrich_kegg <- myenrichCluster(
      object = exp.ds$select_modules_cw, type = "KEGG",
      padjCutoff = input$func_kegg_padj, topn = kegg_topn(),
      all_genes = colnames(exp.ds$table2), seed = 2025825,
      stringdb = stringdb, string_org = string_org, string_version = string_version
    )
    exp.ds$enrich_kegg <- exp.ds$enrich_kegg %>%
      mutate(Description = paste0("KEGG: ", Description))
  }
  
  shinyjs::show("creat_step2")
  
  exp.ds$cluster_info <- exp.ds$select_modules_cw$wide.res %>%
    dplyr::select(cluster, modulecol) %>%
    distinct() %>%
    arrange(cluster) %>%
    mutate(group = paste0("C", cluster))
  
  exp.ds$original_cluster_colors <- exp.ds$cluster_info$modulecol
  
  print(exp.ds$original_cluster_colors)
  
  output$func_choose_ui <- renderUI({
    radioGroupButtons(
      inputId = "func_choose",
      choices = c(`<i class='fa fa-rectangle-list'></i>` = "heatmap",
                  `<i class='fa fa-chart-bar'></i>` = "bar_plot"),
      selected = "heatmap", justified = TRUE, width = "50", status = "info"
    )
  })
  exp.ds$enrich_anno <- TRUE

  shinyjs::show("func_download_ui")
  shinyjs::show("func_download2_ui")
  removeModal()
})

# modules2changecolor <- reactive({
#   if ("All" %in% input$func_module_select) {
#     modules <- unique(exp.ds$moduleColors)
#     modules <- modules[modules != "grey"]
#   } else {
#     modules <- input$func_module_select
#   }
#   modules
# })
# observe({
#   updateSelectInput(session, "func_heatmap_module_color", choices = modules2changecolor())
# })
# modules2changeblack <- reactive({
#   if (is.null(input$func_heatmap_module_color) || length(input$func_heatmap_module_color) == 0) {
#     return(character(0))
#   }
  
#   if ("All" %in% input$func_heatmap_module_color) {
#     modules <- unique(exp.ds$moduleColors)
#     modules[modules != "grey"]
#     print("All modules present with black")
#   } else {
#     input$func_heatmap_module_color
#     print(input$func_heatmap_module_color)
#     print("present with black")
#   }
# })

box_style <- reactive({
  if (input$func_heatmap_box_style == "link") {
    "anno_link"
  } else {
    "anno_block"
  }
})
anno_size <- reactive({
  as.numeric(input$func_heatmap_anno_size)
})
# plot_type <- reactive({if(input$func_heatmap_plot_type == "heatmap") {"heatmap"} else {"both"}})
x_expand <- reactive({
  as.numeric(input$func_bar_x_expand)
})
facet_ncol <- reactive({
  as.numeric(input$func_bar_facet_ncol)
})
text_size <- reactive({
  as.numeric(input$func_bar_text_size)
})

generateFuncHeatmap <- function(anno_data, go_col, kegg_data = NULL, kegg_col = NULL) {
  go_col <- go_col
  kegg_col <- kegg_col
  
  if(input$func_heatmap_module_color){
    go_col <- "black" 
    if (!is.null(kegg_data)) {
      kegg_col <- "black" 
    }
  }

  myvisCluster(
    object = exp.ds$select_modules_cw,
    # plot.type = "both",
    column_names_rot = 45,
    show_row_dend = FALSE,
    annoTerm.data = anno_data,
    go.col = go_col,
    annoKegg.data = kegg_data,
    kegg.col = kegg_col,
    ctAnno.col = exp.ds$original_cluster_colors,
    line.side = "left",
    word_wrap = FALSE,
    add_new_line = TRUE,
    by.go = box_style(),
    by.kegg = if(!is.null(kegg_data)) box_style() else NULL,
    go.size = anno_size(),
    kegg.size = if(!is.null(kegg_data)) anno_size() else NULL,
    markGenes = if (!is.null(exp.ds$genelist)) exp.ds$genelist else NULL
  )
}

generateFuncBar <- function(df, separate = TRUE) {
  df <- df %>%
    mutate(color = exp.ds$group_to_color_map[group],
           ratio = round(ratio, 2),
           Description = paste0(Description, " (", ratio, ")"),
           number = as.numeric(sub("C", "", group)),
           type = ifelse(separate & grepl("KEGG", Description), "KEGG", "GO"),
           group_type = if(separate) interaction(group, type, sep = "_") else NULL)
  
  df <- df[order(df$number), ]
  
  if(separate) {
    color_labeller <- as.list(setNames(df$color, df$group_type))
    
    p <- ggplot(df, aes(x = ratio, y = reorder(Description, ratio), fill = color)) +
      geom_col(alpha = 0.6, width = 0.7) +
      scale_fill_identity() +
      geom_text(aes(label = Description), position = position_stack(vjust = 0), 
                hjust = -0.01, color = "black", fontface = "bold", size = text_size()) +
      labs(x = "Ratio", y = "Pathway") +
      theme_classic() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", size = 12, face = "bold"),
        legend.position = "none",
        panel.spacing = unit(1, "lines")
      ) +
      facet_wrap(~group_type, scales = "free", ncol = facet_ncol(),
                 labeller = labeller(group_type = function(x) color_labeller[x])) +
      scale_x_continuous(expand = expansion(mult = c(0, x_expand())))
  } else {
    p <- ggplot(df, aes(x = ratio, y = reorder(Description, ratio), fill = color)) +
      geom_col(alpha = 0.6, width = 0.7) +
      scale_fill_identity() +
      geom_text(aes(label = Description), position = position_stack(vjust = 0),
                hjust = -0.01, color = "black", fontface = "bold", size = text_size()) +
      labs(x = "Ratio", y = "Pathway") +
      theme_classic() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", size = 12, face = "bold"),
        legend.position = "none",
        panel.spacing = unit(1, "lines")
      ) +
      facet_wrap(~color, scales = "free", ncol = facet_ncol()) +
      scale_x_continuous(expand = expansion(mult = c(0, x_expand())))
  }
  
  return(p)
}

output$func_bar_plot_ui <- renderUI({
  req(input$func_bar_width, input$func_bar_height)
  plotOutput("func_bar_plot", width = input$func_bar_width, height = input$func_bar_height)
})

output$func_heatmap_ui <- renderUI({
  req(input$func_heatmap_width, input$func_heatmap_height)
  plotOutput("func_heatmap", width = input$func_heatmap_width, height = input$func_heatmap_height) 
})

# plot ####
observe({
  req(exp.ds$enrich_anno)
  
  if ((exp.ds$enrich_mode == "GO") ||
      (exp.ds$enrich_mode == "KEGG")) {
    
    if(input$func_go_category == "all"){
      exp.ds$moduel_anno <- 
        rbind(exp.ds$enrich_bp %>% group_by(group) %>% slice_head(n = go_topn()),
              exp.ds$enrich_mf %>% group_by(group) %>% slice_head(n = go_topn()),
              exp.ds$enrich_cc %>% group_by(group) %>% slice_head(n = go_topn()),
              exp.ds$enrich_kegg %>% group_by(group) %>% slice_head(n = kegg_topn()))
    } else {
      exp.ds$moduel_anno <- 
        rbind(exp.ds$enrich_go %>% group_by(group) %>% slice_head(n = go_topn()),
              exp.ds$enrich_kegg %>% group_by(group) %>% slice_head(n = kegg_topn())) 
    }
    exp.ds$group_to_color_map <- setNames(exp.ds$cluster_info$modulecol, exp.ds$cluster_info$group)
    exp.ds$final_color <- exp.ds$group_to_color_map[exp.ds$moduel_anno$group]
    
    output$func_heatmap <- renderPlot({
      generateFuncHeatmap(exp.ds$moduel_anno, exp.ds$final_color)
    }, width = function() input$func_heatmap_width, height = function() input$func_heatmap_height,
    execOnResize = TRUE)
    
  } else if (exp.ds$enrich_mode == "GO&KEGG") {
    if(input$func_heatmap_show_mode == "separate") {
      if(input$func_go_category == "all"){
      exp.ds$moduel_anno <- 
        rbind(exp.ds$enrich_bp %>% group_by(group) %>% slice_head(n = go_topn()),
              exp.ds$enrich_mf %>% group_by(group) %>% slice_head(n = go_topn()),
              exp.ds$enrich_cc %>% group_by(group) %>% slice_head(n = go_topn()))
      } else {
        exp.ds$moduel_anno <- exp.ds$enrich_go
      }
      exp.ds$group_to_color_map <- setNames(exp.ds$cluster_info$modulecol, exp.ds$cluster_info$group)
      exp.ds$go_final_color <- exp.ds$group_to_color_map[exp.ds$moduel_anno$group]
      exp.ds$kegg_final_color <- exp.ds$group_to_color_map[exp.ds$enrich_kegg$group]
      
      output$func_heatmap <- renderPlot({
        generateFuncHeatmap(
          exp.ds$moduel_anno, exp.ds$go_final_color, exp.ds$enrich_kegg, exp.ds$kegg_final_color
        )
      }, width = function() input$func_heatmap_width, height = function() input$func_heatmap_height,
      execOnResize = TRUE)
    } else {
      if(input$func_go_category == "all"){
        exp.ds$moduel_anno <- 
          rbind(exp.ds$enrich_bp %>% group_by(group) %>% slice_head(n = go_topn()),
                exp.ds$enrich_mf %>% group_by(group) %>% slice_head(n = go_topn()),
                exp.ds$enrich_cc %>% group_by(group) %>% slice_head(n = go_topn()),
                exp.ds$enrich_kegg %>% group_by(group) %>% slice_head(n = kegg_topn()))
      } else {
        exp.ds$moduel_anno <- 
          rbind(exp.ds$enrich_go %>% group_by(group) %>% slice_head(n = go_topn()),
                exp.ds$enrich_kegg %>% group_by(group) %>% slice_head(n = kegg_topn())) 
      }
      
      exp.ds$group_to_color_map <- setNames(exp.ds$cluster_info$modulecol, exp.ds$cluster_info$group)
      exp.ds$final_color <- exp.ds$group_to_color_map[exp.ds$moduel_anno$group]
      
      output$func_heatmap <- renderPlot({
        generateFuncHeatmap(exp.ds$moduel_anno, exp.ds$final_color)
      }, width = function() input$func_heatmap_width, height = function() input$func_heatmap_height,
      execOnResize = TRUE)
    }
  } 
  
  output$func_bar_plot <- renderPlot({
    if(exp.ds$enrich_mode == "GO&KEGG" && input$func_heatmap_show_mode == "separate") {
      df <- rbind(exp.ds$moduel_anno, exp.ds$enrich_kegg)
    } else {
      df <- exp.ds$moduel_anno
    }
    separate <- (input$func_bar_show_mode == "separate")
    
    generateFuncBar(df, separate = separate)
  }, width = input$func_bar_width, input$func_bar_height)
  
})

##############################
# download
observeEvent(input$func_update, {
  output$func_download <- downloadHandler(
    filename = function() { paste0("05_FuncHeatmap.pdf") },
    content = function(file) {
      plot_width_in <- input$func_heatmap_width / 96
      plot_height_in <- input$func_heatmap_height / 96
      
      pdf(file, height = plot_height_in, width = plot_width_in, onefile = F)
      if((exp.ds$enrich_mode == "GO") ||
         (exp.ds$enrich_mode == "KEGG")) {
        generateFuncHeatmap(exp.ds$moduel_anno, exp.ds$final_color)
      } else if(exp.ds$enrich_mode == "GO&KEGG") {
        if(input$func_heatmap_show_mode == "separate") {
          generateFuncHeatmap(exp.ds$enrich_go, exp.ds$go_final_color,
                              exp.ds$enrich_kegg, exp.ds$kegg_final_color)
        } else {
          generateFuncHeatmap(exp.ds$moduel_anno, exp.ds$final_color)
        }
      }
      dev.off()
    }
  )
  
  output$func_download2 <- downloadHandler(
    filename = function() { paste0("05_Bar_plot.png") },
    content = function(file) {
      png(file,
          width = input$func_bar_width,
          height = input$func_bar_height,
          units = "px", res = 96)
      
      df <- exp.ds$moduel_anno
      separate <- (input$func_bar_show_mode == "separate")
      print(generateFuncBar(df, separate = separate))
      dev.off()
    }
  )
})






