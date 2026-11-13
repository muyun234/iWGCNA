# server
hub.ds <- reactiveValues(
  gs_mm_module = NULL,
  gs_mm_trait = NULL,
  hub_module = NULL,
  hub_trait = NULL,
  hub_gene_info = NULL
)

generate_gsmm_df <- function(exp.ds, trait, module){
  
  selected_df <- as.data.frame(exp.ds$phen[, colnames(exp.ds$phen) == trait])
  names(selected_df) <- trait
  
  modNames = substring(names(exp.ds$MEs_col), 3)
  geneTraitSignificance <- as.data.frame(stats::cor(exp.ds$table2, selected_df, use = "p"))
  GSPvalue <- as.data.frame(WGCNA::corPvalueStudent(as.matrix(geneTraitSignificance), exp.ds$nSamples))
  names(geneTraitSignificance) <- names(selected_df)
  names(GSPvalue) <- names(selected_df)
  
  column <- match(module, modNames)
  if (is.na(column)) {
    stop("Selected module '", module, "' was not found. Available modules: ", paste(modNames, collapse = ", "))
  }
  moduleGenes_name <- exp.ds$Gene2module$GID[exp.ds$Gene2module$Module == module]
  
  moduleGenes_MM <- exp.ds$KME[moduleGenes_name, column, drop = FALSE]
  moduleGenes_GS <- geneTraitSignificance[moduleGenes_name, 1, drop = FALSE]
  
  genes <- rownames(moduleGenes_MM)
  idx <- match(genes, rownames(moduleGenes_GS))
  
  hubgene_df <- cbind(
    moduleGenes_MM,
    moduleGenes_GS[idx, , drop = FALSE]
  )
  rownames(hubgene_df) <- genes
  colnames(hubgene_df) <- c("MM", "GS")
  
  # hubgene_df[] <- abs(hubgene_df)
  hubgene_df$absMM <- abs(hubgene_df$MM)
  hubgene_df$absGS <- abs(hubgene_df$GS)
  
  return(hubgene_df)
}

generate_ppi_network <- function(exp.ds, hub_gene_info, hub_ppi_level) {
  gene_df <- data.frame(gene = hub_gene_info$Gene)

  map_path <- paste0("stringdb/", exp.ds$orgnum, "_map.RDS")
  print(paste("map path: ", map_path))
  
  alias <- readRDS(map_path)
  
  mapped_all <- gene_df |>
    dplyr::inner_join(alias, by = "gene") |>
    dplyr::distinct(gene, STRING_id)
  
  ambiguous_genes <- mapped_all |>
    dplyr::count(gene, name = "n_STRING_ids") |>
    dplyr::filter(n_STRING_ids > 1)
  
  offline_mapped_genes <- mapped_all |>
    dplyr::anti_join(ambiguous_genes, by = "gene")
  
  unmapped_genes <- gene_df |>
    dplyr::anti_join(alias, by = "gene")
  
  print(paste("Input genes:", nrow(gene_df)))
  print(paste("Uniquely mapped genes:", nrow(offline_mapped_genes)))
  print(paste("Ambiguous genes:", nrow(ambiguous_genes)))
  print(paste("Unmapped genes:", nrow(unmapped_genes)))
  
  ppi_path <- paste0("stringdb/", exp.ds$orgnum , "_", hub_ppi_level, "_ppi.RDS")
  print(paste("map path: ", ppi_path))
  ppi <- readRDS(ppi_path)
  
  hits <- offline_mapped_genes$STRING_id
  
  ppi_sub <- ppi[
    ppi$protein1 %in% hits &
      ppi$protein2 %in% hits,
  ]
  colnames(ppi_sub) <- c("from", "to", "combined_score")

  edges <- as.data.frame(ppi_sub) %>%
    dplyr::transmute(from = as.character(from), to = as.character(to)) %>%
    dplyr::filter(!is.na(from), !is.na(to), from != to) %>%
    dplyr::transmute(protein1 = pmin(from, to), protein2 = pmax(from, to)) %>%
    dplyr::distinct(protein1, protein2)
  
  degree_counts <- data.frame(STRING_id = c(edges$protein1, edges$protein2)) %>%
    dplyr::count(STRING_id, name = "degree")
  
  protein_degree <- data.frame(
    STRING_id = unique(as.character(offline_mapped_genes$STRING_id))
  ) %>%
    dplyr::left_join(degree_counts, by = "STRING_id") %>%
    dplyr::mutate(degree = dplyr::coalesce(degree, 0L))
  
  node_info <- offline_mapped_genes %>%
    dplyr::left_join(protein_degree, by = "STRING_id")
  
  return(list(
    node_info = node_info
  ))
}

#########################
# select module and trait
s_mod <- reactive({
  sorted_mod <- rownames(exp.ds$moduleTraitPvalue)
  sort(sorted_mod)
})

s_trait <- reactive({
  sort(colnames(exp.ds$moduleTraitPvalue))
})

observe({
  updateSelectInput(session, "hub_module", choices = s_mod())
  updateSelectInput(session, "hub_trait", choices = s_trait())
})

observe({
  req(input$demo_on, exp.ds$moduleTraitPvalue)
  if(input$demo_on){
    updateSelectInput(session, "hub_module", selected = "salmon")
    updateSelectInput(session, "hub_trait", choices = "Stage")
    updateSliderInput(session, "hub_MMcut", value = 0.8)
    updateSliderInput(session, "hub_GScut", value = 0.2) 
  }
})

observeEvent(
  list(input$hub_module, input$hub_trait),
  {
    # 清除物件
    hub.ds$hubgene_df <- NULL
    hub.ds$hubgene_cor <- NULL
    hub.ds$hubgene_corp <- NULL
    hub.ds$hubgene_plot <- NULL
    hub.ds$hubgene_filter <- NULL
    hub.ds$wgcna_hubgene_df <- NULL
    hub.ds$hub_gene_info <- NULL
    hub.ds$modTOM <- NULL
    hub.ds$ppi_hub_gene_info <- NULL
    hub.ds$ppi_list <- NULL
    hub.ds$wgcna_edge_df <- NULL
    hub.ds$wgcna_node_df <- NULL
    hub.ds$wgcna_g <- NULL
    hub.ds$connect_ppi_df <- NULL
    
    # 重設畫圖狀態
    hub.ds$select_hubgene <- FALSE
    hub.ds$plot_connect_mm <- FALSE
    hub.ds$plot_wgcna_network <- FALSE
    hub.ds$plot_connect_ppi <- FALSE
    
    # UI
    shinyjs::hide("hub_gs2mm_plot_all_ui")
    shinyjs::hide("hub_gs2mmdt_ui")
    shinyjs::hide("hub_ppi_update_ui")
    shinyjs::hide("hub_connect_mm_plot_all_ui")
    shinyjs::hide("hub_connect_mm_dt_ui")
    shinyjs::hide("hub_wgcna_plot_ui")
    shinyjs::hide("hub_connect_ppi_plot_ui")
    shinyjs::hide("hub_connect_ppi_dt_ui")
    shinyjs::show("hub_gsmm_run_ui")
    shinyjs::show("hub_run_ui")
  },
  ignoreInit = TRUE,
  priority = 100
)

#########################
# GS MM plot ####
observeEvent(input$hub_gsmm_run, {
  if (is.null(exp.ds$phen)) {
    shinyalert("Error", "Please upload a trait file", type = "error")
    updateTabItems(session, "steps", "interested")
    return()
  }
  if(input$hub_module == "") {
    shinyalert("Please select a gene group", type = "warning")
    return()
  }
  if(input$hub_trait == "") {
    shinyalert("Please select a trait", type = "warning")
    return()
  }
  
  showModal(modalDialog(
    title = "Preparing hub gene data...",
    easyClose = FALSE,
    footer = NULL
  ))
  shinyjs::show("hub_choose_ui")
  # updateRadioGroupButtons(session, "hub_choose", selected = "mm_gs")
  
  hub.ds$select_hubgene <- FALSE
  
  if (is.null(exp.ds$KME)) {
    print("calculate KME")
    exp.ds$KME <- WGCNA::signedKME(exp.ds$table2, exp.ds$MEs_col, outputColumnName="")
  }

  hub.ds$gs_mm_module <- as.character(input$hub_module)
  hub.ds$gs_mm_trait <- as.character(input$hub_trait)
  
  # DEBUG
  print(paste("Selected module:", hub.ds$gs_mm_module))
  print(paste("Selected trait:", hub.ds$gs_mm_trait))
  
  hubgene_df <- generate_gsmm_df(exp.ds, hub.ds$gs_mm_trait, hub.ds$gs_mm_module)
  hub.ds$hubgene_df <- hubgene_df
  
  corExpr <- parse(text = paste("stats::cor", "(hub.ds$hubgene_df$absMM, hub.ds$hubgene_df$absGS ", 
                                WGCNA::prepComma("use = 'p'"), ")"))
  hubgene_cor <- signif(eval(corExpr), 2)

  corp <- WGCNA::corPvalueStudent(
    hubgene_cor, sum(is.finite(hub.ds$hubgene_df$absMM) & is.finite(hub.ds$hubgene_df$absGS))
  )
  
  hub.ds$hubgene_cor <- hubgene_cor
  hub.ds$hubgene_corp <- corp
  
  hub.ds$select_hubgene <- TRUE
  shinyjs::show("hub_gs2mm_plot_all_ui")
  shinyjs::show("hub_gs2mmdt_ui")
  shinyjs::hide("hub_gsmm_run_ui")
  
  removeModal()
})

observe({
  req(input$hub_MMcut, input$hub_GScut, hub.ds$select_hubgene)
  
  hub.ds$MMcut <- as.numeric(input$hub_MMcut)
  hub.ds$GScut <- as.numeric(input$hub_GScut)
  
  hub.ds$hubgene_filter <- hub.ds$hubgene_df %>%
    dplyr::filter(absMM >= hub.ds$MMcut & absGS >= hub.ds$GScut)
  
  hubgene_plot <- hub.ds$hubgene_df %>%
    mutate(color = ifelse(absMM >= hub.ds$MMcut & absGS >= hub.ds$GScut, hub.ds$gs_mm_module, "grey"))
  hubgene_plot$Gene <- rownames(hubgene_plot)
  hub.ds$hubgene_plot <- hubgene_plot
  # print(head(hub.ds$hubgene_plot))
  # hub.ds$hubgene_plot$Gene <- rownames(hub.ds$hubgene_plot)
})

output$hub_gsmm_text <- renderUI({
  req(hub.ds$MMcut, hub.ds$GScut, hub.ds$select_hubgene)
  tags$div(
    style = "font-size: 18px; font-weight: bold;padding-top: 10px;",
    "Highlighting genes with MM ≥ ",
    tags$span(hub.ds$MMcut, style = "color: red;"),
    " and GS ≥ ",
    tags$span(hub.ds$GScut, style = "color: red;")
  )
})

output$hub_gs2mm_plot_ui <- renderUI({
  req(input$hub_gsmmplot_width, input$hub_gsmmplot_height)
  plotlyOutput("hub_gs2mmplot", width = input$hub_gsmmplot_width, height = input$hub_gsmmplot_height) 
})

output$hub_gs2mmplot <- renderPlotly({
  req(hub.ds$hubgene_plot, hub.ds$select_hubgene)
  p <- ggplot(hub.ds$hubgene_plot, aes(x = absMM, y = absGS)) +
    geom_point(aes(fill = color, text = Gene), shape = 21, color = "black", size = 2) +
    geom_smooth(method = 'lm', color = "royalblue") +
    xlab(paste0("Absolute Module membership (|MM|) in ", hub.ds$gs_mm_module, " module")) +
    ylab(paste0("Absolute Gene significance (|GS|) for ", hub.ds$gs_mm_trait)) +
    ggtitle(paste0("Absolute module membership vs. gene significance\n", "cor = ", 
                   signif(hub.ds$hubgene_cor, 2), ", p = ", signif(hub.ds$hubgene_corp, 2))) +
    theme_prism(border = T) +
    theme(legend.position = "none") +
    scale_fill_identity() 
  
  ggplotly(p, width = input$hub_gsmmplot_width, height = input$hub_gsmmplot_height) %>%
    plotly::layout(
      xaxis = list(showline = TRUE, linewidth = 2, linecolor = "black", mirror = TRUE),
      yaxis = list(showline = TRUE, linewidth = 2, linecolor = "black", mirror = TRUE),
      plot_bgcolor = "white", paper_bgcolor = "white", margin = list(t = 100)
    )
})

output$hub_kMEhub <- renderDT({
  req(hub.ds$hubgene_filter)
  df <- hub.ds$hubgene_filter
  df$Gene <- rownames(df)
  df <- df %>%
    dplyr::select(Gene, MM, GS) %>%
    dplyr::rename(
      "Module Membership (MM)" = MM,
      "Gene Significance (GS)" = GS
    )
  
  numeric_cols <- sapply(df, is.numeric)
  df[numeric_cols] <- round(df[numeric_cols], 2)
  
  datatable(df, class = 'compact norap', selection = "none", rownames = FALSE, 
            #caption = tags$caption("content", style = "caption-side: bottom; color:grey; font-size:10; text-align:left;"),
            filter = list(position = 'top', clear = FALSE, plain = TRUE), 
            options = list(
              pageLength = 5, lengthMenu = c(5, 10, 25), scrollX = TRUE, search = list(regex = TRUE),
              order = list(list(1, "desc"), list(2, "desc")),
              columnDefs = list(
                list(targets = "_all", className = "dt-left"),
                list(targets = c("Module Membership (MM)", "Gene Significance (GS)"), searchable = FALSE)
              )
            ) 
  )
})

observeEvent(input$hub_step6,{
  updateTabItems(session, "steps", "function")
  print("Step 6 button clicked, proceed to 6_Functional Enrichment Analysis.")
})
#####################
# hub network ####
observeEvent(input$hub_run, {
  if (is.null(exp.ds$phen)) {
    shinyalert("Error", "Please upload a trait file", type = "error")
    updateTabItems(session, "steps", "interested")
    return()
  }
  if(input$hub_module == "") {
    shinyalert("Please select a gene group", type = "warning")
    return()
  }
  if(input$hub_trait == "") {
    shinyalert("Please select a trait", type = "warning")
    return()
  }
  if(is.na(exp.ds$orgnum)) {
    shinyalert("Warning", "The STRINGdb database for this organism has not been created yet.", type = "warning")
    return()
  }
  
  showModal(modalDialog(
    title = "Preparing hub gene data...",
    easyClose = FALSE,
    footer = NULL
  ))
  # updateRadioGroupButtons(session, "hub_choose", selected = "hub")
  
  hub.ds$plot_connect_mm <- FALSE
  hub.ds$plot_wgcna_network <- FALSE
  hub.ds$plot_connect_ppi <- FALSE
  shinyjs::hide("hub_ppi_update_ui")

  hub.ds$hub_module <- as.character(input$hub_module)
  hub.ds$hub_trait <- as.character(input$hub_trait)
  
  hubgene_df <- generate_gsmm_df(exp.ds, hub.ds$hub_trait, hub.ds$hub_module)
  hubgene_df$gene <- rownames(hubgene_df)
  hub.ds$wgcna_hubgene_df <- hubgene_df
  
  ## connectivity v.s. MM ####
  # connectivity <- rowSums(modTOM)
  module_connectivity <- exp.ds$connectivity_df |>
    dplyr::filter(Module == hub.ds$hub_module) |>
    dplyr::select(Gene, Connectivity = kWithin)
  
  hub_gene_info <- hub.ds$wgcna_hubgene_df %>%
    dplyr::transmute(Gene = gene, AbsMM = absMM, AbsGS = absGS) %>%
    dplyr::inner_join(module_connectivity, by = "Gene") |>
    dplyr::select(Gene, Connectivity, AbsMM, AbsGS) |>
    dplyr::arrange(dplyr::desc(Connectivity))
  
  hub.ds$hub_gene_info <- hub_gene_info
  hub.ds$plot_connect_mm <- TRUE
  
  ## WGCNA network ####
  probes = names(exp.ds$table2)
  inModule = is.finite(match(exp.ds$moduleColors, hub.ds$hub_module))
  modProbes = probes[inModule]
  
  modTOM = exp.ds$TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  
  hub.ds$modTOM <- modTOM
  
  ## PPI network ####
  shinyjs::show("hub_step1")
  
  ppi_hub_gene_info <- hub.ds$hub_gene_info
  
  hub.ds$ppi_hub_gene_info <- ppi_hub_gene_info
  # hub.ds$ppi_list <- generate_ppi_network(exp.ds, hub.ds$ppi_hub_gene_info, input$hub_ppi_level)
  
  ## mm vs ppi degree ####
  
  shinyjs::show("hub_connect_mm_plot_all_ui")
  shinyjs::show("hub_connect_mm_dt_ui")
  shinyjs::show("hub_wgcna_plot_ui")
  shinyjs::show("hub_ppi_update_ui")
  shinyjs::show("hub_connect_ppi_plot_ui")
  shinyjs::show("hub_connect_ppi_dt_ui")
  shinyjs::hide("hub_run_ui")
  
  removeModal()
})

con_mm_cutoff <- reactive({
  req(hub.ds$hub_gene_info, input$hub_top_cut)
  
  top_n <- as.numeric(input$hub_top_cut)
  
  req(is.finite(top_n), top_n >= 0, top_n <= 100)
  
  quantile(
    hub.ds$hub_gene_info$Connectivity,
    probs = 1 - top_n / 100,
    na.rm = TRUE
  )
})

# Connectivity–MM
connect_mm_data <- reactive({
  req(hub.ds$hub_gene_info)
  
  df <- hub.ds$hub_gene_info
  top_n <- as.numeric(input$hub_top_cut)
  top_lab <- paste0("Top ", top_n, "% by connectivity")
  
  df$Group <- factor(
    ifelse(
      df$Connectivity > con_mm_cutoff(),
      top_lab,
      "Other genes"
    ),
    levels = c(top_lab, "Other genes")
  )
  
  df
})

# wgcna network
observe({
  req(hub.ds$hub_gene_info, hub.ds$modTOM, hub.ds$wgcna_hubgene_df, input$hub_wgcna_hub_num)
  
  showModal(modalDialog(
    title = "Building WGCNA network...",
    easyClose = FALSE,
    footer = NULL
  ))
  
  hub.ds$plot_wgcna_network <- FALSE
  hub.ds$wgcna_hub_num <- as.integer(input$hub_wgcna_hub_num)
  
  genes_to_display <- hub.ds$hub_gene_info$Gene[
    hub.ds$hub_gene_info$Connectivity > con_mm_cutoff()
  ]
  genes_to_label <- head(genes_to_display, hub.ds$wgcna_hub_num)
  
  print("WGCNA network node: ")
  print(genes_to_display)
  
  selected_TOM <- hub.ds$modTOM[genes_to_display, genes_to_display, drop = FALSE]
  selected_node_attr <- hub.ds$wgcna_hubgene_df[genes_to_display, , drop = FALSE]
  
  cyt = WGCNA::exportNetworkToCytoscape(selected_TOM, weighted = TRUE, threshold = 0, nodeAttr = selected_node_attr)
  
  hub.ds$wgcna_edge_df <- cyt[[1]]
  hub.ds$wgcna_node_df <- cyt[[2]]
  
  g <- igraph::graph_from_data_frame(
    d = hub.ds$wgcna_edge_df[, c("fromNode", "toNode", "weight")],
    directed = FALSE, vertices = data.frame(name = genes_to_display)
  )
  
  node_info <- data.frame(
    gene = igraph::V(g)$name,
    degree = as.numeric(igraph::degree(g))
  )
  
  wgcna_g <- tidygraph::as_tbl_graph(g) %>%
    tidygraph::activate(nodes) %>%
    dplyr::left_join(node_info, by = c("name" = "gene")) %>%
    dplyr::left_join(selected_node_attr, by = c("name" = "gene")) %>%
    dplyr::left_join(hub.ds$hub_gene_info %>% dplyr::select(Gene, Connectivity), by = c("name" = "Gene")) %>%
    dplyr::mutate(label = ifelse(name %in% genes_to_label, name, ""), MM = absMM, GS = absGS)
  
  hub.ds$wgcna_g <- wgcna_g
  
  hub.ds$plot_wgcna_network <- TRUE
  
  removeModal()
})

# mm vs ppi degree
ppi_network <- reactive({
  req(hub.ds$ppi_hub_gene_info, input$hub_ppi_level)
  
  showModal(modalDialog(
    title = "Building PPI network...",
    easyClose = FALSE,
    footer = NULL
  ))
  
  on.exit(removeModal(), add = TRUE)
  
  result <- generate_ppi_network(exp.ds, hub.ds$ppi_hub_gene_info, input$hub_ppi_level)
  
  result
})

observe({
  req(hub.ds$ppi_hub_gene_info, input$hub_ppi_degree_percent)
  
  hub.ds$plot_connect_ppi <- FALSE
  
  ppi <- ppi_network()
  node_info <- ppi$node_info
  
  df <- hub.ds$ppi_hub_gene_info %>%
    dplyr::left_join(node_info, by = c("Gene" = "gene"))
  df$degree[is.na(df$degree)] <- 0
  
  top_connect <- df$Gene[df$Connectivity > con_mm_cutoff()]
  
  n_top <- ceiling(nrow(df) * as.numeric(input$hub_ppi_degree_percent) / 100)
  top_ppi <- df %>%
    dplyr::filter(degree > 0) %>%
    dplyr::arrange(dplyr::desc(degree), Gene) %>%
    dplyr::slice_head(n = n_top) %>%
    dplyr::pull(Gene)
  
  top_both <- base::intersect(top_connect, top_ppi)
  
  df$type <- "Other"
  df$type[df$Gene %in% top_connect] <- "Top Connectivity"
  df$type[df$Gene %in% top_ppi] <- "Top Degree"
  df$type[df$Gene %in% top_both] <- "Top Connectivity & Degree"
  
  hub.ds$ppi_list <- ppi
  hub.ds$connect_ppi_df <- df
  hub.ds$plot_connect_ppi <- TRUE
})

## plot ####
output$hub_connect_mm_plot_ui <- renderUI({
  req(input$hub_connect_mm_plot_width, input$hub_connect_mm_plot_height)
  
  plotlyOutput(
    "hub_connect_mm_plot",
    width = paste0(input$hub_connect_mm_plot_width, "px"),
    height = paste0(input$hub_connect_mm_plot_height, "px")
  )
})

hub_connect_mm_ggplot <- reactive({
  df <- connect_mm_data()
  
  ggplot(
    df,
    aes(
      x = AbsMM, y = Connectivity, color = Group,
      text = paste0(
        "Gene: ", Gene,
        "<br>MM: ", signif(AbsMM, 2),
        "<br>GS: ", signif(AbsGS, 2),
        "<br>Intramodular Connectivity: ", signif(Connectivity, 2)
      )
    )
  ) +
    geom_point(size = 2) +
    scale_color_manual(
      values = setNames(c("red", "black"), levels(df$Group)), drop = FALSE
    ) +
    labs(
      title = paste0(
        "Hub Gene Identification - ",
        hub.ds$hub_module, " Module"
      ),
      x = "Module Membership (MM)",
      y = "Intramodular Connectivity (kWithin)",
      color = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
})

output$hub_connect_mm_plot <- renderPlotly({
  req(hub.ds$plot_connect_mm)
  
  ggplotly(
    hub_connect_mm_ggplot(),
    tooltip = "text",
    width = input$hub_connect_mm_plot_width,
    height = input$hub_connect_mm_plot_height
  ) %>%
    plotly::layout(
      title = list(font = list(size = 16)),
      xaxis = list(
        title = list(
          text = "Module Membership (MM)",
          font = list(size = 14)
        ),
        tickfont = list(size = 12),
        showline = TRUE, linewidth = 2,
        linecolor = "black", mirror = TRUE
      ),
      yaxis = list(
        title = list(
          text = "Intramodular Connectivity (kWithin)",
          font = list(size = 14)
        ),
        tickfont = list(size = 12),
        showline = TRUE, linewidth = 2,
        linecolor = "black", mirror = TRUE
      ),
      plot_bgcolor = "white",
      paper_bgcolor = "white",
      margin = list(t = 100)
    )
})

output$hub_connect_mm_dt <- renderDT({
  req(hub.ds$plot_connect_mm)
  
  df <- connect_mm_data()
  rownames(df) <- NULL
  
  df <- df %>%
    dplyr::rename(
      "Intramodular Connectivity (kWithin)" = Connectivity,
      "MM" = AbsMM,
      "GS" = AbsGS,
      "Highlight Status" = Group
    )

  df_update <- sapply(df, is.numeric)
  df[df_update] <- round(df[df_update], 2)
  
  datatable(df, class = 'compact norap', selection = "none", rownames = FALSE, 
            # filter = list(position = 'top', clear = FALSE, plain = TRUE), 
            options = list(
              pageLength = 5, lengthMenu = c(5, 10, 25), scrollX = TRUE, search = list(regex = TRUE),
              order = list(list(2, "desc"), list(3, "desc")),
              columnDefs = list(
                list(targets = "_all", className = "dt-left")
              )    
            )
  )
})

output$hub_wgcna_network_ui <- renderUI({
  req(input$hub_wgcna_width, input$hub_wgcna_height)
  
  plotOutput(
    "hub_wgcna_network_plot",
    width = paste0(input$hub_wgcna_width, "px"),
    height = paste0(input$hub_wgcna_height, "px")
  )
})

observe({
  req(hub.ds$plot_wgcna_network, hub.ds$wgcna_g)
  color <- "#66C6FF"
  
  if (igraph::gorder(hub.ds$wgcna_g) == 0) {
    output$hub_wgcna_network_plot <- renderPlot({
      ggplot() + 
        annotate("text", x = 0.5, y = 0.5,
                 label = paste0("No nodes left after degree filtering.\n",
                                "Try lowering Degree Threshold.\n"),
                 size = 6) +
        theme_void()
    },  width = input$hub_wgcna_width, height = input$hub_wgcna_height)
  } else {
    num <- sum(hub.ds$wgcna_g %>% tidygraph::activate(nodes) %>% dplyr::pull(label) != "")
    
    output$hub_wgcna_network_plot <- renderPlot({
      ggraph(hub.ds$wgcna_g, layout = "fr") +
        geom_edge_link(alpha = 0.3, colour = "grey8") +
        geom_node_point(aes(size = Connectivity, fill = GS, alpha = MM),
                        shape = 21, colour = "black") +
        ggplot2::scale_size_continuous(range = c(4, 12), name = "Connectivity") + 
        geom_node_text(aes(label = label), repel = TRUE, size = input$hub_wgcna_tex_size) +
        scale_fill_gradient(low = "grey40", high = color) +
        theme_void(base_size = 14) +
        labs(
          title = "WGCNA network of module genes",
          subtitle = paste("Top", num, "hub genes labeled")
        )
    }, width = input$hub_wgcna_width, height = input$hub_wgcna_height)
  }
})

output$hub_connect_ppi_plot_ui <- renderUI({
  req(input$hub_connect_ppi_width, input$hub_connect_ppi_height)
  
  plotOutput(
    "hub_connect_ppi_plot",
    width = paste0(input$hub_connect_ppi_width, "px"),
    height = paste0(input$hub_connect_ppi_height, "px")
  )
})

observe({
  req(hub.ds$plot_connect_ppi)
  
  size <- 2
  color <- hub.ds$hub_module
  text_size <- as.numeric(input$hub_connect_ppi_text)
  
  df <- hub.ds$connect_ppi_df
  
  output$hub_connect_ppi_plot <- renderPlot({
    q <- ggplot(
      df,
      aes(x = Connectivity, y = degree, shape = type)
    ) +
      ylab("PPI Degree") +
      geom_point(
        data = subset(df, type != "Top Connectivity & Degree"),
        size = size,
        fill = color,
        color = color
      ) +
      geom_point(
        data = subset(df, type == "Top Connectivity & Degree"),
        size = size * 1.5,
        fill = color,
        color = "black",
        stroke = 0.4
      ) +
      ggrepel::geom_text_repel(
        data = subset(df, type == "Top Connectivity & Degree"),
        aes(label = Gene),
        size = text_size,
        color = "black",
        box.padding = 0.5,
        point.padding = 0.3,
        max.overlaps = Inf,
        seed = 123,
        show.legend = FALSE
      ) +
      scale_shape_manual(
        values = c(
          "Other" = 21,
          "Top Degree" = 24,
          "Top Connectivity" = 25,
          "Top Connectivity & Degree" = 23
        )
      ) +
      theme_classic()
    
    q
  })
})

output$hub_connect_ppi_dt <- renderDT({
  req(hub.ds$connect_ppi_df)
  df <- hub.ds$connect_ppi_df

  df <- df %>%
    dplyr::select(Gene, Connectivity, degree, everything()) %>%
    dplyr::rename(
      "Intramodular Connectivity (kWithin)" = Connectivity,
      "MM" = AbsMM,
      "GS" = AbsGS,
      "PPI Degree" = degree,
      "Type" = type
    )
  
  df_update <- sapply(df, is.numeric)
  df[df_update] <- round(df[df_update], 2)
  
  datatable(df, class = 'compact norap', selection = "none", rownames = FALSE, 
            # filter = list(position = 'top', clear = FALSE, plain = TRUE), 
            options = list(
              pageLength = 5, lengthMenu = c(5, 10, 25), scrollX = TRUE, search = list(regex = TRUE),
              order = list(list(1, "desc"), list(2, "desc")),
              columnDefs = list(
                list(targets = "_all", className = "dt-left")
              )    
            )
  )
})

observeEvent(input$hub_step6_2, {
  updateTabItems(session, "steps", "function")
  print("Step 6 button clicked, proceed to 6_Functional Enrichment Analysis.")
})
#####################
# download ####
observeEvent(input$hub_gsmm_run, {
  output$hub_download <- renderUI({
    screenshotButton(
      label = "Download Plot",
      download = TRUE, 
      server_dir = tempdir(), 
      id = "hub_gs2mmplot", 
      filename = "05_GSMMPlot")
  })
  
  output$hub_download2 <- downloadHandler(
    filename = function() {
      paste("05_kMEhub.csv")
    },
    content = function(file) {
      
      df <- hub.ds$hubgene_filter
      df$Gene <- rownames(df)
      df <- df %>%
        dplyr::select(MM, GS, Gene) %>%
        dplyr::rename(
          "Module Membership (MM)" = MM,
          "Gene Significance (GS)" = GS
        )
        
      numeric_cols <- sapply(df, is.numeric)
      df[numeric_cols] <- round(df[numeric_cols], 2)
      write.table(x = df, file = file, sep = ",", quote = FALSE, row.names = FALSE)
    }
  ) 
})

observeEvent(input$hub_run, {
  output$hub_download6 <- downloadHandler(
    filename = function() {
      paste0("05_ConnectMMPlot_", hub.ds$hub_module, ".png")
    },
    content = function(file) {
      ggplot2::ggsave(
        filename = file, plot = hub_connect_mm_ggplot(), device = "png",
        width = 10, height = 7,  units = "in", dpi = 300, bg = "white"
      )
    },
    contentType = "image/png"
  )
  
  output$hub_download7 <- downloadHandler(
    filename = function() {
      paste("05_ConnetMM.csv")
    },
    content = function(file) {
      
      df <- hub.ds$hub_gene_info
      rownames(df) <- NULL
      
      df <- df %>%
        dplyr::rename(
          "Intramodular Connectivity (kWithin)" = Connectivity,
          "MM" = AbsMM,
          "GS" = AbsGS,
          "Highlight Status" = Group
        )

      df_update <- sapply(df, is.numeric)
      df[df_update] <- round(df[df_update], 2)
      
      write.table(x = df, file = file, sep = ",", quote = FALSE, row.names = FALSE)
    }
  ) 
  
  output$hub_download3 <- downloadHandler(
    filename = function() { "05_Hubwgcnanetwork.pdf" },
    content = function(file) {
      plot_width_in  <- input$hub_wgcna_width / 96
      plot_height_in <- input$hub_wgcna_height / 96

      pdf(file, height = plot_height_in, width = plot_width_in, onefile = FALSE)
      color <- "#66C6FF"
      
      if (igraph::gorder(hub.ds$wgcna_g) == 0) {
        print(
          ggplot() + 
            annotate("text", x = 0.5, y = 0.5,
                     label = paste0("No nodes left after degree filtering.\n",
                                    "Try lowering Degree Threshold.\n"),
                     size = 6) +
            theme_void()
        )
      } else {
        num <- sum(hub.ds$wgcna_g %>% tidygraph::activate(nodes) %>% dplyr::pull(label) != "")
        
        print(
          ggraph(hub.ds$wgcna_g, layout = "fr") +
            geom_edge_link(alpha = 0.3, colour = "grey8") +
            geom_node_point(aes(size = degree, fill = GS, alpha = MM),
                            shape = 21, colour = "black") +
            geom_node_text(aes(label = label), repel = TRUE, size = input$hub_wgcna_tex_size) +
            scale_size(range = c(4, 12)) +
            scale_fill_gradient(low = "grey40", high = color) +
            theme_void(base_size = 14) +
            labs(
              title = "WGCNA network of module genes",
              subtitle = paste("Top", num, "hub genes labeled")
            )
        )
      }
      dev.off()
    }
  )
  
  output$hub_download5 <- renderUI({
    screenshotButton(
      label = "Download Plot",
      download = TRUE, 
      server_dir = tempdir(), 
      id = "hub_connect_ppi_plot", 
      filename = "05_ConnectPPI_Plot")
  })
  
  output$hub_download8 <- downloadHandler(
    filename = function() {
      paste("05_ConnetPPI.csv")
    },
    content = function(file) {
      
      df <- hub.ds$connect_ppi_df
      
      df <- df %>%
        dplyr::select(Gene, Connectivity, degree, everything()) %>%
        dplyr::rename(
          "Intramodular Connectivity (kWithin)" = Connectivity,
          "MM" = ModuleMembership,
          "GS" = GeneSignificance,
          "PPI Degree" = degree,
          "Type" = type
        )
      
      df_update <- sapply(df, is.numeric)
      df[df_update] <- round(df[df_update], 2)
      
      write.table(x = df, file = file, sep = ",", quote = FALSE, row.names = FALSE)
    }
  ) 
})
###################
## information ####
addPopover(
  id = "kmecut",
  options = list(
    content = paste0(
      "Filters genes by their correlation with the module eigengene. 
      Higher values select genes that better represent the module."
    ),
    title = "More information",
    placement = "right",
    trigger = "click",
    html = TRUE  # Ensure HTML content is parsed
  )
)

addPopover(
  id = "gscut",
  options = list(
    content = paste0(
      "Filters genes by their association with the selected trait. 
      Higher values select genes more strongly associated with the trait."
    ),
    title = "More information",
    placement = "right",
    trigger = "click",
    html = TRUE  # Ensure HTML content is parsed
  )
)

addPopover(
  id = "thres",
  options = list(
    content = paste0(
      "Weight refers to the connectivity between all genes in the gene group",br(),
      "you can select the appropriate weight threshold based on the 
      weight of the gene group you are interested in."
    ),
    title = "More information",
    placement = "right",
    trigger = "click",
    html = TRUE  # Ensure HTML content is parsed
  )
)
