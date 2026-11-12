myenrichCluster <- function(object = NULL, type = c("BP", "MF", "CC", "KEGG"), OrgDb = NULL, id.trans = TRUE,
                            fromType = "SYMBOL", toType = c("ENTREZID"), all_genes = NULL, readable = TRUE, 
                            organism = "hsa", padjCutoff = 0.05, topn = 5, add.gene = FALSE, use_internal_data = FALSE,
                            seed = 2025825, stringdb = FALSE, string_org = NULL, string_version = NULL) {
  type <- match.arg(type)
  
  if(stringdb){
    string_db <- STRINGdb::STRINGdb$new(version = string_version, species = string_org, score_threshold = 400)
    
    if(type == "BP"){
      string_type <- "Process"
    } else if(type == "MF"){
      string_type <- "Function"
    } else if(type == "CC"){
      string_type <- "Component"
    } else {
      string_type <- "KEGG"
    }
  }
  
  if (!stringdb && (is.null(all_genes) || length(all_genes) == 0L)) {
    stop("Please provide background genes via `all_genes`.")
  }
  
  enrich.data <- object$wide.res
  
  if(stringdb == FALSE){
    if (id.trans) {
      universe <- clusterProfiler::bitr(
        all_genes, fromType = fromType, toType = toType, OrgDb = OrgDb
      )[[toType]]
    } else {
      universe <- all_genes
    }
    universe <- unique(universe)
  }
  
  result <- purrr::map_df(seq_along(unique(enrich.data$cluster)), function(x) {
    
    # filter
    tmp <- enrich.data %>%
      dplyr::filter(cluster == unique(enrich.data$cluster)[x])
    
    # enrich
    if(stringdb == FALSE){
      if(id.trans == TRUE){
        gene.ent <- clusterProfiler::bitr(
          tmp$gene, fromType = fromType, toType = toType, OrgDb = OrgDb
        )
        target.gene = unlist(gene.ent[,c(toType)])
      }else{
        target.gene = tmp$gene
      }
    } else {
      target.gene <- unique(as.character(tmp$gene))
    }
    
    if (length(target.gene) == 0L) return(NULL)
    
    if(type %in% c("BP","MF","CC")){
      if(stringdb == FALSE){
        set.seed(seed)
        ego <- clusterProfiler::enrichGO(
          gene = target.gene, keyType = toType, OrgDb = OrgDb, ont = type,
          pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = readable,
          universe = universe
        )
      } else {
        ego <- string_db$get_enrichment(target.gene, category = string_type)
      }
      
    } else if (type %in% "KEGG") {
      
      if(stringdb == FALSE){
        set.seed(seed)
        ego <- clusterProfiler::enrichKEGG(
          gene = target.gene, keyType = "kegg", organism = organism,
          pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff  = 1,
          use_internal_data = use_internal_data, universe = universe
        )
        
        # transform gene id
        if(readable == TRUE){
          ego <- clusterProfiler::setReadable(ego, OrgDb = OrgDb, keyType = toType)
        }else{
          ego <- ego
        }
      } else{
        ego <- string_db$get_enrichment(target.gene, category = string_type)
      }
    }
    
    # to data.frame
    enrich_res <- data.frame(ego)
    if(nrow(enrich_res) > 0){
      if(stringdb == FALSE){
        df <- as.data.frame(ego) %>%
          dplyr::filter(p.adjust < padjCutoff) %>%
          dplyr::mutate(group = paste0("C", unique(enrich.data$cluster)[x])) %>%
          dplyr::arrange(dplyr::desc(FoldEnrichment))
        
        if (nrow(df) == 0) return(NULL)
        
        df <- df %>%
          dplyr::mutate(
            ratio = vapply(
              strsplit(as.character(GeneRatio), "/", fixed = TRUE),
              function(z) as.numeric(z[1]) / as.numeric(z[2]) * 100,
              numeric(1)
            )
          )
      } else {
        n_input <- length(target.gene)
        
        df <- as.data.frame(ego) %>%
          dplyr::mutate(
            group = paste0("C", unique(enrich.data$cluster)[x]),
            Description = description,
            p.adjust = fdr,
            FoldEnrichment = NA_real_,
            Count = number_of_genes,
            GeneRatio = Count / n_input,
            ratio = GeneRatio * 100
          ) %>%
          dplyr::filter(p.adjust < padjCutoff) %>%
          dplyr::arrange(p.adjust)
        
        
        if (nrow(df) == 0) return(NULL)
      }
      
      # whether save all res
      if(!is.null(topn)){
        top <- if (length(topn) == 1) topn else topn[x]
        df <- dplyr::slice_head(df, n = top)
      }
      
      
      if (isTRUE(add.gene)) {
        df1 <- df %>%
          dplyr::select(group, Description, p.adjust, FoldEnrichment, Count, ratio, geneID)
      } else {
        df1 <- df %>%
          dplyr::select(group, Description, p.adjust, FoldEnrichment, Count, ratio)
      }
      
      # results
      return(df1)
    }else{
      return(NULL)
    }
  })
  
  if (nrow(result) == 0L) {
    result <- data.frame(
      group = character(),
      Description = character(),
      p.adjust = numeric(),
      FoldEnrichment = numeric(),
      Count = numeric(),
      ratio = numeric()
    )
    
    if (isTRUE(add.gene)) {
      result$geneID <- character()
    }
  }
  
  return(result)
}

myvisCluster <- function(object = NULL,
                         ht.col.list = list(col_range = c(-2, 0, 2), col_color = c("#08519C", "white", "#A50F15")),
                         border = TRUE, mline.col = "#CC3333", ctAnno.col = NULL, set.md = "median",
                         textbox.pos = c(0.5, 0.8), textbox.size = 8,
                         panel.arg = c(2, 0.25, 4, "grey90", NA),
                         annoTerm.data = NULL, annoTerm.mside = "right", termAnno.arg = c("grey95", "grey50"),
                         add.bar = FALSE, bar.width = 8, textbar.pos = c(0.8, 0.8),
                         go.col = NULL, go.size = NULL, by.go = "anno_link",
                         annoKegg.data = NULL, annoKegg.mside = "right", keggAnno.arg = c("grey95", "grey50"),
                         add.kegg.bar = FALSE, 
                         kegg.col = NULL, kegg.size = NULL, by.kegg = "anno_link",
                         word_wrap = TRUE, add_new_line = TRUE,
                         add.box = FALSE, boxcol = NULL, box.arg = c(0.1, "grey50"),
                         add.point = FALSE, point.arg = c(19, "orange", "orange", 1),
                         add.line = TRUE, line.side = "right",
                         markGenes = NULL, markGenes.side = "right", genes.gp = c("italic", 10, NA),
                         term.text.limit = c(10, 18), mulGroup = NULL, lgd.label = NULL,
                         show_row_names = FALSE, subgroup.anno = NULL, annnoblock.text = TRUE,
                         annnoblock.gp = c("white", 8), add.sampleanno = TRUE,
                         sample.group = NULL, sample.col = NULL, sample.order = NULL,
                         cluster.order = NULL, HeatmapAnnotation = NULL, column.split = NULL, cluster_columns = FALSE,
                         pseudotime_col = NULL, row_annotation_obj = NULL, bar_bary = "padj", ...) {
  ComplexHeatmap::ht_opt(message = FALSE)
  
  if(is.null(ht.col.list[["col_range"]])){
    col_range = c(-2,0,2)
  }else{
    col_range = ht.col.list[["col_range"]]
  }
  
  if(is.null(ht.col.list[["col_color"]])){
    col_color = c("#08519C", "white", "#A50F15")
  }else{
    col_color = ht.col.list[["col_color"]]
  }
  
  col_fun = circlize::colorRamp2(col_range, col_color)
  
  data <- data.frame(object$wide.res,check.names = FALSE) %>%
    dplyr::arrange(as.numeric(as.character(cluster)))
  
  # prepare matrix, type = wgcna
  mat <- data %>%
    dplyr::arrange(as.numeric(as.character(cluster))) %>%
    dplyr::select(-gene,-cluster,-modulecol)
  
  rownames(mat) <- data$gene
  
  if(!is.null(sample.order)){
    mat <- mat[,sample.order]
  }
  
  # split info
  cl.info <- data.frame(table(data$cluster)) %>%
    dplyr::mutate(Var1 = as.numeric(as.character(Var1))) %>%
    dplyr::arrange(Var1)
  cluster.num <- nrow(cl.info)
  
  subgroup <- lapply(1:nrow(cl.info),function(x){
    nm <- rep(as.character(cl.info$Var1[x]),cl.info$Freq[x])
    paste("C",nm,sep = '')
  }) %>% unlist()
  
  # cluster orders
  if(!is.null(cluster.order)){
    subgroup <- factor(subgroup,levels = paste("C",cluster.order,sep = ""))
    cluster_row_slices = FALSE
  }else{
    cluster_row_slices = TRUE
  }
  
  # plot
  if(is.null(sample.group)){
    sample.info = colnames(mat)
    
    # split columns
    if(is.null(HeatmapAnnotation)){
      column_split = NULL
    }else{
      column_split = column.split
    }
  }else{
    sample.info = sample.group
    
    # split columns
    column_split = factor(sample.group,levels = unique(sample.group))
  }
  
  # sample colors
  # if no pseudotime_col
  
  sample.info <- factor(sample.info,levels = unique(sample.info))
  # sample colors
  if(is.null(sample.col)){
    scol <- circlize::rand_color(n = length(sample.info))
    names(scol) <- sample.info
  }else{
    scol <- sample.col
    names(scol) <- sample.info
  }
  
  # top anno
  if(add.sampleanno == TRUE){
    if(is.null(HeatmapAnnotation)){
      topanno = ComplexHeatmap::HeatmapAnnotation(
        sample = sample.info, col = list(sample = scol), 
        gp = grid::gpar(col = "white"),
        show_legend = TRUE,
        show_annotation_name = FALSE
      )
    }else{
      topanno = HeatmapAnnotation
    }
  }else{
    topanno = NULL
  }
  
  # bar annotation for clusters
  if(is.null(ctAnno.col)){
    colanno <- circlize::rand_color(n = cluster.num)
  }else{
    colanno <- ctAnno.col
  }
  
  names(colanno) <- as.character(cl.info$Var1)
  align_to = split(1:nrow(mat), subgroup)
  
  anno.block <- ComplexHeatmap::anno_block(
    align_to = align_to,
    panel_fun = function(index, nm) {
      # npos = as.numeric(unlist(strsplit(nm,split = "C"))[2])
      # rect
      cluster_id <- sub("^C", "", nm)
      grid::grid.rect(gp = grid::gpar(fill = unname(colanno[cluster_id]),col = NA))
      # text
      if(annnoblock.text == TRUE){
        grid::grid.text(
          label = paste("n:",length(index),sep = ''), rot = 90, 
          gp = grid::gpar(col = annnoblock.gp[1], fontsize = as.numeric(annnoblock.gp[2]))
        )
      }
    }, which = "row")
  
  # gene annotation for heatmap
  # whether mark your genes on plot
  if(!is.null(markGenes)){
    # all genes
    rowGene <- rownames(mat)
    
    # tartget gene
    annoGene <- markGenes
    
    # add color for gene
    gene.col <- data %>%
      dplyr::select(gene,cluster) %>%
      dplyr::filter(gene %in% annoGene)
    
    # change method
    changed_clusters <- cl.info$Var1
    
    gene.col <- purrr::map_df(seq_along(changed_clusters), function(i) {
      cl <- changed_clusters[i]
      gene.col |> 
        dplyr::filter(cluster == cl) |> 
        dplyr::mutate(col = unname(colanno[as.character(cl)]))
    })
    
    gene.col <- gene.col[match(annoGene,gene.col$gene),]
    
    if(is.na(genes.gp[3])){
      gcol = gene.col$col
    }else{
      gcol = genes.gp[3]
    }
    
    # get target gene index
    index <- match(annoGene,rowGene)
    
    # some genes annotation
    geneMark = ComplexHeatmap::anno_mark(
      at = index, labels = annoGene, which = "row", side = markGenes.side, 
      labels_gp = grid::gpar(fontface = genes.gp[1], fontsize = as.numeric(genes.gp[2]), col = gcol)
    )
  }else{
    geneMark = NULL
  }
  
  # final annotation for heatmap
  right_annotation = ComplexHeatmap::rowAnnotation(gene = geneMark,cluster = anno.block)
  
  # return plot according to plot type
  if(ncol(mat) > 20){
    show_column_names = FALSE
  }else{
    show_column_names = TRUE
  }
  
  lgd_list <- NULL
  
  # heatmap + line
  rg = range(mat)
  anno_ggplot2 <- NULL
  
  # panel_fun for line plot
  panel_fun = function(index, nm) {
    
    # whether add boxplot
    if(add.box == TRUE & add.line != TRUE){
      xscale = c(-0.1,1.1)
    }else{
      # xscale = c(0,1)
      xscale = c(-0.1,1.1)
      panel_scale <- c(0.1,0.9)
    }
    
    grid::pushViewport(grid::viewport(xscale = xscale, yscale = c(0,1)))
    grid::grid.rect()
    
    
    if(is.null(mulGroup)){
      mulGroup <- ncol(mat)
      
      # ================ calculate group columns index
      seqn <- data.frame(st = 1,sp = ncol(mat))
    }else{
      mulGroup <- mulGroup
      
      grid::grid.lines(x = c(0,1),y = rep(0.5,2),
                       gp = grid::gpar(col = "black",lty = "dashed"))
      
      # ================ calculate group columns index
      cu <- cumsum(mulGroup)
      seqn <- data.frame(st = c(1,cu[1:(length(cu) - 1)] + 1),
                         sp = c(cu[1],cu[2:length(cu)]))
    }
    
    
    # ===========================================================================
    # multiple lines
    lapply(1:nrow(seqn), function(x){
      tmp <- seqn[x,]
      tmpmat <- mat[index, c(tmp$st:tmp$sp)]
      
      # choose method
      if(set.md == "mean"){
        mdia <- colMeans(tmpmat)
      }else if(set.md == "median"){
        mdia <- apply(tmpmat, 2, stats::median)
      }else{
        message("supply mean/median !")
      }
      
      # boxplot xpos
      pos = scales::rescale(1:ncol(tmpmat),to = c(0,1))
      
      # boxcol
      if(is.null(boxcol)){
        boxcol <- rep("grey90",ncol(tmpmat))
      }else{
        boxcol <- boxcol
      }
      
      # boxplot grobs
      if(add.box == TRUE){
        lapply(1:ncol(tmpmat), function(x){
          ComplexHeatmap::grid.boxplot(
            scales::rescale(tmpmat[,x], to = c(0,1), from = c(rg[1] - 0.5,rg[2] + 0.5)),
            pos = pos[x], direction = "vertical", box_width = as.numeric(box.arg[1]),
            outline = FALSE, gp = grid::gpar(col = box.arg[2],fill = boxcol[x])
          )
        })
      }
      
      # points grobs
      if(add.point == TRUE){
        grid::grid.points(
          x = scales::rescale(1:ncol(tmpmat),to = panel_scale),
          y = scales::rescale(mdia,to = c(0,1),from = c(rg[1] - 0.5,rg[2] + 0.5)),
          pch = as.numeric(point.arg[1]),
          gp = grid::gpar(fill = point.arg[2],col = point.arg[3]),
          size = grid::unit(as.numeric(point.arg[4]), "char")
        )
      }
      
      # lines grobs
      if(add.line == TRUE){
        grid::grid.lines(
          x = scales::rescale(1:ncol(tmpmat),to = panel_scale),
          y = scales::rescale(mdia,to = c(0,1),from = c(rg[1] - 0.5,rg[2] + 0.5)),
          gp = grid::gpar(lwd = 3,col = mline.col[x])
        )
      }
    })
    
    
    # get gene numbers
    grid.textbox <- utils::getFromNamespace("grid.textbox", "ComplexHeatmap")
    
    text <- paste("Gene Size:",nrow(mat[index, ]),sep = ' ')
    grid.textbox(text,x = textbox.pos[1],y = textbox.pos[2],
                 gp = grid::gpar(fontsize = textbox.size,
                                 fontface = "italic",
                                 ...))
    
    grid::popViewport()
  }
  
  # whether annotate subgroups
  if(!is.null(subgroup.anno)){
    align_to = split(1:nrow(mat), subgroup)
    align_to = align_to[subgroup.anno]
  }else{
    align_to = subgroup
  }
  
  # anno link annotation
  anno = ComplexHeatmap::anno_link(
    align_to = align_to, which = "row", panel_fun = panel_fun, 
    size = grid::unit(as.numeric(panel.arg[1]), "cm"),
    gap = grid::unit(as.numeric(panel.arg[2]), "cm"),
    width = grid::unit(as.numeric(panel.arg[3]), "cm"), 
    side = line.side, link_gp = grid::gpar(fill = panel.arg[4],col = panel.arg[5])
  )
  
  # add go term annotations
  # load term info
  termanno <- annoTerm.data
  colnames(termanno) <- c("id", "term", "padj", "foldenrichment", "count", "ratio")
  
  # term colors
  if(is.null(go.col)){
    gocol <- circlize::rand_color(n = nrow(termanno))
  }else{
    gocol <- go.col
  }
  
  # term text size
  if(is.null(go.size)){
    gosize <- rep(12,nrow(termanno))
  }else{
    if(go.size == "padj"){
      # loop for re-scaling pvalue
      purrr::map_df(unique(termanno$id),function(x){
        tmp <- termanno %>%
          dplyr::filter(id == x) %>%
          dplyr::mutate(size = scales::rescale(-log10(padj),to = term.text.limit))
      }) -> termanno.tmp
      
      gosize <- termanno.tmp$size
    }else{
      gosize <- go.size
    }
  }

  # add to termanno
  termanno <- termanno %>%
    dplyr::ungroup() %>%
    dplyr::mutate(col = gocol, fontsize = gosize)
  
  # to list
  lapply(1:length(unique(termanno$id)), function(x){
    tmp = termanno[which(termanno$id == unique(termanno$id)[x]),]
    df <- data.frame(
      text = tmp$term, col = tmp$col, fontsize = tmp$fontsize, 
      fontfamily = "sans",   # 指定字型
      fontface = "plain" 
    )
    return(df)
  }) -> term.list
  
  # add names
  names(term.list) <- unique(termanno$id)
  
  term.list <- lapply(term.list, function(df) {
    result <- lapply(seq_len(nrow(df)), function(i) {
      term_text <- as.character(df$text[i])
      exdent <- if (isTRUE(startsWith(trimws(term_text), "KEGG"))) 12L else 7L
      
      lines <- strwrap(term_text, width = 45, exdent = exdent)
      
      out <- df[rep(i, length(lines)), , drop = FALSE]
      out$text <- lines
      out
    })
    
    do.call(rbind, result)
  })
  
  # whether annotate subgroups
  if(!is.null(subgroup.anno)){
    align_to2 = split(seq_along(subgroup), subgroup)
    align_to2 = align_to2[subgroup.anno]
    
    term.list = term.list[subgroup.anno]
  }else{
    align_to2 = subgroup
    term.list = term.list
  }
  
  textbox = ComplexHeatmap::anno_textbox(
    align_to2, term.list, word_wrap = word_wrap, add_new_line = add_new_line, max_width = unit(100, "mm"),
    side = annoTerm.mside, text_space = grid::unit(2, "mm"), line_space = grid::unit(1.5, "mm"),
    background_gp = grid::gpar(fill = termAnno.arg[1], col = termAnno.arg[2]),
    by = by.go
  )
  
  # GO bar anno function
  anno_gobar <- function(data = NULL, bar.width = 0.1, align_to = NULL, panel.arg = panel.arg, ...){
    # process data
    if(bar_bary == "padj"){
      data <- data %>%
        dplyr::mutate(bary = -log10(padj))
    }else{
      data <- data %>%
        dplyr::mutate(bary = ratio)
    }
    
    ComplexHeatmap::anno_zoom(
      align_to = align_to, which = "row",
      panel_fun = function(index,nm){
        grid::pushViewport(grid::viewport(xscale = c(0,1),yscale = c(0,1)))
        grid::grid.rect()
        # sub data
        tmp <- data %>%
          dplyr::filter(id == nm)
        
        grid::grid.segments(
          x0 = rep(0,nrow(tmp)),
          x1 = scales::rescale(rev(tmp$bary),to = c(0.1,0.9)),
          y0 = scales::rescale(1:nrow(tmp),to = c(0.1,0.9)),
          y1 = scales::rescale(1:nrow(tmp),to = c(0.1,0.9)),
          gp = grid::gpar(lwd = bar.width, col = rev(tmp$col), lineend = "butt")
        )
        
        # add cluster name
        grid.textbox <- utils::getFromNamespace("grid.textbox", "ComplexHeatmap")
        text <- nm
        grid.textbox(
          text,
          x = textbar.pos[1],y = textbar.pos[2],
          gp = grid::gpar(fontsize = textbox.size,
                          fontface = "italic",
                          col = unique(tmp$col),
                          ...))
        grid::popViewport()
      },
      size = grid::unit(as.numeric(panel.arg[1]), "cm"),
      gap = grid::unit(as.numeric(panel.arg[2]), "cm"),
      width = grid::unit(as.numeric(panel.arg[3]), "cm"),
      side = "right", link_gp = grid::gpar(fill = termAnno.arg[1],col = termAnno.arg[2]),
      ...
    )
  }
  
  # ================================
  # bar anno
  baranno = anno_gobar(
    data = termanno,
    align_to = align_to2,
    panel.arg = panel.arg,
    bar.width = bar.width
  )
  # ============================================================================================
  # whether add bar annotation
  # ============================================================================================
  if(add.bar == TRUE){
    baranno
  }else{
    baranno = NULL
  }
  
  # whether add kegg term annotations
  if(!is.null(annoKegg.data)){
    # load term info
    termanno <- annoKegg.data
    colnames(termanno) <- c("id", "term", "padj", "foldenrichment", "count", "ratio")
    
    # term colors
    if(is.null(kegg.col)){
      gocol <- circlize::rand_color(n = nrow(termanno))
    }else{
      gocol <- kegg.col
    }
    
    # term text size
    if(is.null(kegg.size)){
      gosize <- rep(12,nrow(termanno))
    }else{
      if(kegg.size == "padj"){
        # loop for re-scaling pvalue
        purrr::map_df(unique(termanno$id),function(x){
          tmp <- termanno %>%
            dplyr::filter(id == x) %>%
            dplyr::mutate(size = scales::rescale(-log10(padj),to = term.text.limit))
        }) -> termanno.tmp
        
        gosize <- termanno.tmp$size
      }else{
        gosize <- kegg.size
      }
    }
    
    # add to termanno
    termanno <- termanno %>%
      dplyr::ungroup() %>%
      dplyr::mutate(col = gocol, fontsize = gosize)
    
    # to list
    lapply(1:length(unique(termanno$id)), function(x){
      tmp = termanno[which(termanno$id == unique(termanno$id)[x]),]
      df <- data.frame(
        text = tmp$term, col = tmp$col, fontsize = tmp$fontsize, 
        fontfamily = "sans",   # 指定字型
        fontface = "plain" 
      )
      return(df)
    }) -> term.list
    
    # add names
    names(term.list) <- unique(termanno$id)
    
    term.list <- lapply(term.list, function(df) {
      result <- lapply(seq_len(nrow(df)), function(i) {
        lines <- strwrap(
          as.character(df$text[i]),
          width = 45,
          exdent = 12
        )
        
        out <- df[rep(i, length(lines)), , drop = FALSE]
        out$text <- lines
        out
      })
      
      do.call(rbind, result)
    })
    
    # whether annotate subgroups
    if(!is.null(subgroup.anno)){
      align_to2 = split(seq_along(subgroup), subgroup)
      align_to2 = align_to2[subgroup.anno]
      
      term.list = term.list[subgroup.anno]
    }else{
      align_to2 = subgroup
      term.list = term.list
    }
    
    # anno_textbox
    textbox.kegg = ComplexHeatmap::anno_textbox(
      align_to2, term.list, word_wrap = word_wrap, add_new_line = add_new_line, max_width = unit(100, "mm"),
      side = annoKegg.mside, text_space = grid::unit(2, "mm"), line_space = grid::unit(1.5, "mm"),
      background_gp = grid::gpar(fill = keggAnno.arg[1], col = keggAnno.arg[2]),
      by = by.kegg
    )
    
    # KEGG bar anno function
    anno_keggbar <- function(data = NULL, bar.width = 0.1, align_to = NULL, panel.arg = panel.arg, ...){
      # process data
      if(bar_bary == "padj"){
        data <- data %>%
          dplyr::mutate(bary = -log10(padj))
      }else{
        data <- data %>%
          dplyr::mutate(bary = ratio)
      }
      
      ComplexHeatmap::anno_zoom(
        align_to = align_to, which = "row",
        panel_fun = function(index,nm){
          grid::pushViewport(grid::viewport(xscale = c(0,1),yscale = c(0,1)))
          grid::grid.rect()
          # sub data
          tmp <- data %>%
            dplyr::filter(id == nm)
          grid::grid.segments(
            x0 = rep(0,nrow(tmp)),
            x1 = scales::rescale(rev(tmp$bary),to = c(0.1,0.9)),
            y0 = scales::rescale(1:nrow(tmp),to = c(0.1,0.9)),
            y1 = scales::rescale(1:nrow(tmp),to = c(0.1,0.9)),
            gp = grid::gpar(lwd = bar.width, col = rev(tmp$col), lineend = "butt")
          )
          # add cluster name
          grid.textbox <- utils::getFromNamespace("grid.textbox", "ComplexHeatmap")
          text <- nm
          grid.textbox(
            text, x = textbar.pos[1], y = textbar.pos[2],
            gp = grid::gpar(fontsize = textbox.size, fontface = "italic", col = unique(tmp$col), ...)
          )
          grid::popViewport()
        },
        size = grid::unit(as.numeric(panel.arg[1]), "cm"),
        gap = grid::unit(as.numeric(panel.arg[2]), "cm"),
        width = grid::unit(as.numeric(panel.arg[3]), "cm"),
        side = "right",
        link_gp = grid::gpar(fill = keggAnno.arg[1], col = keggAnno.arg[2]),
        ...
      )
    }
    
    # bar anno
    baranno.kegg = anno_keggbar(
      data = termanno, align_to = align_to2, panel.arg = panel.arg, bar.width = bar.width
    )
    # whether add bar annotation
    if(add.kegg.bar == TRUE){
      baranno.kegg
    }else{
      baranno.kegg = NULL
    }
  }else{
    textbox.kegg = NULL
    baranno.kegg = NULL
  }
  
  # final row annotations
  if(line.side == "right"){
    if(markGenes.side == "right"){
      right_annotation2 = ComplexHeatmap::rowAnnotation(
        gene = geneMark, cluster = anno.block, line = anno, anno_ggplot2 = anno_ggplot2,
        textbox = textbox, bar = baranno, textbox.kegg = textbox.kegg, baranno.kegg = baranno.kegg
      )
      
      if(!is.null(row_annotation_obj)){
        left_annotation <- row_annotation_obj
      }else{
        left_annotation <- NULL
      }
    }else{
      right_annotation2 = ComplexHeatmap::rowAnnotation(
        cluster = anno.block, line = anno, anno_ggplot2 = anno_ggplot2, textbox = textbox,
        bar = baranno, textbox.kegg = textbox.kegg, baranno.kegg = baranno.kegg
      )
      
      if(!is.null(row_annotation_obj)){
        left_annotation <- do.call(
          ComplexHeatmap::rowAnnotation, modifyList(list(gene = geneMark), row_annotation_obj)
        )
      }else{
        left_annotation <- ComplexHeatmap::rowAnnotation(gene = geneMark)
      }
    }
    
  } else {
    if(markGenes.side == "right"){
      right_annotation2 = ComplexHeatmap::rowAnnotation(
        gene = geneMark, cluster = anno.block, anno_ggplot2 = anno_ggplot2,
        textbox = textbox, bar = baranno, textbox.kegg = textbox.kegg, baranno.kegg = baranno.kegg
      )
      
      if(!is.null(row_annotation_obj)){
        left_annotation <- do.call(
          ComplexHeatmap::rowAnnotation, modifyList(list(line = anno), row_annotation_obj)
        )
      }else{
        left_annotation <- ComplexHeatmap::rowAnnotation(line = anno)
      }
    }else{
      right_annotation2 = ComplexHeatmap::rowAnnotation(
        cluster = anno.block, anno_ggplot2 = anno_ggplot2, textbox = textbox,
        bar = baranno, textbox.kegg = textbox.kegg, baranno.kegg = baranno.kegg
      )
      
      if(!is.null(row_annotation_obj)){
        left_annotation <- do.call(
          ComplexHeatmap::rowAnnotation, modifyList(list(gene = geneMark,line = anno), row_annotation_obj)
        )
      }else{
        left_annotation <- ComplexHeatmap::rowAnnotation(line = anno,gene = geneMark)
      }
    }
  }
  
  # save
  
  htf <- ComplexHeatmap::Heatmap(
    as.matrix(mat), name = "Z-score", cluster_columns = cluster_columns,
    show_row_names = show_row_names, border = border, column_split = column_split,
    top_annotation = topanno, right_annotation = right_annotation2,
    left_annotation = left_annotation,
    column_names_side = "top", show_column_names = show_column_names,
    row_split = subgroup, cluster_row_slices = cluster_row_slices,
    col = col_fun, ...
  )
  
  # draw lines legend
  if(is.null(mulGroup)){
    ComplexHeatmap::draw(htf, merge_legend = TRUE, annotation_legend_list = lgd_list)
  }else{
    if(is.null(lgd.label)){
      lgd.label <- paste("group", 1:length(mulGroup), sep = '')
    }else{
      lgd.label <- lgd.label
    }
    
    lgd_list2 = ComplexHeatmap::Legend(
      labels = lgd.label, type = "lines", legend_gp = grid::gpar(col = mline.col, lty = 1)
    )
    
    if(!is.null(lgd_list)){
      lgd_list_com <- ComplexHeatmap::packLegend(lgd_list,lgd_list2)
    }else{
      lgd_list_com <- lgd_list2
    }
    
    ComplexHeatmap::draw(htf, annotation_legend_list = lgd_list_com, merge_legend = TRUE)
  }
}
