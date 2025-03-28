create_lineage_heatmaps <- function(heat.2.df, 
                                  # heat.2.name= heat.2.name,
                                  legend.title.font = legend.title.font,
                                  legend.label.font = legend.label.font,
                                  annot.title.side = annot.title.side,
                                  show.sample.names = FALSE,
                                  rows.fs= rows.fs,
                                  cols.fs = 5
                                  ) {
  
  # B-cell maturation
  bcell_col_fun <- colorRamp2(c(-4, 0, 4), c("#f2d4f7", "white", "#7b3294"))
  
  # Lineage bias
  lineage_col_fun <- colorRamp2(c(-4, 0, 4), c("#92c5de", "white", "#f4a582" ))
  
  
  
 # Build the cellsorted scaled values
  #=====================================
  heatmap.2.A <- Heatmap(
    
    as.matrix(heat.2.df$cell.state),
    
    name = "Cell-population z-score",
    
    na_col = "white",
    
    heatmap_height = unit(3, "in"),
    
    #col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    #col = colorRamp2(c(-3, 0, 3), c("#7b3294", "white", "#d6604d")),
    
    show_column_names = show.sample.names,
    column_names_side =  "top",
    
    show_row_names = TRUE,
    
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    
    column_gap = unit(5, "mm"),
    row_names_side = "right", 
    row_names_max_width = unit(10, "cm"),
    
    # column_order = samples.order.mod,  # match oncoprint order!
    # column_split = split.cols.order,   # optional: match groupings
    
    row_names_gp = gpar(fontsize = rows.fs, fontface="bold"), # gene-names and percent (if not prc_gp is defined above)
    row_title_gp = gpar(fontsize =rows.fs+3, fontface = "bold"),
    
    column_names_gp = gpar(fontsize = cols.fs), 
    
    top_annotation = NULL,
    show_heatmap_legend = TRUE,
    
    heatmap_legend_param = list(
      title = "Cell-type enrichment\n(Relative Z-score)",
      
      title_gp = gpar(fontsize = legend.title.font, fontface = "bold"),
      title_position = annot.title.side, 
      labels_gp = gpar(fontsize = legend.label.font),
      
      grid_height= unit(1, "cm"), # size of the box in the legends
      grid_width= unit(1, "cm"),
      legend_height = unit(5, "cm"),
      legend_direction = "horizontal")
  )
  #===============================
  #  B-cell maturation index
  #===============================
  
  B.dev <- as.matrix(heat.2.df$dev.index)
  rownames(B.dev) <- "B-cell maturation index"
    
  heatmap.2.B <- Heatmap(
    
    B.dev, # B_dev_index
    
    # heatmap_height = unit(0.8, "in"),
    
    na_col = "white",
    
    # col = bcell_col_fun,
    
    # name = heat.2.name[2],
    
    show_column_names = show.sample.names,
    show_row_names = TRUE,
    row_names_side = "right", 

    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    
    column_gap = unit(5, "mm"),
    
    # row_title = "B-cell maturation index",
    # row_title_side= "right",
    
    # column_order = samples.order.mod,  # match oncoprint order!
    # column_split = split.cols.order,   # optional: match groupings
    
    ### FUTURTE DEV : add this: https://github.com/jokergoo/ComplexHeatmap/issues/116
    
    row_names_gp = gpar(fontsize = rows.fs, fontface="bold",col="blue"), # gene-names and percent (if not prc_gp is defined above)
    row_title_gp = gpar(fontsize =rows.fs+3, fontface = "bold"),
    
    row_names_max_width = unit(10, "cm"),
    
    
    column_names_gp = gpar(fontsize = cols.fs), 
    
    top_annotation = NULL,
    show_heatmap_legend = TRUE,
    
    heatmap_legend_param = list(
      title = "B-cell maturation\nLower: Early/progenitor-like,\nHigher: Mature/differentiated",
      
      title_gp = gpar(fontsize = legend.title.font, fontface = "bold"),
      title_position = annot.title.side, 
      labels_gp = gpar(fontsize = legend.label.font),
      
      grid_height= unit(1, "cm"), # size of the box in the legends
      grid_width= unit(1, "cm"),
      legend_height = unit(5, "cm"),
      legend_direction = "horizontal")
  )
  #===============================
  # lineage bias index
  #===============================
  ling.indx <- as.matrix(heat.2.df$ling.index)
  rownames(ling.indx) <- "Lineage bias index"
    
  heatmap.2.C <- Heatmap(
    
    ling.indx, # lineage_index
    # heatmap_height = unit(0.8, "in"),
    
    na_col = "white",
    
    # col = lineage_col_fun,
    
    # name = heat.2.name[2],
    
    show_column_names = show.sample.names,
    show_row_names = TRUE,
    row_names_side = "right",
    row_names_max_width = unit(10, "cm"),
    
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    
    # row_title = "Lineage bias index",
    # row_title_side= "right",
    
    column_gap = unit(5, "mm"),

    # column_order = samples.order.mod,  # match oncoprint order!
    # column_split = split.cols.order,   # optional: match groupings
    
    row_names_gp = gpar(fontsize = rows.fs, fontface="bold", col="blue"), # gene-names and percent (if not prc_gp is defined above)
    row_title_gp = gpar(fontsize =rows.fs+3, col="blue",fontface = "bold"),
    
    column_names_gp = gpar(fontsize = cols.fs), # sample-id fs
    
    top_annotation = NULL,
    show_heatmap_legend = TRUE,
    
    heatmap_legend_param = list(
      title = "Lineage bias\nLower: Myeloid-biased,\nHigher: Lymphoid-biased",
      
      title_gp = gpar(fontsize = legend.title.font, fontface = "bold"),
      title_position = annot.title.side, 
      labels_gp = gpar(fontsize = legend.label.font),
      
      grid_height= unit(1, "cm"), # size of the box in the legends
      grid_width= unit(1, "cm"),
      legend_height = unit(5, "cm"),
      legend_direction = "horizontal")
  )
  
  return(list(heat.1 = heatmap.2.A,
              heat.2 = heatmap.2.B,
              heat.3 = heatmap.2.C))
}