draw_basic_oncoprint <-  function(M, EFFECT, alter_fun, 
                                  list.colors, top_annotation, 
                                  
                                  heatmap.legend.list,
                                  annot.legend.list,
                                  
                                  heatmap.legend.side= "right",
                                  annot.legend.side= "bottom",
                                  
                                  column_order= NULL,
                                  right.w= 13, 
                                  LABS, 
                                  
                                  saveFile= NA,
                                  w=3200, h=1800,
                                  num.rows.heatmap.lgd= NA,

                                  font.obj= NA, 
                                  
                                  fig.title= NULL,
                                  show.border= TRUE, show.sample.names= TRUE){

  
  h1 <- top_annotation
  
  
  simple.ht <- oncoPrint(M, get_type = function(x) strsplit(x, ";")[[1]],
                         
                         alter_fun = alter_fun, col = append(list.colors$mut.colors, list.colors$cyto.colors),
                         
                         column_order = column_order,
                         
                         remove_empty_columns = FALSE,
                         
                         show_column_names = show.sample.names,
                         
                         # === Gene barplots on the left ====
                         
                         top_annotation = h1,
                         
                         right_annotation = rowAnnotation(row_bar = anno_oncoprint_barplot(type= NULL,
                                                                                           border= show.border, 
                                                                                           axis_param = list(side= "top", gp= gpar(fontsize= font.obj$barplot.font, fontface="bold"))),
                                                          annotation_width= unit(right.w,"cm")),   ## controls the width of the row.barplots
                         
                         #show_row_barplot = TRUE, # obsolete param
                         #row_barplot_width = unit(right.w, "cm"), # obsolete param
                         
                         split= LABS,
                         
                         # ==========================================
                         # ==========================================
                         
                         # === Title ====
                         column_title = fig.title,
                         column_title_gp = gpar(fontsize = font.obj$fig.title.font, fontface = "bold"), # title font-size
                         gap = unit(10, "mm"),
                         
                         # === Column/Sample names ====
                         column_names_gp = gpar(cex=1, col= "black", fontsize = font.obj$cols.font, fontface="bold"), #default size = 18
                         column_names_max_height= unit(15,"cm") , # adjust this to control the name of samples (col names)
                         
                         # === Percent ====
                         pct_gp=gpar(fontsize = font.obj$pct.font, fontface = "bold", col="black"), # specific control over percentage info on the left (add col="blue" to change colors)
                         row_names_gp = gpar(fontsize = font.obj$rows.font, fontface="bold"), # gene-names and percent (if not prc_gp is defined above)
                         row_title_gp = gpar(fontsize =font.obj$rows.font, col="blue",fontface = "bold"),
                         
                         # === Percent ====
                         # layer_fun = function(j, i, x, y, w, h, fill) {
                         #   ind_mat = restore_matrix(j, i, x, y)
                         #   ind = unique(c(ind_mat[2, ], ind_mat[, 3]))
                         #   grid.points(x[ind], y[ind], pch = 16, size = unit(4, "mm"))},
                         
                         # === Legend ====
                         
                         heatmap_legend_param = list(title = "Alterations", at = EFFECT$variants,
                                                     labels = EFFECT$labels,
                                                     title_gp = gpar(fontsize = font.obj$legend.title.font, fontface="bold"),
                                                     # title_position = main.legend.pos, 
                                                     title_position = "topleft", 
                                                     labels_gp = gpar(fontsize = font.obj$legend.label.font),
                                                     grid_height= unit(1, "cm"),
                                                     nrow=num.rows.heatmap.lgd,
                                                     grid_width= unit(1, "cm"),
                                                     legend_height = unit(20, "cm"))) 

##======================================================
## Draw simple.ht ====
##======================================================

  jpeg(saveFile, width=w, height=h, pointsize =14, res = 100)
  
  if (heatmap.legend.side== annot.legend.side){
    to.merge.param= TRUE
    } else {
    to.merge.param= FALSE
    }
  
    draw(simple.ht, split= LABS,  
         merge_legend = to.merge.param, 
         
         heatmap_legend_list = heatmap.legend.list,
         heatmap_legend_side = heatmap.legend.side, 
         
         annotation_legend_list= annot.legend.list,
         annotation_legend_side = annot.legend.side
    )
  
  dev.off()
  
  return(simple.ht)
}


