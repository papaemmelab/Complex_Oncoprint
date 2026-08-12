generate_surv_heat <- function(heatmap.mat, OS, heat.name = NULL, 
                                  barplot.fs = 10, 
                                  right.w= 2,
                                  annot.grp.fs= 10,
                                  annotation_name_side = "top"
                                  ){

  
  if (!(all(rownames(heatmap.mat)==rownames(OS)) & all(colnames(heatmap.mat)==colnames(OS)))) {
    stop("\n OS order is different from mat!")
  }
  
  #=================================================================
  # Right heatmap of OS 
  #=================================================================
  
  # browser()
  
  heat.2 <- ComplexHeatmap::Heatmap(OS, name = heat.name, 
                                    na_col = "grey",
                                    cluster_rows = FALSE,
                                    cluster_columns = FALSE)
  
  # browser()
  
  ha <- rowAnnotation(OS_mutated_pts = anno_boxplot(
    OS,
    
    axis_param = list(side= "top",
                      gp= gpar(fontsize= barplot.fs, fontface="bold")),
    box_width = 0.6,
    width = unit(right.w, "cm"),  # ✅ This controls the sidebar width
    
    outline = FALSE,
    gp= gpar(fill = "lightgrey", fontsize= barplot.fs, fontface="bold"),
    
    fun = function(x) {
      x <- na.omit(x)
      if (length(x) < 1) return(numeric(0))  # gracefully skip
      x
    }
  ),
  annotation_name_rot = 0,
  annotation_name_side = annotation_name_side, ### this controls where you see OS_time labels
  annotation_name_gp= gpar(fontsize= annot.grp.fs, fontface="bold", col="black") ### this controls OS_time labels etc
  
  )
  
  return(ha)
}