prepare_TOP_annotation <- function(list.colors,show.border,axis.side,barplot.font, legend.title.font, top.w, show.survival=FALSE,...){
  
  
  if (show.survival) {
    h1 = HeatmapAnnotation(column_bar = anno_oncoprint_barplot(type= NULL,
                                                               border= show.border, # do you want the top-barplot to have a border?
                                                               axis_param = list(side = axis.side, 
                                                                                 # side = "right", 
                                                                                 #labels = c("zero", "half", "one"),
                                                                                 # at = c(0, 0.5, 1), 
                                                                                 # labels_rot = 45,
                                                                                 gp= gpar(fontsize= barplot.font, fontface="bold"))), 
                           Survival= anno_points(surv.df$survival.time, 
                                                 
                                                 pch = surv.df$pch,
                                                 
                                                 size = unit(5, "mm"),
                                                 
                                                 gp = gpar(col = surv.df$status.col),
                                                 
                                                 axis = TRUE, 
                                                 
                                                 axis_param = list(side = axis.side, 
                                                                   gp= gpar(fontsize= barplot.font, fontface="bold"))),
                           
                           simple_anno_size = unit(1, "cm"), height = unit(top.w, "cm"),
                           annotation_name_gp= gpar(fontsize= legend.title.font, fontface="bold", col="blue"), 
                           annotation_name_offset = c(survival = "0.7cm"),
                           annotation_height =c(20,20), # this controls the height of the response/etc annotation that is added to the columns. However, in order to use mutiple features (e.g., response/celltype/etc) you have to use c(20,20,..) otherwise this generates error
                           gap = unit(c(3,3), "mm")) # this controls the gap between multiple annotation heatbars (for example, the space btw response and patient.id bars) / columns
    
    
    # show_annotation_name= TRUE,
    # annotation_legend_param = list(#title = legend.tit.df,
    #                                 title_gp = gpar(fontsize = legend.title.font, fontface="bold"),
    #                                 title_position = "topcenter", 
    #                                 labels_gp = gpar(fontsize = legend.label.font),
    #                                 grid_height= unit(2, "cm"),
    #                                 grid_width= unit(2, "cm"),
    #                                 nrow= num.rows.annot.lgd,
    #                                 legend_height = unit(2, "cm")
    
    
  }  else {
    
    cat(paste0("\nStart HeatMapAnnotation (line 570)...\n"))
    
    h1 = HeatmapAnnotation(column_bar = anno_oncoprint_barplot(type= NULL,
                                                               border= show.border, # do you want the top-barplot to have a border?
                                                               axis_param = list(side = axis.side, 
                                                                                 # side = "right", 
                                                                                 #labels = c("zero", "half", "one"),
                                                                                 # at = c(0, 0.5, 1), 
                                                                                 # labels_rot = 45,
                                                                                 gp= gpar(fontsize= barplot.font, fontface="bold"))), 
                           
                           simple_anno_size = unit(1, "cm"), height = unit(top.w, "cm"),
                           annotation_name_gp= gpar(fontsize= legend.title.font, fontface="bold", col="blue"), 
                           # annotation_name_offset = c(survival = "0.7cm"),
                           annotation_height = unit(c(20), "mm"), # this controls the height of the response/etc annotation that is added to the columns. However, in order to use mutiple features (e.g., response/celltype/etc) you have to use c(20,20,..) otherwise this generates error
                           gap = unit(c(3), "mm")) # this controls the gap between multiple annotation heatbars (for example, the space btw response and patient.id bars) / columns
    
  }
  
  return(h1)
}