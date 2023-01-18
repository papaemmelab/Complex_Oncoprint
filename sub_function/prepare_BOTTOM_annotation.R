prepare_BOTTOM_annotation <- function(df, list.my.cols,legend.title.font,legend.label.font, annot.title.side, num.rows.annot.lgd, show.annot.legend, ribbon.size, show.individuals= FALSE){
  
    
    cat(paste0("\nPrepare bottom annotation ...\n"))
    
    # names(list.my.cols)[names(list.my.cols) == "new.banner.col"] <- toupper(banner.name[1])
    
    # anno_oncoprint_barplot(type = NULL, which = c("column", "row"),
    #                        bar_width = 0.6, axis = TRUE,
    #                        axis_param = if(which == "column") default_axis_param("column") else list(side = "top", labels_rot = 0),
    #                        width = NULL, height = NULL, border = FALSE)
    
    ##=========================================================================================
    ## Sort the df of bottom annotation according to the original simple.ht sample order   ====
    ##=========================================================================================
    
    rownames(df) <- df$TARGET_NAME
    # df<-df[new.column_order,]
    rownames(df) <- NULL
    df$TARGET_NAME <- NULL
    
    
    # if (!is.null(response.order)){
    #   df$RESPONSE <- factor(df$RESPONSE, levels= response.order)
    # }   
    
    # response.title.pos <- annot.title.side
    
    ########################################################
    #### Define Bottom Annotation obj (e.g., RESPONSE) ====
    ########################################################
    
    if (!(show.individuals)){
      df$INDIVIDUAL.ID <- NULL
    }
    
    colnames(df)[colnames(df) == 'INDIVIDUAL.ID'] <- 'Patient.ID'
    names(list.my.cols)[names(list.my.cols) == 'INDIVIDUAL.ID'] <- 'Patient.ID'
    
    colnames(df) <- str_to_title(colnames(df)) 
    names(list.my.cols) <- str_to_title(names(list.my.cols))
    
    h2 = HeatmapAnnotation(df = df , name= "TEST", #df = data.frame(PATIENTS = pts), col= list(PATIENTS = col.assign), 
                           col = list.my.cols,
                           na_col = "grey",
                           simple_anno_size = unit(ribbon.size, "cm"), # size of the ribbon
                           annotation_height =c(20,20), # this controls the height of the response/etc annotation that is added to the columns. However, in order to use mutiple features (e.g., response/celltype/etc) you have to use c(20,20,..) otherwise this generates error
                           # gap = unit(c(5,5), "mm"), # this controls the gap between multiple annotation heatbars (for example, the space btw response and patient.id bars)
                           gap = unit(rep(5,ncol(df)),"mm"),
                           show_annotation_name= rep(TRUE,ncol(df)),
                           show_legend = as.logical(show.annot.legend),
                           #show_annotation_name= show.annot.legend, 
                           annotation_name_offset = unit(20, "mm"),
                           gp = gpar(col = "black"),
                           annotation_name_gp= gpar(fontsize = legend.title.font, fontface= "bold", col="blue"),
                           annotation_legend_param = list(#title = legend.tit.df,
                             title_gp = gpar(fontsize = legend.title.font, fontface="bold"),
                             # title_position = annot.title.side, 
                             title_position = annot.title.side, 
                             
                             labels_gp = gpar(fontsize = legend.label.font),
                             grid_height= unit(1, "cm"), # size of the box in the legends
                             grid_width= unit(1, "cm"),
                             nrow= num.rows.annot.lgd,
                             legend_height = unit(20, "cm")
                           )
    )

  return(list(h2=h2,
              df.updated=df,
              list.my.cols.updated= list.my.cols))
  
}