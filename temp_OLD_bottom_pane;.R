prepare_BOTTOM_annotation <- function(df, ## lookup table
                                      
                                      list.my.cols,legend.title.font,legend.label.font, annot.title.side, 
                                      num.rows.annot.lgd, show.annot.legend, 
                                      ribbon.size, 
                                      
                                      show.ALL= FALSE,
                                      show.MPN= FALSE,
                                      
                                      banner.name= NULL, 
                                      
                                      show.individuals= FALSE,
                                      
                                      banner.label.col = "#5b859e"){
  
  library(stringr)
  
  cat(paste0("\nPrepare bottom annotation ...\n"))
  
  MN <- met.brewer("Monet", type = "discrete")
  RD <- met.brewer("Redon", type = "discrete")
  FH <- met.brewer("Isfahan1", type = "discrete") 
  DM <- met.brewer("Demuth", type = "discrete") 
  HK1 <- met.brewer("Hokusai1", type = "discrete") 
  HK3 <- met.brewer("Hokusai3", type = "discrete") 
  DeR <- met.brewer("Derain", type = "discrete")
  TP <- met.brewer("Tiepolo", type = "discrete")
  LK <- met.brewer("Lakota", type = "discrete")  
  
  # browser()
  
  # names(list.my.cols)[names(list.my.cols) == "new.banner.col"] <- toupper(banner.name[1])
  
  # anno_oncoprint_barplot(type = NULL, which = c("column", "row"),
  #                        bar_width = 0.6, axis = TRUE,
  #                        axis_param = if(which == "column") default_axis_param("column") else list(side = "top", labels_rot = 0),
  #                        width = NULL, height = NULL, border = FALSE)
  
  ##=========================================================================================
  ## Sort the df of bottom annotation according to the original simple.ht sample order   ====
  ##=========================================================================================
  
  banner.name <- toupper(banner.name)
  
  if (show.ALL) {
    banner.name <- c(banner.name, c("FINAL_SUBTYPE", "CNV.WGS.CNVKIT.RHO", "RNA.EE"))
    df <- df %>% dplyr::select(all_of(unique(c("TARGET_NAME", (banner.name)))))
  } else if (show.MPN) {
    banner.name <- c(banner.name, c("Complex.Karyotype"))
    df <- df %>% dplyr::select(all_of(unique(c("TARGET_NAME", (banner.name)))))
  } else {
    df <- df %>% dplyr::select(all_of(c("TARGET_NAME", (banner.name))))
  }
  
  banner.name <- unique(banner.name)
  
  rownames(df) <- df$TARGET_NAME
  # df<-df[new.column_order,]
  rownames(df) <- NULL
  df$TARGET_NAME <- NULL
  
  colnames(df)[colnames(df) == 'INDIVIDUAL.ID'] <- 'Patient.ID'
  names(list.my.cols)[names(list.my.cols) == 'INDIVIDUAL.ID'] <- 'Patient.ID'
  
  if (!(show.individuals)){
    df$INDIVIDUAL.ID <- NULL
  }
  
  # colnames(df) <- str_to_title(colnames(df)) 
  # names(list.my.cols) <- str_to_title(names(list.my.cols))
  # if (!is.null(response.order)){
  #   df$RESPONSE <- factor(df$RESPONSE, levels= response.order)
  # }   
  
  # response.title.pos <- annot.title.side
  
  ################################################################################
  #### Change the UPPSER-CASE extra-annot titles to capitalize the 1st letter ====
  ################################################################################
  
  EXCEPT <- c("Patient.ID", "FINAL_SUBTYPE", "CNV.WGS.CNVKIT.RHO", "RNA.EE", "Complex.Karyotype")
  
  banner.name <- ifelse(banner.name %in% EXCEPT, banner.name, str_to_title(banner.name))
  names(list.my.cols) <- ifelse(names(list.my.cols) %in% EXCEPT,
                                names(list.my.cols),
                                str_to_title(names(list.my.cols)))
  
  colnames(df) <- ifelse(colnames(df) %in% EXCEPT, colnames(df), str_to_title(colnames(df)))
  
  ########################################################
  #### Define Bottom Annotation obj (e.g., RESPONSE) ====
  ########################################################
  
  if (show.MPN){
    
    df$Complex.Karyotype <- factor(df$Complex.Karyotype, 
                                   levels = c("complex", "not complex", "not available"))
    
    list.my.cols$Complex.Karyotype <- c(
      "complex" = "firebrick",
      "not complex" = MN[3],
      "not available" ="lightgray"
      
    )
    
    names(list.my.cols)[names(list.my.cols) == 'INDIVIDUAL.ID'] <- 'Patient.ID'
    
    h2 = HeatmapAnnotation(df = df %>% dplyr::select(all_of(banner.name)), name= "TEST", #df = data.frame(PATIENTS = pts), col= list(PATIENTS = col.assign), 
                           col = list.my.cols,
                           na_col = "darkgrey",
                           simple_anno_size = unit(ribbon.size, "cm"), # size of the ribbon
                           annotation_height =c(20,20), # this controls the height of the response/etc annotation that is added to the columns. However, in order to use mutiple features (e.g., response/celltype/etc) you have to use c(20,20,..) otherwise this generates error
                           # gap = unit(c(5,5), "mm"), # this controls the gap between multiple annotation heatbars (for example, the space btw response and patient.id bars)
                           gap = unit(rep(5,ncol(df)),"mm"),
                           show_annotation_name= rep(FALSE,ncol(df)), ################# CHANGD THIS to FALSE for MPN
                           show_legend = as.logical(show.annot.legend),
                           #show_annotation_name= show.annot.legend, 
                           annotation_name_offset = unit(20, "mm"),
                           gp = gpar(col = "black"),
                           annotation_name_gp= gpar(fontsize = legend.title.font, fontface= "bold", col= banner.label.col), # blue
                           annotation_legend_param = list(#title = legend.tit.df,
                             title_gp = gpar(fontsize = legend.title.font, fontface="bold"),
                             title_position = "leftcenter",  ##### annot.title.side
                             
                             labels_gp = gpar(fontsize = legend.label.font),
                             grid_height= unit(1, "cm"), # size of the box in the legends
                             grid_width= unit(1, "cm"),
                             nrow= num.rows.annot.lgd,
                             legend_height = unit(20, "cm"),
                             annotation_name_align = TRUE,
                             
                             annotation_legend_param = list(title = "Karyotype Complexity",
                                                            Karyotype.Complexity = list(title = "Karyotype Complexity", title_gp = gpar(fontsize = 12)), # Modify legend for "Group" title
                                                            title_gp = gpar(fontsize = legend.title.font, fontface="bold"),
                                                            title_position = "left-center", 
                                                            
                                                            labels_gp = gpar(fontsize = legend.label.font),
                                                            grid_height= unit(1, "cm"), # size of the box in the legends
                                                            grid_width= unit(1, "cm"),
                                                            nrow= num.rows.annot.lgd,
                                                            legend_height = unit(5, "cm"),
                                                            legend_direction = "horizontal")
                           ))
    
  } else if (show.ALL) {
    
    banner.labels <- banner.name
    banner.labels <- gsub("CNV.WGS.CNVKIT.RHO", "WGS.RHO", banner.labels, ignore.case = TRUE)
    
    # browser()
    
    col_fun = colorRamp2(
      c(0, 60, 100),                 # shifted midpoint from 50 → 60
      c(DM[9], MN[4], "#cc6c4a")     # softened darkest color (previously "#af4f2f")
    ) # DM[9] is unclassified color
    
    
    list.my.cols$CNV.WGS.CNVKIT.RHO <- col_fun
    list.my.cols$RNA.EE <- col_fun
    
    # colnames(df) <- toupper(colnames(df))
    # browser()
    
    h2 = HeatmapAnnotation(df = df %>% dplyr::select(all_of(banner.name)) , name= "TEST", #df = data.frame(PATIENTS = pts), col= list(PATIENTS = col.assign), 
                           col = list.my.cols,
                           na_col = "white", # Set missing values to white
                           simple_anno_size = unit(ribbon.size, "cm"), # size of the ribbon
                           annotation_height =c(20,20), # this controls the height of the response/etc annotation that is added to the columns. However, in order to use mutiple features (e.g., response/celltype/etc) you have to use c(20,20,..) otherwise this generates error
                           gap = unit(rep(5,ncol(df)),"mm"), # this controls the gap between multiple annotation heatbars (for example, the space btw response and patient.id bars)
                           
                           show_annotation_name= rep(TRUE,ncol(df)),
                           annotation_label = banner.labels,
                           
                           # show_legend = as.logical(show.annot.legend),
                           show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                           
                           annotation_name_offset = unit(20, "mm"),
                           
                           # gp = gpar(col = "black"), ## this adds a line to the annotaton blocks not nice for heatmap of purity
                           
                           annotation_name_gp= gpar(fontsize = legend.title.font, fontface= "bold", col= banner.label.col),
                           
                           annotation_legend_param = list(WGS.RHO.RNA = list(title = "HOOOOPOOOO", title_gp = gpar(fontsize = 12)), # Modify legend for "Group" title CNV.WGS.CNVKIT.RHO
                                                          
                                                          # CNV.WGS.CNVKIT.RHO = list(title = "WGS.Purity/RNA.EE"),
                                                          title_gp = gpar(fontsize = legend.title.font, fontface="bold"),
                                                          title_position = annot.title.side, 
                                                          
                                                          labels_gp = gpar(fontsize = legend.label.font),
                                                          grid_height= unit(1, "cm"), # size of the box in the legends
                                                          grid_width= unit(1, "cm"),
                                                          nrow= num.rows.annot.lgd,
                                                          legend_height = unit(5, "cm"),
                                                          legend_direction = "horizontal"
                           )
                           
    )
    
  } else {
    
    # test colors if not shown properly
    #----------------------------------------------
    # source("./sub_function/color_alpha_test.R")
    # color_alpha_test(list.my.cols$RESPONSE)
    #----------------------------------------------
    # browser()
    
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
                           
                           annotation_name_offset = unit(20, "mm"), ### KEEP the label of annotation bar with a bit offset
                           
                           gp = gpar(col = "black"),
                           
                           annotation_name_gp= gpar(fontsize = legend.title.font, fontface= "bold", col= banner.label.col), #"blue"
                           
                           annotation_legend_param = list(#title = legend.tit.df,
                             title_gp = gpar(fontsize = legend.title.font, fontface="bold"),
                             title_position = annot.title.side, 
                             labels_gp = gpar(fontsize = legend.label.font),
                             grid_height= unit(1, "cm"), # size of the box in the legends
                             grid_width= unit(1, "cm"),
                             nrow= num.rows.annot.lgd,
                             legend_height = unit(20, "cm")
                           )
    )
  }
  
  
  
  return(list(h2=h2,
              df.updated=df,
              list.my.cols.updated= list.my.cols))
  
}