prepare_TOP_annotation <- function(list.colors,show.border,axis.side,barplot.font, legend.title.font, top.w, 
                                   show.purity= FALSE, purity.df= NULL, show.survival=FALSE, show.blast= FALSE,show.latest.blast= FALSE, 
                                   lookup.table= lookup.table,  ...){
  
  library(circlize)
  library(viridis)
  library(MetBrewer)
  
  MN <- met.brewer("Monet", type = "discrete")
  RD <- met.brewer("Redon", type = "discrete")
  FH <- met.brewer("Isfahan1", type = "discrete") 
  DM <- met.brewer("Demuth", type = "discrete") 
  HK1 <- met.brewer("Hokusai1", type = "discrete") 
  HK3 <- met.brewer("Hokusai3", type = "discrete") 
  DeR <- met.brewer("Derain", type = "discrete")
  TP <- met.brewer("Tiepolo", type = "discrete")
  LK <- met.brewer("Lakota", type = "discrete")
  
  ####################################
  ####################################
  
  if (show.blast){
    
    # Update the annotations with the shared ylim
    h1 <- HeatmapAnnotation(column_bar = anno_oncoprint_barplot(type= NULL,
                                                                        border= show.border, # do you want the top-barplot to have a border?
                                                                        axis_param = list(side = axis.side, 
                                                                                          # side = "right", 
                                                                                          #labels = c("zero", "half", "one"),
                                                                                          # at = c(0, 0.5, 1), 
                                                                                          # labels_rot = 45,
                                                                                          gp= gpar(fontsize= barplot.font, fontface="bold"))), 
      BM.BLASTS = anno_points(lookup.table$BM...BLASTS,  ylim = c(0, max(lookup.table$BM...BLASTS, na.rm = TRUE)+5),
                              size = unit(5,"mm"),
                        width = unit(1, "cm"),
                        gp = gpar(col = ifelse(lookup.table$BM...BLASTS > 10, RD[1], RD[6])),
                        height = unit(1, "cm")),
      
      Heme.BLASTS = anno_points(lookup.table$HEME...BLASTS, ylim = c(0, max(lookup.table$HEME...BLASTS, na.rm = TRUE)+5),
                                size = unit(5,"mm"),
                                
                              width = unit(1, "cm"),
                              gp = gpar(col = ifelse(lookup.table$BM...BLASTS > 10, RD[1], RD[6])),
                              height = unit(1, "cm")),
      
      
      simple_anno_size = unit(5, "cm"), height = unit(top.w, "cm"), ## not affecting annot point size, idk
      
      annotation_name_gp= gpar(fontsize= legend.title.font, fontface="bold", col="blue"), 
      
      show_legend = c(TRUE, FALSE, FALSE),
      
      # annotation_name_offset = c(survival = "0.7cm"),
      
      annotation_height = unit(c(20), "mm"), # this controls the height of the response/etc annotation that is added to the columns. However, in order to use mutiple features (e.g., response/celltype/etc) you have to use c(20,20,..) otherwise this generates error
      
      gap = unit(c(3), "mm")) 
  
  } else if (show.latest.blast){
      
      # Update the annotations with the shared ylim
      h1 <- HeatmapAnnotation(column_bar = anno_oncoprint_barplot(type= NULL,
                                                                  border= show.border, # do you want the top-barplot to have a border?
                                                                  axis_param = list(side = axis.side, 
                                                                                    # side = "right", 
                                                                                    #labels = c("zero", "half", "one"),
                                                                                    # at = c(0, 0.5, 1), 
                                                                                    # labels_rot = 45,
                                                                                    gp= gpar(fontsize= barplot.font, fontface="bold"))), 
                              BLAST = anno_points(lookup.table$LATEST.BLAST,  ylim = c(0, max(lookup.table$LATEST.BLAST, na.rm = TRUE)+5),
                                                      size = unit(5,"mm"),
                                                      width = unit(1, "cm"),
                                                      gp = gpar(col = ifelse(lookup.table$LATEST.BLAST > 20, RD[1], RD[6])),
                                                      height = unit(1, "cm"),
                                                  axis_param = list(side = axis.side, 
                                                                    gp= gpar(fontsize= barplot.font, fontface="bold"))),

                              
                              simple_anno_size = unit(5, "cm"), height = unit(top.w, "cm"), ## not affecting annot point size, idk
                              
                              annotation_name_gp= gpar(fontsize= legend.title.font, fontface="bold", col="blue"), 
                              
                              show_legend = c(TRUE, FALSE),
                              
                              # annotation_name_offset = c(survival = "0.7cm"),
                              
                              annotation_height = unit(c(20), "mm"), # this controls the height of the response/etc annotation that is added to the columns. However, in order to use mutiple features (e.g., response/celltype/etc) you have to use c(20,20,..) otherwise this generates error
                              
                              gap = unit(c(3), "mm")) 
      
    }
  
  ####################################
  ####################################
  
  else if (show.survival) {
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
    
    
  } else if (show.purity){

      col_fun = colorRamp2(c(0, 50, 100), c("blue", "white", "#af4f2f"))

      h1 = HeatmapAnnotation(column_bar = anno_oncoprint_barplot(type= NULL,
                                                                 border= show.border, # do you want the top-barplot to have a border?
                                                                 axis_param = list(side = axis.side, 
                                                                                   # side = "right", 
                                                                                   #labels = c("zero", "half", "one"),
                                                                                   # at = c(0, 0.5, 1), 
                                                                                   # labels_rot = 45,
                                                                                   gp= gpar(fontsize= barplot.font, fontface="bold"))), 
                             TBA.Survival = anno_points(purity.df %>% dplyr::select(CNV.WGS.CNVKIT.RHO),
                                                        axis_param = list(side = axis.side, 
                                                                          gp= gpar(fontsize= barplot.font, fontface="bold"))),
                             
                             # CNVkit.Purity  = purity.df %>% pull(CNV.WGS.CNVKIT.RHO),  
                             # 
                             # RNA.EE  = purity.df %>% pull(RNA.EE),  
                             # 
                             # col = list(CNVkit.Purity = col_fun,
                             #            RNA.EE= col_fun), 
                             # 
                             # annotation_legend_param = list(CNVkit.Purity = list(title = "Purity/EE",
                             #                                                     labels_gp = gpar(col = "black", fontsize = barplot.font),
                             #                                                     title_gp = gpar(col = "black", fontsize = barplot.font, fontface="bold"))
                             #                                ),
                             

                             
                             simple_anno_size = unit(1, "cm"), height = unit(top.w, "cm"),
                             
                             annotation_name_gp= gpar(fontsize= legend.title.font, fontface="bold", col="blue"), 
                             
                             show_legend = c(TRUE, FALSE, FALSE),
                             
                             # annotation_name_offset = c(survival = "0.7cm"),
                             
                             annotation_height = unit(c(20), "mm"), # this controls the height of the response/etc annotation that is added to the columns. However, in order to use mutiple features (e.g., response/celltype/etc) you have to use c(20,20,..) otherwise this generates error
                             
                             gap = unit(c(3), "mm")) # this controls the gap between multiple annotation heatbars (for example, the space btw response and patient.id bars) / columns
    } else {
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