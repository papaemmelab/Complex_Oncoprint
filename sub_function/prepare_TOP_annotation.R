prepare_TOP_annotation <- function(list.colors,show.border,  
                                   
                                   barplot.font, 
                                   legend.title.font, 
                                   
                                   top.w, 
                                   
                                   show.ALL= FALSE, 
                                   purity.df= NULL, 
                                   
                                   show.survival=FALSE, 
                                   show.blast= FALSE,
                                   show.MPN= FALSE, 
                                   
                                   axis.side = "left",
                                   top.annot.axis.side = "left",  ## for now keep it left, u can adapt it in future
                                   top.annotation_name_side = "left",
                                   
                                   banner.label.col = "#5b859e",
                                   
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
  
  bad.col <-  "red"
  good.col <- MN[1] 
  
  legend.label.font <- barplot.font  #### <<< fix this in future if u want to control the legend font for top plot
  
  ####################################
  ####################################

  #== ALL options ====
  
  show.MPN.historic = FALSE 
  
  if (show.MPN.historic){
    
    # browser()
    col_fun = colorRamp2(c(0, 50, 100), c("blue", "white", "#af4f2f"))
    
    h1 = HeatmapAnnotation(column_bar = anno_oncoprint_barplot(type= NULL,
                                                               border= show.border, # do you want the top-barplot to have a border?
                                                               axis_param = list(side = axis.side, 
                                                                                 # side = "right", 
                                                                                 #labels = c("zero", "half", "one"),
                                                                                 # at = c(0, 0.5, 1), 
                                                                                 # labels_rot = 45,
                                                                                 gp= gpar(fontsize= barplot.font, fontface="bold"))), 
                           
                           BM_Blast_Baseline = anno_points(lookup.table$BM_BLASTS_BASELINE,  ylim = c(0, max(lookup.table$BM_BLASTS_BASELINE, na.rm = TRUE)+5),
                                                           size = unit(5,"mm"),
                                                           width = unit(1, "cm"),
                                                           gp = gpar(col = ifelse(lookup.table$BM_BLASTS_BASELINE < 20, bad.col, FH[1])), 
                                                           height = unit(1, "cm"),
                                                           axis_param = list(side = top.annot.axis.side, 
                                                                             gp= gpar(fontsize= barplot.font, fontface="bold"))),
                                                           
                           # BM_Blasts = anno_points(lookup.table %>% dplyr::select(BM_BLASTS_BASELINE),
                           #                         axis_param = list(side = axis.side, 
                           #                                           #gp= gpar(fontsize= barplot.font, fontface="bold")),
                           # 
                           #                         border = TRUE),
                           
                           
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
                           
                           annotation_name_gp= gpar(fontsize= legend.title.font, fontface="bold", col= banner.label.col), #blue
                           annotation_name_side= top.annotation_name_side, ## brings TBA.blast label to left
                           annotation_name_rot= 0,
                           annotation_name_offset = unit(20, "mm"),
                           
                           
                           show_legend = c(TRUE, FALSE, FALSE),
                           
                           # annotation_name_offset = c(survival = "0.7cm"),
                           
                           annotation_height = unit(c(20), "mm"), # this controls the height of the response/etc annotation that is added to the columns. However, in order to use mutiple features (e.g., response/celltype/etc) you have to use c(20,20,..) otherwise this generates error
                           
                           gap = unit(c(3), "mm")) # this controls the gap between multiple annotation heatbars (for example, the space btw response and patient.id bars) / columns
    
    #== MPN options ====
    
  } else if (show.MPN){
    
    browser()
    
    blast_col <- c("<20"=DM[6], ">=20"=FH[1])

    # Update the annotations with the shared ylim
    
    h1 <- HeatmapAnnotation(column_bar = anno_oncoprint_barplot(type= NULL,
                                                                border= show.border, # do you want the top-barplot to have a border?
                                                                axis_param = list(side = axis.side, 
                                                                                  # side = "right", 
                                                                                  #labels = c("zero", "half", "one"),
                                                                                  # at = c(0, 0.5, 1), 
                                                                                  # labels_rot = 45,
                                                                                  gp= gpar(fontsize= barplot.font, fontface="bold"))), 
                            
                            BLAST = anno_points(
                              lookup.table$BLAST,
                              ylim  = c(0, max(lookup.table$BLAST, na.rm = TRUE) + 5),
                              size  = unit(3, "mm"),
                              width = unit(1, "cm"),
                              height = unit(1, "cm"),
                              pch   = 16,
                              gp    = gpar(col = ifelse(lookup.table$BLAST >= 20, blast_col[2], blast_col[1])),
                              axis_param = list(
                                side = top.annot.axis.side,
                                gp   = gpar(fontsize = barplot.font, fontface = "bold")
                              ),
                              legend_param = list(   ### currently not working and am setting legend outside directly
                                title  = "Blast",
                                at     = c("<20", ">=20"),
                                labels = c("<20", ">=20"),
                                type   = "points",
                                pch    = 5,
                                legend_gp = gpar(col = c(blast_col[1], blast_col[2]))
                              )
                            ),
                            
                            simple_anno_size = unit(5, "cm"), height = unit(top.w, "cm"), ## not affecting annot point size, idk
                            
                            annotation_name_gp= gpar(fontsize= legend.title.font, fontface="bold", col= banner.label.col), # blue top annot names
                            
                            annotation_name_side= top.annotation_name_side, ## brings TBA.blast label to left
                            
                            annotation_name_rot= 0,
                            
                            annotation_name_offset = unit(3, "cm"),
                            
                            show_legend = c(TRUE, TRUE),
                            
                            # annotation_name_offset = c(survival = "0.7cm"),
                            
                            annotation_height = unit(c(15), "mm"), # this controls the height of the response/etc annotation that is added to the columns. However, in order to use mutiple features (e.g., response/celltype/etc) you have to use c(20,20,..) otherwise this generates error
                            
                            gap = unit(c(5,5), "mm")) 
    
    
  } else if (show.blast){
    
    # browser()

    blast_col <- c("≤10"=DM[6], ">10"=FH[1])
    
    # Update the annotations with the shared ylim
    
    h1 <- HeatmapAnnotation(column_bar = anno_oncoprint_barplot(type= NULL,
                                                                        border= show.border, # do you want the top-barplot to have a border?
                                                                        axis_param = list(side = axis.side, 
                                                                                          # side = "right", 
                                                                                          #labels = c("zero", "half", "one"),
                                                                                          # at = c(0, 0.5, 1), 
                                                                                          # labels_rot = 45,
                                                                                          gp= gpar(fontsize= barplot.font, fontface="bold"))), 
      BM.Blast = anno_points(
        lookup.table$BM...BLASTS,
        ylim  = c(0, max(lookup.table$BM...BLASTS, na.rm = TRUE) + 5),
        size  = unit(5, "mm"),
        width = unit(1, "cm"),
        height = unit(1, "cm"),
        pch   = 16,
        gp    = gpar(col = ifelse(lookup.table$BM...BLASTS > 10,blast_col[1], blast_col[2])),
        axis_param = list(
          side = top.annot.axis.side,
          gp   = gpar(fontsize = barplot.font, fontface = "bold")
        ),
        legend_param = list(
          title  = "BM Blast",
          at     = c("≤10", ">10"),
          labels = c("≤10", ">10"),
          type   = "points",
          pch    = 16,
          legend_gp = gpar(col = c(FH[1], DM[6]))
        )
      ),
                            
                            
                          
      
      Heme.Blast = anno_points(lookup.table$HEME...BLASTS, ylim = c(0, max(lookup.table$HEME...BLASTS, na.rm = TRUE)+5),
                                size = unit(5,"mm"),
                                
                              width = unit(1, "cm"),
                              gp = gpar(col = ifelse(lookup.table$HEME...BLASTS > 10, FH[1], DM[6])),
                              height = unit(1, "cm"),
                              axis_param = list(side = top.annot.axis.side, 
                                                gp= gpar(fontsize= barplot.font, fontface="bold"))),
      
      
      simple_anno_size = unit(5, "cm"), height = unit(top.w, "cm"), ## not affecting annot point size, idk
      
      annotation_name_gp= gpar(fontsize= legend.title.font, fontface="bold", col=banner.label.col), # blue top annot names
      annotation_name_side= top.annotation_name_side, ## brings TBA.blast label to left
      annotation_name_rot= 0,
      
      annotation_name_offset = unit(20, "mm"),
      
      show_legend = TRUE,
      
      # annotation_name_offset = c(survival = "0.7cm"),
      
      annotation_height = unit(c(20), "mm"), # this controls the height of the response/etc annotation that is added to the columns. However, in order to use mutiple features (e.g., response/celltype/etc) you have to use c(20,20,..) otherwise this generates error
      
      gap = unit(c(5,5), "mm")) 
    
    #== Add survival (need DEV) ================
    
  } else if (show.survival) {
    
    stop("\n DEVELOP show.survival in TOP annot")
    
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
                           
                           simple_anno_size = unit(1, "cm"), 
                           height = unit(top.w, "cm"),
                           annotation_name_gp= gpar(fontsize= legend.title.font, fontface="bold", col= banner.label.col), #blue
                           
                           annotation_name_offset = unit(20, "mm"),
                           
                           # annotation_name_offset = c(survival = "0.7cm"),
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
    #== Default TOP annot ================
    
  }  else {
      h1 = HeatmapAnnotation(column_bar = anno_oncoprint_barplot(type= NULL,
                                                                 border= show.border, # do you want the top-barplot to have a border?
                                                                 axis_param = list(side = axis.side, 
                                                                                   # side = "right", 
                                                                                   #labels = c("zero", "half", "one"),
                                                                                   # at = c(0, 0.5, 1), 
                                                                                   # labels_rot = 45,
                                                                                   gp= gpar(fontsize= barplot.font, fontface="bold"))), 
                             
                             simple_anno_size = unit(1, "cm"), height = unit(top.w, "cm"),
                             annotation_name_gp= gpar(fontsize= legend.title.font, fontface="bold", col= banner.label.col), #blue 
                             annotation_name_offset = unit(20, "mm"),
                             show_legend = TRUE,
                             
                             # annotation_name_offset = c(survival = "0.7cm"),
                             annotation_height = unit(c(20), "mm"), # this controls the height of the response/etc annotation that is added to the columns. However, in order to use mutiple features (e.g., response/celltype/etc) you have to use c(20,20,..) otherwise this generates error
                             gap = unit(c(3), "mm")) # this controls the gap between multiple annotation heatbars (for example, the space btw response and patient.id bars) / columns
      
    }
   
  return(h1)
  
}

