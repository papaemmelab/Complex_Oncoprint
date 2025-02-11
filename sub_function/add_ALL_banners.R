add_ALL_banners <- function(list.my.cols, show.annot.legend, list.ht.colors, lookup.table){
  
  highANY2.col <- list.ht.colors$ALL.SUBTYPE[names(list.ht.colors$ALL.SUBTYPE) %in% unique(lookup.table$HIGH_ANY2_TYPE)]
  
  list.my.cols$HIGH_ANY2_TYPE <- highANY2.col
  
  show.annot.legend <- c(show.annot.legend, "TRUE")
  
  
  
  RNA.col <- list.ht.colors$ALL.SUBTYPE[names(list.ht.colors$ALL.SUBTYPE) %in% unique(lookup.table$RNA_SUBTYPE)]
  
  list.my.cols$RNA_SUBTYPE <- RNA.col
  
  show.annot.legend <- c(show.annot.legend, "TRUE")
  
  
  
  DNA.col <- list.ht.colors$ALL.SUBTYPE[names(list.ht.colors$ALL.SUBTYPE) %in% unique(lookup.table$DNA_SUBTYPE)]
  
  list.my.cols$DNA_SUBTYPE <- DNA.col
  
  show.annot.legend <- c(show.annot.legend, "TRUE")
  
  
  gender.col <- list.ht.colors$GENDER[names(list.ht.colors$GENDER) %in% unique(lookup.table$GENDER)]
  
  list.my.cols$GENDER <- gender.col
  
  show.annot.legend <- c(show.annot.legend, "TRUE")
  
  #=========================
  purity.col <- list.ht.colors$PURITY[names(list.ht.colors$PURITY) %in% unique(lookup.table$PURITY)]
  
  purity.col <- c(purity.col, "white" = "#FFFFFF")
  
  list.my.cols$PURITY <- purity.col
  
  
  show.annot.legend <- c(show.annot.legend, "TRUE")
  
  ee.col <- list.ht.colors$EE[names(list.ht.colors$EE) %in% unique(lookup.table$EE)]
  
  ee.col <- c(ee.col, "white" = "#FFFFFF")
  
  list.my.cols$EE <- ee.col
  
  show.annot.legend <- c(show.annot.legend, "TRUE")
  #=========================
  
  fin.subtype.col <- list.ht.colors$ALL.SUBTYPE[names(list.ht.colors$ALL.SUBTYPE) %in% unique(lookup.table$FINAL_SUBTYPE)]
  
  list.my.cols$FINAL_SUBTYPE <- fin.subtype.col
  
  show.annot.legend <- c(show.annot.legend, "TRUE")
  
  
  soc.col <- list.ht.colors$ALL.SUBTYPE[names(list.ht.colors$ALL.SUBTYPE) %in% unique(lookup.table$SOC_SUBTYPE)]
  
  list.my.cols$SOC_SUBTYPE <- soc.col
  
  show.annot.legend <- c(show.annot.legend, "TRUE")
  
  # extended.diag.col <- list.ht.colors$ALL.SUBTYPE[names(list.ht.colors$ALL.SUBTYPE) %in% unique(lookup.table$EXTENDED_DIAGNOSTIC_SUBTYPE)]
  
  extended.diag.col <- list.ht.colors$ALL.SUBTYPE[names(list.ht.colors$ALL.SUBTYPE) %in% unique(lookup.table$EXT.DX_SUBTYPE)]
  
  # list.my.cols$EXTENDED_DIAGNOSTIC_SUBTYPE <- extended.diag.col
  
  list.my.cols$EXT.DX_SUBTYPE <- extended.diag.col
  
  show.annot.legend <- c(show.annot.legend, "TRUE")
  
  return(list(updated.list.my.cols=list.my.cols, 
              updated.show.annot.legend= show.annot.legend))
}