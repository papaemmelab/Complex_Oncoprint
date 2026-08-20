add_ALL_banners <- function(list.my.cols, show.annot.legend, list.ht.colors, lookup.table){
  
  # match_default_colors() matches values against the palette ignoring case/
  # spacing/punctuation, so e.g. "Bcr-Abl1" still resolves to "BCR-ABL1"'s color.
  
  highANY2.col <- match_default_colors(list.ht.colors$ALL.SUBTYPE, lookup.table$HIGH_ANY2_TYPE)
  
  list.my.cols$HIGH_ANY2_TYPE <- highANY2.col
  
  show.annot.legend <- c(show.annot.legend, "TRUE")
  
  
  
  RNA.col <- match_default_colors(list.ht.colors$ALL.SUBTYPE, lookup.table$RNA_SUBTYPE)
  
  list.my.cols$RNA_SUBTYPE <- RNA.col
  
  show.annot.legend <- c(show.annot.legend, "TRUE")
  
  
  
  DNA.col <- match_default_colors(list.ht.colors$ALL.SUBTYPE, lookup.table$WGS_SUBTYPE)
  
  list.my.cols$DNA_SUBTYPE <- DNA.col
  list.my.cols$WGS_SUBTYPE <- DNA.col
  
  show.annot.legend <- c(show.annot.legend, "TRUE")
  
  
  gender.col <- match_default_colors(list.ht.colors$GENDER, lookup.table$GENDER)
  
  list.my.cols$GENDER <- gender.col
  
  show.annot.legend <- c(show.annot.legend, "TRUE")
  
  #=========================
  purity.col <- match_default_colors(list.ht.colors$PURITY, lookup.table$PURITY)
  
  purity.col <- c(purity.col, "white" = "#FFFFFF")
  
  list.my.cols$PURITY <- purity.col
  
  
  show.annot.legend <- c(show.annot.legend, "TRUE")
  
  ee.col <- match_default_colors(list.ht.colors$EE, lookup.table$EE)
  
  ee.col <- c(ee.col, "white" = "#FFFFFF")
  
  list.my.cols$EE <- ee.col
  
  show.annot.legend <- c(show.annot.legend, "TRUE")
  #=========================
  
  fin.subtype.col <- match_default_colors(list.ht.colors$ALL.SUBTYPE, lookup.table$FINAL_SUBTYPE)
  
  list.my.cols$FINAL_SUBTYPE <- fin.subtype.col
  
  show.annot.legend <- c(show.annot.legend, "TRUE")
  
  
  soc.col <- match_default_colors(list.ht.colors$ALL.SUBTYPE, lookup.table$SOC_SUBTYPE)
  
  list.my.cols$SOC_SUBTYPE <- soc.col
  
  show.annot.legend <- c(show.annot.legend, "TRUE")
  
  extended.diag.col <- match_default_colors(list.ht.colors$ALL.SUBTYPE, lookup.table$EXT.DX_SUBTYPE)
  
  list.my.cols$EXT.DX_SUBTYPE <- extended.diag.col
  
  show.annot.legend <- c(show.annot.legend, "TRUE")
  
  return(list(updated.list.my.cols=list.my.cols, 
              updated.show.annot.legend= show.annot.legend))
}