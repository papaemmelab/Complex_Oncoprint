prepare_COMPLEX_aes <- function(data, M, highlight.events, df, list.my.cols, 
                                show.multis, show.another.banner, 
                                show.response = FALSE, # kept for backward-compat with generate_complex_oncoprint_MRD.R/_paired.R; unused here, RESPONSE is now a generic banner
                                show.individuals= FALSE, show.individuals.legend= FALSE,
                                legend.title.font, legend.label.font, 
                                annot.title.side, 
                                num.rows.annot.lgd, show.annot.legend, 
                                ribbon.size= 1, 
                                banner.name, 
                                rows.font,
                                split.cols.by,
                                show.ALL= FALSE,
                                show.MPN= FALSE,
                                banner.label.col = "#5b859e",
                                legend.height = NULL,
                                banner.case.exceptions = NULL, # default exceptions set in prepare_BOTTOM_annotation : c("Patient.ID", "MRD_subtype", "CNV.WGS.CNVKIT.RHO", "RNA.EE", "Complex.Karyotype")
                                na_col = NULL) {
  
  MN <- met.brewer("Monet", type = "discrete")
  RD <- met.brewer("Redon", type = "discrete")
  FH <- met.brewer("Isfahan1", type = "discrete") 
  DM <- met.brewer("Demuth", type = "discrete") 
  HK1 <- met.brewer("Hokusai1", type = "discrete") 
  HK3 <- met.brewer("Hokusai3", type = "discrete") 
  DeR <- met.brewer("Derain", type = "discrete")
  TP <- met.brewer("Tiepolo", type = "discrete")
  LK <- met.brewer("Lakota", type = "discrete")
  CAS1 <- met.brewer("Cassatt1", type = "discrete")
  CAS2 <- met.brewer("Cassatt2", type = "discrete")
  
  # browser()
  
  #############################
  #### Tag multis (dots) ====
  #############################
  
  if (show.multis){
    cat(paste0("\nStart multi.hit Oncoprint preparation...\n"))
    
    multi.hits <- data %>% dplyr::group_by(TARGET_NAME, GENE) %>% dplyr::mutate(N= n()) %>% dplyr::filter(N>1) %>% dplyr::select(TARGET_NAME, GENE) %>% unique()
    
    multi.hits <- data.frame(multi.hits)
    
    if (nrow(multi.hits)>0){
      for (k in 1: nrow(multi.hits)){
        M[as.character(multi.hits$GENE[k]), as.character(multi.hits$TARGET_NAME[k])] <- paste0(M[as.character(multi.hits$GENE[k]), as.character(multi.hits$TARGET_NAME[k])], "multi_hit",";", collapse = "")
      }
    }
  }
  
    #############################
    #### Prepare banner list ====
    #############################
  
    # Remove INDIVIDUAL.ID from lookup (df) before adding banners
    
    # if ((show.another.banner) | (show.response) | (show.individuals) ) {
    if ((show.another.banner) | (show.individuals) ) {
      
      if (!show.individuals){
        df$INDIVIDUAL.ID <- NULL
      }
    
    #############################
    #### BOTTOM ANNOTATION ====
    #############################
    
    source(file.path("./sub_function/prepare_BOTTOM_annotation.R"))
    BotAnnot <- prepare_BOTTOM_annotation(df, list.my.cols,
                                          legend.title.font,legend.label.font,
                                          annot.title.side, num.rows.annot.lgd, show.annot.legend, 
                                          ribbon.size= ribbon.size, 
                                          banner.name= banner.name, 
                                          show.individuals= show.individuals,
                                          show.ALL= show.ALL,
                                          show.MPN= show.MPN,
                                          banner.label.col= banner.label.col,
                                          banner.case.exceptions= banner.case.exceptions,
                                          legend.height = legend.height,
                                          na_col = na_col
                                          )
    
    #### Get bottom.Annot features ====
    
    h2 <- BotAnnot$h2
    df <- BotAnnot$df.updated
    list.my.cols <- BotAnnot$list.my.cols.updated
    
  } else  {
    h2 = NULL # for example, you do not have any added bottom annotation but still like to see multis
    
  }

  ###############################################################
  # == Add split cols (needs UPDATE) ==== 
  ##############################################################

  if (!is.null(split.cols.by)){
    
    split.cols.by = toupper(split.cols.by)
    
    # a complicated select based on dynamic col-name that is passed in "split.cols.by"
    # first select the dynamic col from lookup and then choose the order based on the sample-names in M.
    # The final class must be numeric for proper depiction
    
    split.cols.order <-  as.numeric(as.factor(lookup.table[[split.cols.by]][match(colnames(M), lookup.table$TARGET_NAME)]))
    
    # split.cols.order <- as.numeric(lookup.table$RESPONSE.ELN.R1[match(colnames(M), lookup.table$TARGET_NAME)])
  } else {
    split.cols.order <- NULL
  }
  
  
  # split.cols.order <- as.numeric(lookup.table$RESPONSE.ELN.R1[match(colnames(M), lookup.table$TARGET_NAME)])
  
  ###############################################################
  # == Add Highlights ==== 
  ##############################################################
  
  # Rows to highlight
  if (!is.null(highlight.events)){
    
    myRows <- intersect(highlight.events, rownames(M))
    
    # Set stylings for row names and make our selected rows unique
    row_idx <- which(rownames(M) %in% myRows)
    fontsizes <- rep(rows.font, nrow(M))
    fontfaces <- rep("bold", nrow(M))
    fontcolors <- rep("black", nrow(M))
    
    fontsizes[row_idx] <- rows.font+2
    # fontcolors[row_idx] <- "#175f5d" #FH[7] # B-ALL
    fontcolors[row_idx] <- "#7b2c36" # MPN EHA
    
    # Set up fill colors for rows
    fill.colors <- rep("white", nrow(M)) # Default color
    # fill.colors[row_idx] <- "#E8F5E9"        # Highlight color # B-ALL
    fill.colors[row_idx] <- "#fff1cc"
    
  } else {
    fontsizes <- rows.font
    fontfaces <- "bold"
    fontcolors <- "black"
    fill.colors <- "white"
  }
  
  # Create text annotation object for displaying row names
  #==========================================================
  
  rowAnno <- rowAnnotation(rows = anno_text(rownames(M), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors, fill= fill.colors)))
  
  # col_hclust = hclust(dist(matrix(rnorm(nrow(M)*ncol(M)), ncol(M))))
  
  #=====================================
  # Generate MPN-specific aethetics
  #=====================================
  
  if (show.MPN){
    names(h2) <- "Complex karyotype Status"
  }
  
  return(list(BotAnnot= h2,
              rowAnno= rowAnno,
              split.cols.order= split.cols.order
              ))
}
