prepare_BOTTOM_annotation <- function(df, ## lookup table
                                      
                                      list.my.cols,legend.title.font,legend.label.font, annot.title.side, 
                                      num.rows.annot.lgd, show.annot.legend, 
                                      ribbon.size, 
                                      
                                      show.ALL= FALSE,
                                      show.MPN= FALSE,

                                      banner.name= NULL, 
                                      
                                      show.individuals= FALSE,
                                      show.individuals.legend = FALSE,
                                      
                                      legend.height = NULL,
                                      
                                      # legend.direction = NULL,  ## currently sticking to default

                                      na_col = NULL,
                                      
                                      banner.label.col = "#5b859e"){
  
    library(stringr)
  
    cat(paste0("\nPrepare bottom annotation ...\n"))
      
      # palettes you use elsewhere (unchanged)
      MN  <- met.brewer("Monet", type = "discrete")
      RD  <- met.brewer("Redon", type = "discrete")
      FH  <- met.brewer("Isfahan1", type = "discrete") 
      DM  <- met.brewer("Demuth", type = "discrete") 
      HK1 <- met.brewer("Hokusai1", type = "discrete") 
      HK3 <- met.brewer("Hokusai3", type = "discrete") 
      DeR <- met.brewer("Derain", type = "discrete")
      TP  <- met.brewer("Tiepolo", type = "discrete")
      LK  <- met.brewer("Lakota", type = "discrete")  
      
      ################################################################################
      ## Prepare banner.name/df with requested extra banners -----
      ################################################################################

      banner.name <- toupper(banner.name)
      
      if (show.ALL) {
        banner.name <- c(banner.name, c("FINAL_SUBTYPE", "CNV.WGS.CNVKIT.RHO", "RNA.EE"))
        df <- df %>% dplyr::select(all_of(unique(c("TARGET_NAME", banner.name))))
      } else if (show.MPN) {
        # banner.name <- c(banner.name, c("Complex.Karyotype"))
        df <- df %>% dplyr::select(all_of(unique(c("TARGET_NAME", banner.name))))
      } else {
        df <- df %>% dplyr::select(all_of(c("TARGET_NAME", banner.name)))
      }
      
      banner.name <- unique(banner.name)
      
      #----------------------------------------
      # == update lookup.table cols (df) ----
      #----------------------------------------
      
      rownames(df) <- df$TARGET_NAME
      rownames(df) <- NULL
      df$TARGET_NAME <- NULL
      
      # keep individual id rules first, before any title-casing
      
      colnames(df)[colnames(df) == "INDIVIDUAL.ID"] <- "Patient.ID"
      banner.name <- gsub("INDIVIDUAL.ID","Patient.ID",banner.name)
      names(list.my.cols)[names(list.my.cols) == "INDIVIDUAL.ID"] <- "Patient.ID"
      if (!show.individuals) {
        df$INDIVIDUAL.ID <- NULL
      }
      
      ###################################################################
      ## Title-case labels **everywhere** except these cols  ----
      #  (banner.name, df colnames, list.my.cols names)
      ###################################################################
      
      EXCEPT <- c("Patient.ID", "FINAL_SUBTYPE", "CNV.WGS.CNVKIT.RHO", "RNA.EE", "Complex.Karyotype","INDIVIDUAL.ID")
      title_except <- function(x, except = EXCEPT) ifelse(x %in% except, x, stringr::str_to_title(x))
      
      banner.name          <- title_except(banner.name)
      names(list.my.cols)  <- title_except(names(list.my.cols))
      colnames(df)         <- title_except(colnames(df))
      
      #####################################################################
      ## Mode-specific data prep that affects colors/levels, not layout
      #####################################################################
      # browser()
      
      # MPN: factor levels for Complex.Karyotype (color is added to heatmap_colors but check)
      if ("complex.karyotype" %in% tolower(colnames(df))){
          df$Complex.karyotype <- factor(df$Complex.karyotype,
                                         levels = c("complex", "not complex", "not available"))
      }
      
      # ALL: continuous col_fun for RHO and EE, and optional relabeling for legend labels only
      if (show.ALL) {
        col_fun <- colorRamp2(
          c(0, 60, 100),
          c(DM[9], MN[4], "#cc6c4a")
        )
        # assign continuous color functions to these fields in col list
        list.my.cols$`CNV.WGS.CNVKIT.RHO` <- col_fun
        list.my.cols$`RNA.EE`             <- col_fun
      }
      
      # # Make a label vector if user wants to override the printed labels <<<< you can customize banner name labels here!
      # if (is.null(banner.name)) {
      #   # default: only rewrite the label text for display in ALL mode, not the column names
      #   if (show.ALL) {
      #     banner.labels <- banner.name
      #     banner.labels <- gsub("CNV.WGS.CNVKIT.RHO", "WGS.RHO", banner.labels, ignore.case = TRUE)
      #   } else {
      #     banner.labels <- NULL
      #   }
      # }
      
      #####################################################################
      ## Unified controls (backward compatible)
      #####################################################################
      # fallbacks if user did not pass explicit overrides
      
      # browser()
      
      if (is.null(na_col)) {
        na_col <- if (show.ALL) "white" else "darkgrey" 
      }
      
      if (is.null(legend.height)) {
        legend.height <- if (show.ALL || show.MPN) 5 else 20
      }
      
      # if (is.null(legend.direction)) {
      #   legend.direction <- if (show.ALL || show.MPN) "horizontal" else NULL  # NULL means default
      # }
      
      show_names <- if (show.MPN) rep(FALSE, length(banner.name)) else rep(TRUE, length(banner.name))

      #####################################################################
      ## Single HeatmapAnnotation call
      #####################################################################

      show.banner.legends <- as.logical(show.annot.legend)
      
      # browser()
      
      if ((show.individuals) & (!show.individuals.legend)){
        ix = which(banner.name=="Patient.ID")
        show.banner.legends[ix] = FALSE
      }
      
      # browser()
      
      h2 <- HeatmapAnnotation(
        df = df %>% dplyr::select(banner.name),
        name = "my.BottomAnnot",
        col  = list.my.cols,
        na_col = na_col,
        simple_anno_size = unit(ribbon.size, "cm"),
        annotation_height = rep(unit(20, "mm"), length(banner.name)),
        gap = unit(rep(5, ncol(df)), "mm"),
        show_annotation_name = show_names,
        annotation_label = banner.name,              # NULL if not provided
        show_legend = show.banner.legends,
        annotation_name_offset = unit(20, "mm"),
        gp = gpar(col = "black"),
        annotation_name_gp = gpar(fontsize = legend.title.font, fontface = "bold", col = banner.label.col),
        annotation_legend_param = list(
            title_gp = gpar(fontsize = legend.title.font, fontface = "bold"),
            title_position = annot.title.side,
            labels_gp = gpar(fontsize = legend.label.font),
            grid_height = unit(1, "cm"),
            grid_width  = unit(1, "cm"),
            nrow = num.rows.annot.lgd,
            legend_height = unit(legend.height, "cm")
          )
      )
      
      return(list(
        h2 = h2,
        df.updated = df,
        list.my.cols.updated = list.my.cols
      ))
    }
    