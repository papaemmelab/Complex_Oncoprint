prepare_BOTTOM_annotation <- function(df, ## lookup table
                                      
                                      list.my.cols,legend.title.font,legend.label.font, annot.title.side, 
                                      num.rows.annot.lgd, show.annot.legend, 
                                      
                                      ribbon.size= 1, ## the thickness of bottom annotation ribbons in cm 
                                      
                                      show.ALL= FALSE,
                                      show.MPN= FALSE,

                                      banner.name= NULL, 
                                      
                                      show.individuals= FALSE,
                                      show.individuals.legend = FALSE,
                                      
                                      legend.height = NULL,
                                      
                                      # legend.direction = NULL,  ## currently sticking to default

                                      na_col = na_col,
                                      
                                      banner.case.exceptions = NULL, # default kept: c("Patient.ID", "MRD_subtype", "CNV.WGS.CNVKIT.RHO", "RNA.EE", "Complex.Karyotype"),
                                      
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
        df <- df %>% dplyr::select(all_of(unique(c("TARGET_NAME", banner.name))))
      } else if (show.MPN) {
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
      
      EXCEPT <- unique(c(banner.case.exceptions, c("Patient.ID", "CNV.WGS.CNVKIT.RHO", "CNV.WGS.ACE.RHO", "RNA.EE", "Complex.Karyotype")))
      
      # banner.name is already forced to toupper() above, so exceptions are matched case-insensitively;
      # the matched EXCEPT entry's own casing (not x's) is used as the displayed/kept name
      title_except <- function(x, except = EXCEPT) {
        idx <- match(toupper(x), toupper(except))
        ifelse(!is.na(idx), except[idx], stringr::str_to_title(x))
      }
      
      banner.name          <- title_except(banner.name)
      names(list.my.cols)  <- title_except(names(list.my.cols))
      colnames(df)         <- title_except(colnames(df))
      
      #####################################################################
      ## Mode-specific data prep that affects colors/levels, not layout
      #####################################################################

      # MPN: factor levels for Complex.Karyotype (color is added to heatmap_colors but check)
      if ("Complex.Karyotype" %in% colnames(df)){
          df$Complex.Karyotype <- factor(df$Complex.Karyotype,
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
        list.my.cols$`CNV.WGS.ACE.RHO`    <- col_fun
        list.my.cols$`RNA.EE`             <- col_fun
      }
      
      #####################################################################
      ## Unified controls (backward compatible)
      #####################################################################
      # fallbacks if user did not pass explicit overrides
      
      if (is.null(legend.height)) {
        legend.height <- if (show.ALL || show.MPN) 2 else 1
      }
      
      # if (is.null(legend.direction)) {
      #   legend.direction <- if (show.ALL || show.MPN) "horizontal" else NULL  # NULL means default
      # }
      
      show_names <- if (show.MPN) rep(FALSE, length(banner.name)) else rep(TRUE, length(banner.name))

      #####################################################################
      ## Single HeatmapAnnotation call
      #####################################################################

      show.banner.legends <- as.logical(show.annot.legend)
      
      if ((show.individuals) & (!show.individuals.legend)){
        ix = which(banner.name=="Patient.ID")
        show.banner.legends[ix] = FALSE
      }
      
      if (all("CNV.WGS.CNVKIT.RHO" %in% banner.name, "CNV.WGS.ACE.RHO" %in% banner.name)){  ## do not show 2 rho legends
        ix = which(banner.name=="CNV.WGS.CNVKIT.RHO")
        show.banner.legends[ix] = FALSE
      }

      if ("Final_subtype" %in% banner.name & ("WGS_subtype" %in% banner.name | "RNA_subtype" %in% banner.name)){  ## do not show 2 rho legends
        ix = which(banner.name %in% c("RNA_subtype", "WGS_subtype"))
        show.banner.legends[ix] = FALSE
      }
      
      h2 <- HeatmapAnnotation(
        df = df %>% dplyr::select(banner.name),
        name = "my.BottomAnnot",
        col  = list.my.cols,
        na_col = na_col,
        simple_anno_size = unit(ribbon.size, "cm"),
        gap = unit(1, "mm"),
        show_annotation_name = show_names,
        annotation_label = banner.name,              # NULL if not provided
        show_legend = show.banner.legends,
        annotation_name_offset = unit(20, "mm"),
        annotation_name_side = "right",
        annotation_name_rot = 0,
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
    