generate_complex_oncoprint <-  function(muts= muts, cnvs= NULL, svs= NULL ,  # ******* define variants DFs [mut is required]
                                           
                                           cnvs.order= NULL, svs.order= NULL, muts.order= NULL, patients.order= NULL,   # ******* allows pre-defined orders
                                           
                                           surval.data= NULL, show.survival= FALSE, # currently under development
                                           
                                           show.response= FALSE, response.order= NULL, # ******* allows pre-defined orders
                                        
                                           show.another.banner=FALSE, banner.name= NULL, 
                                           
                                           show.ALL= FALSE, ## added specifically for ALL prj. keep it as temp for other adaptations
                                        
                                           show.individuals= FALSE, show.individuals.legend= FALSE,  
                                           
                                           lookup.table= NULL, # ******* pass lookup.table 
                                           
                                           show.sample.names = TRUE, show.border= FALSE, show.multis= FALSE, rem.empty= TRUE, # ******* what params to show in legend?
                                           
                                           split.cols.by = NULL, 
                                           
                                           heatmap.legend.side= "right",
                                           
                                           mut.legend.title.side= "topleft", # ******* HINT: this can only be topleft/topcenter/ etc. otherwise error
                                           
                                           num.rows.heatmap.lgd= NULL, # ******* HEATMAP.legend 
                                           
                                           annot.legend.side= "bottom", 
                                           
                                           annot.title.side= "topleft", # this can only be topleft/topcenter/ etc. otherwise error
                                           
                                           num.rows.annot.lgd= NULL,  # ******* ANNOT.legend 
                                           
                                           min.freq= 1,
                                        
                                           prior.min.Freq= NULL, # For title-ONLY; used when you filter the data in advance but still want to add what you had fileted based on in your title 
                                        
                                           show.title= TRUE, 
                                           
                                           title.str= NULL, 
                                           
                                           save.path= NULL, # ******* title and save path
                                           
                                           save.name= NULL,
                                           
                                           cols.font= 18, rows.font= 18, pct.font= 16, 
                                           legend.label.font= 10, legend.title.font= 14, 
                                           fig.title.font= 18,  barplot.font= 10,  
                                           
                                           multis.dot.size = 0.8, #****FONTs: row.groupname.font is the same as rows.font
                                           
                                           right.w= 13, top.w= 8 , ribbon.size= 2, w=3200, h=1800,  #**** Sizes of barplots and fig 
                                           
                                           axis.side= "left"){
  
  ## must be main branch
  
  
  library(randomcoloR)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(stringr)
  library(plyr)
  library(dplyr)
  
  graphics.off()
  
  #==================================================================================
  #   Written by Noushin Farnoud, Jul 2018. Last Update Dec June 2020  ====
  #----------------------------------------------------------------------------------
  #   The main function to plot the histogram of MUTATIONs [required], CNVs [optional] and Structural Varianrs (SVs) [optional]
  #   mut.order/mut.order/cytogenetics.order
  # === Input Variant info ====
  #==========================================   
  #   muts [required]                     List of mutations (required columns : TARGET_NAME, EFFECT (e.g., missense), GENE (the order is not important, you can have additional cols (e.g., stopgain)))
  #   cnvs [optional]                     List of CNVs (required columns : TARGET_NAME, EFFECT (AMP/DEL/..), VAR_ID (e.g., del(2q)) (the order is not important, you can have additional cols)) [default = NULL] 
  #   svs [optional]                      List of SVs (required columns : TARGET_NAME, EFFECT (e.g., fusion), VAR_ID (the row names you like to use for the SV) (the order is not important, you can have additional cols)) [default = NULL] 
  # === User-defined patient/gene order ====
  #==========================================
  #   muts/cnvs/svs.order [optional]      default row-order of genes/events is based on clustering the data, otherwise specify your desired order for each data type [default = NULL] 
  #   patients.order [optional]           default col-order of patients is based on clustering the data, otherwise specify your order [default = NULL] 
  # === Add annotation rubbons for response/etc ====
  #=================================================
  #   show.sample.source [optional]       Add annotation bar to highlight source of the sample.
  #   show.response [optional]            Add annotation bar for response.
  #   show.individuals [optional]         Add annotation bar to highlight samples that belong to the same patient (useful for dataset with timeline data for patients).
  #   show.individuals.legend             Do you want to add a legend for patients? (only used when show.individuals is set to TRUE) [default= FALSE]
  #   lookup.table                        If any annotation bar is set to on, you must pass a table that sumamrizes sample-feature properties (e.g., TARGET-NAME/RESPONSE)
  # === Control display features ====
  #==========================================
  #   show.sample.names                   Add sample names as the column names [default= TRUE]
  #   show.border [optional]              Add a box around the frequency barplots [default= FALSE]
  #   show.multis [default = FALSE]       If set on, a dot will be displayed on grid elements (gene-sample pair) that have >1 variant. ***NOTE: this currently affects the clustering.
  #   rem.empty                           Remove samples (columns) that have no variant from the oncoprint.
  # === Main Heatmap legend params ====
  #====================================
  #   heatmap.legend.side                 The side that the main mutation-legend is displayed [default= right]
  #   mut.legend.title.side               The position of the mutation-legend title [default= topleft]
  #   num.rows.heatmap.lgd                Number of rows for the mutation-legend 
  # === Annotation ribbon(s) legend params ====
  #=============================================
  #   annot.legend.side                   The side that the legend for the optional added annotation bar(s) (for response, disease, or cell.type) are displayed [default= bottom]
  #   annot.title.side                    The side that the annotation bar legend titles are displayed [default= leftcenter]
  #   num.rows.annot.lgd                  Number of rows for annotation bar legend(s)
  # === Control oncoprint title and display ====
  #==============================================
  #   min.freq                            Only applicable for MUTATIONs data: only show GENEs that have >= min.freq mutations [default = 1] 
  #   show.title [default= TRUE]          Display the figure title. By default this option is set on and if no added title string is (next option) is defined the figure will have a title that reports the total # variants and samples
  #   title.str                           The optional title of the figure, By default this will be followed by the total number of variants in the dataset and number of samples/patients
  #   save.path                           The directory of the output oncoplot : by default the name of the plot is hardcoded as : save.path/"Heatmap_minFreq_",min.freq,".jpg"
  # === Control Font size ====
  #==============================================
  #   cols.font                           Sample names font size (i.e., columns) [default = 18]
  #   rows.font                           Gene/CNV/SV names font size (i.e., rows) [default = 18]
  #   pct.font                            Font size for the percentage frequency that is shown on the left [default = 16]
  #   legend.label.font                   Heatmap/annotation legend font size [default = 10]
  #   legend.title.font                   Font size for the legend title [default = 14]
  #   fig.title.font                      Oncorpting title font size
  #   barplot.font                        Font size for the axis of the frequency barplots that is shown at the top and right of the plot [default = 10]
  # === Control Figure size ====
  #==============================================
  #   right.w                             Size of the area for the right barplot (to display the gene frequency bar) [default = 13]
  #   top.w                               Size of the area for the top barplot (to display the patients frequency bar) [default = 8]
  #   w/h                                 The width and height of the saved figure [default = 3200/1800]
  # 
  #   Contact Noushin Farnoud (rahnaman@mskcc.org) if you faced any error.
  # 
  #    See also example_Heatmap, test_required_fields.
  #==================================================================================
  
  suppressMessages(library("argparse", quietly = TRUE))
  
  if(!is.data.frame(muts)) {muts= as.data.frame(muts)}
  if (!is.null(cnvs) & !is.data.frame(cnvs)) {cnvs = as.data.frame(cnvs)}
  if (!is.null(svs) & !is.data.frame(svs)) {svs = as.data.frame(svs)}
  
  if (is.null(save.path)){
    save.path <- getwd()
    cat(paste0("\n ***** NOTE: You did not pass 'save.path' param when calling the function. The default path used to save the generated oncoprints is --> ", save.path,"\n\n"))
  }
  
  dir.create(file.path(save.path,"TEMP"), showWarnings=FALSE)
  
  my.params = as.list(match.call(expand.dots=FALSE))
  
  ###############################################################
  # == Test Required cols and contents  ====
  ##############################################################
  
  source(file.path("./sub_function/test_required_fields.R"))
  
  rename_IDs <- test_required_fields(muts= muts,  svs=svs, cnvs=cnvs, show.another.banner= show.another.banner, banner.name= banner.name, show.response= show.response, 
                                     split.cols.by= split.cols.by, show.individuals= show.individuals, lookup.table= lookup.table, annot.title.side= annot.title.side)
  
  muts <- rename_IDs$muts
  cnvs <- rename_IDs$cnvs
  svs <- rename_IDs$svs
  lookup.table <- rename_IDs$lookup.table
  REQ.cols <- rename_IDs$required.cols.lookup
  
  ############################################################
  # == Find a subset of Mutations that have >= min.freq variants
  ############################################################
  
  muts <- muts %>% group_by(GENE) %>% mutate(gene.freq= n()) %>% filter(gene.freq>= min.freq) 
  
  muts <- as.data.frame(muts)
  
  if (nrow(muts)==0){
    cat(paste("\n You have 0 Mutations to show for this data set with min Freq =",min.freq,"\nreturning NULL"))
    return(list(ht.obj = NULL, annotation_legend_list= NULL, heatmap_legend_list= NULL,
                onco.samples= NULL))
  }
  
  ############################################################
  # Filter svs and cnvs based on min.freq
  ############################################################
  # SVs ---
  #=====================
  if (!(is.null(svs))){
    
    svs.test <- svs %>% group_by(GENE) %>% mutate(gene.freq= n()) %>% filter(gene.freq >= min.freq) 
  
    if (nrow(svs.test)==0){
      cat(paste0("\n >>>>>>>  There are no SVs with min.freq you specified; So, plotting svs with at least 1 hit instead..."))
      svs <- svs %>% group_by(GENE) %>% mutate(gene.freq= n()) %>% filter(gene.freq>= 1) 
    } else {svs <- svs.test}
  
  rm(svs.test)
  
  }
  
  # CNVs ---
  #=====================
  if (!(is.null(cnvs))){
    
    cnvs.test <- cnvs %>% group_by(GENE) %>% mutate(gene.freq= n()) %>% filter(gene.freq >= min.freq) 
  
    if (nrow(cnvs.test)==0){
      cat(paste0("\n >>>>>>> There are no CNVs with min.freq you specified; So, plotting cnvs with at least 1 hit instead..."))
      cnvs <- cnvs %>% group_by(GENE) %>% mutate(gene.freq= n()) %>% filter(gene.freq >= 1) 
    } else {cnvs <- cnvs.test}
    
    rm(cnvs.test)
  }
  
  ##########################################
  # Prepare data for complex heatmap  ====
  ##########################################
  
  source(file.path("./sub_function/initialize_data.R"))
  Init.List <- initialize_data(data, muts, cnvs, svs, muts.order, cnvs.order, svs.order, lookup.table, REQ.cols, save.path, my.params)
  
  saveFile.1 <-  Init.List$saveFile.1
  saveFile.2 <-  Init.List$saveFile.2
  my.fonts <-  Init.List$font.obj
  
  data <- Init.List$data
  muts <-  Init.List$muts
  
  gene.list= Init.List$gene.list
  
  ####################################
  # Load colors  ====
  ####################################
  
  cat(paste0("\nLoading default oncopring colors...\n"))
  
  source(file.path("./sub_function/heatmap_colors.R"))
  list.ht.colors <- heatmap_colors()
  
  ###############################################################
  # == Adjust EFFECT to uniform texts  ====
  ##############################################################
  valid.effects <- tolower(rename_IDs$valid.effects)
  
  source(file.path("./sub_function/make_uniform_EFFECT_values.R"))
  data <- make_uniform_EFFECT_values(data, valid.effects)
  
  ###############################################################
  # == Prepare the Heatmap rows and columns  ====
  ##############################################################
  
  # data <- unique(data) # this would eliminate multis for gene mutations of the same type in a sample
  
  # == GENES is the FINAL order of genes that will be the row names
  GENES <- data.frame(genes= gene.list$GENE, EFFECT=gene.list$LAB)
  
  # BE CAREFUL do not unique data : you will loose cases where a gene has multiple variants in the same patient 
  
  SAMPLES = as.data.frame(with(data, table(TARGET_NAME)),stringsAsFactors = FALSE)
  
  ###############################################################
  # == Prepare M and populate matrix of variants  ====
  ##############################################################
  
  source(file.path("./sub_function/prepare_fill_M.R"))
  M.List <- prepare_fill_M(data, SAMPLES, GENES, lookup.table, rem.empty, gene.list)
  
  M <- M.List$M
  gene.order <- M.List$gene.order
  events <- M.List$events
  
  ###############################################################
  # == Define "alter_fun" =====
  ##############################################################
  
  cat(paste0("\nLoading Default ALTER func...\n"))
  
  source(file.path("./sub_function/define_ALTER_fun.R"))
  alter_fun <- define_ALTER_fun(list.ht.colors, multis.dot.size)
  
  ###############################################################
  # == Define Labels for MUT/CNV/... segments  =====
  ##############################################################
  
  EFFECT.all <- list(variants = c("missense","stop_gain","frameshift_indel",
                                  "inframe_indel","splice_site_variant",
                                  "initiator_codon_change",
                                  "complex",
                                  "unknown",  "complex_karyotype",
                                  "amp",     
                                  "del",  "loh",  "inv",      
                                  "fusion",   
                                  "trans", "other_svs","tdup","dup","rearr",
                                  "add","der",
                                  "other_snvs",
                                  "other_cnvs",
                                  "multi_hit","unavailable","normal","karyotypic_abnormal"), 
                     
                     labels= c("Missense","Stop-gain","Frameshift indel",
                               "Inframe indel","Splicing variant",
                               "Initiator_codon change",
                               "Complex", 
                               "Unknown", "Complex karyotype",  
                               "Amplification",
                               "Deletion","cnLOH", "Inversion",
                               "Fusion",
                               "Translocation","Other SVs","Tandem duplication", "Duplication","Rearrangement",
                               "Add.","Der.",
                               "Other mutations",
                               "Other CN alterations",
                               "Multiple variants","Unavailable","Normal","Karyotypic abnormal"))
  
  EFFECT <- list (variants = EFFECT.all[[1]][EFFECT.all[[1]] %in% data$EFFECT],
                  labels = EFFECT.all[[2]][EFFECT.all[[1]] %in% data$EFFECT]
  )
  
  
  LABS <- factor(gene.list$LAB, levels=c("Substitusions/Indels","Cytogenetics","CNVs", "SVs"))
  
  #################################
  # == Top-annotation (1)  ====
  #################################
  
  cat(paste0("\nPrepare Top Annotation...\n"))
  
  # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual' & brewer.pal.info$colorblind==TRUE,]
  # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  df <-  data.frame(TARGET_NAME= colnames(M))
  df$TARGET_NAME <- as.character(df$TARGET_NAME)
  
  surv.df = df
  
  if (show.another.banner | show.response | show.individuals){
    
    df <- merge(df, lookup.table[,REQ.cols], by=c("TARGET_NAME"), all.x = TRUE)
    
    if (any(is.na(df$INDIVIDUAL.ID))) {
      stop("\n***An error occured in merging dataframe of variants with LOOKUP.TABLE. \nYou have at least one sample where INDIVIDUAL.IDs= NA.\n This can occur if the key TARGET_NAME in MUTs and LOOKUP.TABLE are inconsistent!\n")
    }
  }
  
  #################################################################################
  # == PREPARE HeatmapAnnotation Obj (Response, Sample.Source, Individual.ID) ==== 
  #################################################################################
  
  # df$TARGET_NAME <- NULL
  
  #====================================================
  # == Create COLOR palletes for HeatmapAnnotation ====
  #====================================================
  
  list.my.cols <- list()
  
  show.annot.legend <- c()
  
  #============================
  # if showing RESPONSE ====  
  #============================
  
  if (show.response){
    cat(paste0("\nPrepare RESPONSE...\n"))
    
    resp.col <- list.ht.colors$response.colors[names(list.ht.colors$response.colors) %in% unique(lookup.table$RESPONSE)]
    
    list.my.cols$RESPONSE <- resp.col
    
    show.annot.legend <- c(show.annot.legend, "TRUE")
    
  } 
  
  #============================
  # if showing NEW.BANNER ====   
  #============================
  
  if (show.another.banner){
    
    if (show.ALL) {
      
      highANY2.col <- list.ht.colors$ALL.SUBTYPE[names(list.ht.colors$ALL.SUBTYPE) %in% unique(lookup.table$HIGH_ANY2_TYPE)]
      
      list.my.cols$HIGH_ANY2_TYPE <- highANY2.col
      
      show.annot.legend <- c(show.annot.legend, "TRUE")
      
     
      
      RNA.col <- list.ht.colors$ALL.SUBTYPE[names(list.ht.colors$ALL.SUBTYPE) %in% unique(lookup.table$RNA.SUBTYPE)]
      
      list.my.cols$RNA.SUBTYPE <- RNA.col
      
      show.annot.legend <- c(show.annot.legend, "TRUE")
      
      
      
      DNA.col <- list.ht.colors$ALL.SUBTYPE[names(list.ht.colors$ALL.SUBTYPE) %in% unique(lookup.table$WGS.SV.CNA.SUBTYPE)]
      
      list.my.cols$WGS.SV.CNA.SUBTYPE <- DNA.col
      
      show.annot.legend <- c(show.annot.legend, "TRUE")
      
      
      gender.col <- list.ht.colors$GENDER[names(list.ht.colors$GENDER) %in% unique(lookup.table$GENDER)]
      
      list.my.cols$GENDER <- gender.col
      
      show.annot.legend <- c(show.annot.legend, "TRUE")
      
    } else {
      source(file.path("./sub_function/add_new_banner.R"))
      BannerList <- add_new_banner(banner.name, lookup.table, list.my.cols, show.annot.legend)
      list.my.cols= BannerList$list.my.cols
      new.banner.col = BannerList$new.banner.col
      show.annot.legend= BannerList$show.annot.legend
    }

    # list.my.cols= BannerList$list.my.cols
    # new.banner.col = BannerList$new.banner.col
    # show.annot.legend= BannerList$show.annot.legend
    # 
    # list.my.cols$GROUP["Not Available"] <- "#ffffff" # remove later. meant to improve UK-ALL viz
    # new.banner.col["Not Available"] <- "#ffffff"
    # 
    # list.my.cols$GROUP["Hypodiploid"] <- "#fa9fb5" # remove later. meant to improve UK-ALL viz
    # new.banner.col["Hypodiploid"] <- "#fa9fb5"
    # 
    # list.my.cols$NUM.EVIDENCE["5"] <- "#de2d26" # remove later. meant to improve UK-ALL viz
    # new.banner.col["5"] <- "#de2d26"

  }
  
  ####################################
  # == Prepare the survival data ====
  ####################################
  
  if (!is.null(surval.data) & show.survival){
    
    surv.df <- merge(surv.df, surv.info, by=c("TARGET_NAME"))
    
    surv.df$status.col= ifelse(surv.df$Death.status=="1","red","blue")
    
    surv.df$pch= ifelse(surv.df$Death.status=="1",13,16)
    
    # show.annot.legend <- c(show.annot.legend, "TRUE")
  }
  
  ####################################
  # if showing PATIENTS   ====
  ####################################
  
  if (show.individuals){
    
    n2 <- length(unique(lookup.table$INDIVIDUAL.ID))
    
    # color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    # indiv.col <- color[1:n2]
    
    indiv.col <-  distinctColorPalette(n2)
    
    names(indiv.col) <- unique(lookup.table$INDIVIDUAL.ID)
    list.my.cols$INDIVIDUAL.ID <- indiv.col
    
    show.annot.legend <- c(show.annot.legend, "TRUE")
    
    rm(n2)
  } 
  
  ################################################################
  #### Define Top Annotation ====
  ################################################################
  ###### This enforces to have at least ONE RIBBON + SURV.dots
  ###### update it in future if necessary
  #========================================================
  
  source(file.path("./sub_function/prepare_TOP_annotation.R"))
  h1 <- prepare_TOP_annotation(list.colors,show.border,axis.side,barplot.font, legend.title.font, top.w)
  
  # ###############################################################
  # # == Create a legend for survival if show.survival= TRUE ====
  # ##############################################################
  # currently off as it can not merge properly with response
  
  if (show.survival){
    lgd_list = list(
      Legend(labels = c("Dead", "Alive"),
             labels_gp = gpar(fontsize = legend.label.font),
             title_gp = gpar(fontsize = legend.title.font, fontface="bold"),
             title = "Survival Status",
             title_position = "topcenter",
             type = "points", pch = c(13,16), size = unit(1, "cm"),
             background= "white",
             nrow = 2,
             grid_width = unit(1, "cm"),
             grid_height = unit(1, "cm"),
             legend_gp = gpar(fontsize = 23, col = c("red","blue"), size = unit(8, "cm"))
      )
    )
  } else {
    lgd_list= NULL
  }
  ###############################################################
  # == Create a legend for multis if show.multis= TRUE ==== MUST BE AFTER h1
  ##############################################################
  
  if (show.multis){  
    
    ht.list = list(Legend(labels =  c(">1 variant"),
                          labels_gp = gpar(fontsize = legend.label.font),
                          type = "points",
                          pch = 21,
                          size = unit(0.5, "cm"),
                          legend_gp = gpar(col = "black", fill= "#FAEFD1", lwd= 1, fontsize = legend.label.font),
                          background = NULL, grid_height = unit(1, "cm"),
                          grid_width = unit(1, "cm")))
  } else {
    ht.list = NULL
  }
  
  ###############################################################
  # == Set title/params and figure name ==== 
  ##############################################################
  
  cat(paste0("\nSet row/col orders...\n"))
  
  if (is.null(prior.min.Freq)){
    cat(paste("\n **** You chose there was NO prior filtering based on Gene Freq, so for the title I am going with min.freq as-called by the func= ", min.freq))
    Gene.Freq = min.freq
  } else {
    cat(paste("\n **** You chose there was a prior filtering (based on Gene Freq), and I am trusting you with your enetered param and directly use it in the title= ", prior.min.Freq))
    Gene.Freq = prior.min.Freq
  }
  
  if (show.title){
    my.title <- paste0(title.str," \n# Alterations= ", nrow(data),"; # Genes with >= ", Gene.Freq ," mutations = ",length(unique(muts$GENE)),"; # Samples =", ncol(M))
    # my.title <- paste0(title.str," \n# Alterations= ", data %>% nrow(), "; # Distinct SVs with >= ", Gene.Freq ," alterations = ", svs %>% filter(GENE!="NONE") %>% select(GENE) %>% unique() %>% nrow(),
    #                    "; # Samples =", length(unique(data$TARGET_NAME)))
    
  } else {
    my.title <- NULL    
  }
  
  if (is.null(patients.order)){
    column_order = NULL
  } else {
    column_order= as.character(patients.order)
  }
  
  if (!is.null(muts.order) & !is.null(cnvs) & is.null(cnvs.order)){ # is we have cnvs and also decided to pass cnvs.order
    cnvs.order <- as.character(unique(cnvs$GENE))
  }
  
  if (!is.null(muts.order) & !is.null(svs) & is.null(svs.order)){ # is we have cnvs and also decided to pass cnvs.order
    svs.order <- as.character(unique(cnvs$GENE))
  }
  
  if (is.null(muts.order) & is.null(cnvs.order) & is.null(svs.order)){
    row_order = NULL
    
  } else {
    
    if (is.null(cnvs)) {cnvs.order= NULL} # in case you choose to set CNVs= NULL but forget to set the cnvs.order= NULL
    if (is.null(svs)) {svs.order= NULL}
    
    row_order= c(muts.order,cnvs.order,svs.order)
  }
  
  ###############################################################
  # == Generate Simple.ONCOPRINT ==== 
  ##############################################################
  cat(paste0("\nGenerating simple oncoprint...\n"))
  
  num.my.lgd.rows <- num.rows.heatmap.lgd 
  
  
  
  source(file.path("./sub_function/draw_basic_oncoprint.R"))
  
  simple.ht <- draw_basic_oncoprint(M, EFFECT, alter_fun, 
                                    
                                    saveFile= saveFile.1,
                                    
                                    list.colors= list.ht.colors, 
                                    
                                    top_annotation= h1, 
                                    
                                    heatmap.legend.side= heatmap.legend.side,
                                    annot.legend.side= annot.legend.side,
                                    
                                    heatmap.legend.list= ht.list,
                                    annot.legend.list= lgd_list,
                                    
                                    column_order= column_order,
                                    right.w= 13, 
                                    LABS= LABS, 
                                    
                                    font.obj= my.fonts, 
                                    num.rows.heatmap.lgd= num.my.lgd.rows,
                                    
                                    w=w,
                                    h=h,
                                    
                                    fig.title= NULL,  
                                    show.border= TRUE, show.sample.names= TRUE)
  
  ##======================================================
  ## Finished plotting BASIC oncoprint ====
  ##======================================================
  # If no added RESPONSE/ANNOTATIONBAR/etc was selected, 
  #  but wanted show.sample.names=FASLE repeat the basic 
  #  heatmap plot, but with FALSE option.
  ##======================================================
  ##===============================================================================
  ## *** IMPORTANT: Get the sample.order of simple.ht to sort the annotation, UNLESS
  ##                the user has strict patient order in input
  ##===============================================================================
  cat(paste0("\nFetch the order of samples (cols) from simple.ht...\n"))
  
  if (is.null(patients.order)){
    new.column_order <- colnames(M)[column_order(simple.ht)] #this is the order of the simple oncoprint with basic clustering
    
  } else {
    new.column_order <- patients.order
  }
  
  # my.temp.column_order <- colnames(simple.ht@matrix)
  
  #################################################################################################
  #################################################################################################
  #### Start Complex plot. ====
  #################################################################################################
  #################################################################################################
  
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
  #### BOTTOM ANNOTATION ====
  #############################
  
  if ((show.another.banner) | (show.response) | (show.individuals) ) {
    
    if (!show.individuals){
        df$INDIVIDUAL.ID <- NULL
    }
    
    source(file.path("./sub_function/prepare_BOTTOM_annotation.R"))
    BotAnnot <- prepare_BOTTOM_annotation(df, list.my.cols,legend.title.font,legend.label.font,annot.title.side,num.rows.annot.lgd, show.annot.legend, ribbon.size, show.individuals= show.individuals)
    
    h2 <- BotAnnot$h2
    df <- BotAnnot$df.updated
    list.my.cols <- BotAnnot$list.my.cols.updated
    
  } else  {
    h2 = NULL # for example, you do not have any added bottom annotation but still like to see multis
    
  }
  
  # samples.order.mod <- colnames(simple.ht@matrix)
  
  samples.order.mod <- new.column_order
  
  ###############################################################
  ###############################################################
  # == Generate COMPLEX.ONCOPRINT ==== 
  ##############################################################
  
  cat(paste0("\nGenerate Final COMPLEX oncoprint ...\n"))
  
  
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
  
  
  ht <- oncoPrint(M, get_type = function(x) strsplit(x, ";")[[1]],
                  
                  alter_fun = alter_fun, col = append(list.ht.colors$mut.colors, list.ht.colors$cyto.colors),
                  
                  #axis_gp = gpar(fontsize = 8, fontface="bold"), # obsolete param
                  
                  column_order = samples.order.mod,
                  
                  row_order = row_order, #control the order of genes (rows)
                  
                  remove_empty_columns = rem.empty,
                  
                  show_column_names = show.sample.names,
                  
                  column_split= split.cols.order, # this supposed to add a vertical gap between columns based on a selected characteristic of samples (SPLIT col in lookup-table)
                  
                  column_gap = unit(5, "mm"),
                  
                  # === Gene barplots on the left ====
                  
                  bottom_annotation= h2,
                  top_annotation = h1,
                  
                  right_annotation = rowAnnotation(row_bar = anno_oncoprint_barplot(type= NULL,
                                                                                    border= show.border, 
                                                                                    axis_param = list(side= "top", 
                                                                                                      gp= gpar(fontsize= barplot.font, fontface="bold"))),
                                                   annotation_width= unit(right.w,"cm")),   ## controls the width of the row.barplots
                  
                  #show_row_barplot = TRUE, # obsolete param
                  #row_barplot_width = unit(right.w, "cm"), # obsolete param
                  
                  split= LABS,
                  
                  # ==========================================
                  # ==========================================
                  
                  # === Title ====
                  column_title = my.title,
                  column_title_gp = gpar(fontsize = fig.title.font, fontface = "bold"), # title font-size
                  gap = unit(10, "mm"),
                  
                  # === Column/Sample names ====
                  column_names_gp = gpar(cex=1, col= "black", fontsize = cols.font, fontface="bold"), #default size = 18
                  column_names_max_height= unit(20,"cm") , # adjust this to control the name of samples (col names)
                  
                  # === Percent ====
                  pct_gp=gpar(fontsize = pct.font, fontface = "bold", col="black"), # specific control over percentage info on the left (add col="blue" to change colors)
                  row_names_gp = gpar(fontsize = rows.font, fontface="bold"), # gene-names and percent (if not prc_gp is defined above)
                  row_title_gp = gpar(fontsize =rows.font+3, col="blue",fontface = "bold"),
                  
                  # === Legend ====
                  # heatmap_legend_param = gg_list # list of list does not work here!
                  
                  heatmap_legend_param = list(title = "Alterations", at = EFFECT$variants,
                                              labels = EFFECT$labels,
                                              heatmap_legend_list= ht.list,
                                              title_gp = gpar(fontsize = legend.title.font, fontface="bold"),
                                              title_position = mut.legend.title.side,
                                              # title_position= "topleft",
                                              labels_gp = gpar(fontsize = legend.label.font),
                                              grid_height= unit(1, "cm"), # size of the mutation legend color-boxes
                                              nrow=num.rows.heatmap.lgd,
                                              grid_width= unit(1, "cm"),
                                              legend_height = unit(20, "cm"))
  ) 
  
  ##======================================================
  ## Draw simple.ht ====
  ##======================================================
  
    jpeg(saveFile.2, width=w, height=h, pointsize =14, res = 100)
    
  if (heatmap.legend.side== annot.legend.side){
    draw(ht, split= LABS,  merge_legend = TRUE,  heatmap_legend_side = heatmap.legend.side, annotation_legend_side = heatmap.legend.side, annotation_legend_list = lgd_list,
         heatmap_legend_list = ht.list)
  } else {
    draw(ht, split= LABS,  merge_legend = FALSE,  heatmap_legend_side = heatmap.legend.side, annotation_legend_side = annot.legend.side, annotation_legend_list = lgd_list,
         heatmap_legend_list = ht.list)
  }
  
  dev.off()
  
  if (!is.null(split.cols.by)){
    cat(paste("\n*** NOTE ***You can not get the final ordered list of samples (column_order) if you have chosen to split the columns by RESPONSE.\n 
              You can still get the list if you re-run the function and set split.by.response= FASLE. \n---> Future dev."))
    final.sample_order = NULL
  } else {
    final.sample_order <- colnames(M)[column_order(ht)]
  }
  
  
  
  ###############################################################################
  # === if automatic clustering is done, you can use the codes below 
  # ==== to decipher the exact order of clustered samples (add these to the calling code)
  # ==============================================================================
  # col.list <- column_order(ht)
  # htnames <- names(column_order(ht))
  # col.orders <- col.list[[htnames[2]]]
  # sample_order <- colnames(M)[col.orders]
  # =============================================================================
  
  # == ideas for future dev. 
  # draw(ht, padding = unit(c(40, 40), "mm")) 
  
  # == ideas for future dev. 
  # decorate_annotation("RESPONSE", {grid.text("value", unit(-2, "mm"), just = "right")})
  ###############################################################################
  
  cat(paste("\n\nThe file is saved at",saveFile.2,"\n"))
  
  return(list(ht.obj = ht, annotation_legend_list= lgd_list, heatmap_legend_list= ht.list,
              onco.samples= final.sample_order))
  
}