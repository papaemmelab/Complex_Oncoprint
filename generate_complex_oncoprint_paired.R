generate_complex_oncoprint_paired <-  function(muts= muts, cnvs= NULL, svs= NULL ,  # define variants DFs
                                               
                                cnvs.order= NULL, svs.order= NULL, muts.order= NULL, patients.order= NULL,   # define pre-set orders
                                
                                surval.data= NULL,
                                
                                show.response= FALSE, response.order= NULL, show.another.banner=FALSE, banner.name= NULL, show.individuals= FALSE, show.individuals.legend= FALSE, show.survival= FALSE, # ******* user-defined ORDERs
                                
                                lookup.table= NULL, #******* pass lookup.table 
                                
                                show.sample.names = TRUE, show.border= FALSE, show.multis= FALSE, rem.empty= TRUE, # ******* what params to show in legend?
                                
                                heatmap.legend.side= "right", 
                                
                                mut.legend.title.side= "topleft", # this can only be topleft/topcenter/ etc. otherwise error
                                
                                num.rows.heatmap.lgd= NULL, #******* HEATMAP.legend 
                                
                                annot.legend.side= "bottom", 
                                
                                annot.title.side= "topleft", # this can only be topleft/topcenter/ etc. otherwise error
                                
                                num.rows.annot.lgd= NULL,  #******* ANNOT.legend 
                                
                                min.freq= 1, show.title= TRUE, title.str= NULL, save.path= save.path, #******* title and save path
                                
                                save.name= NULL,
                                
                                cols.font= 18, rows.font= 18, pct.font= 16, 
                                legend.label.font= 10, legend.title.font= 14, 
                                fig.title.font= 18,  barplot.font= 10,  
                                
                                multis.dot.size = 0.8, #****FONTs: row.groupname.font is the same as rows.font
                                 
                                right.w= 13, top.w= 8 , robbon.size= 2, w=3200, h=1800,  #**** Sizes of barplots and fig 
                                
                                axis.side= "left"){

  
  library(randomcoloR)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(plyr)
  library(dplyr)
  
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
  
  ###############################################################
  # == Test Required cols and contents  ====
  ##############################################################
  
  source(file.path("./sub_function/test_required_fields.R"))

  rename_IDs <- test_required_fields(muts= muts,  svs=svs, cnvs=cnvs, show.another.banner= show.another.banner, banner.name= banner.name, show.response= show.response, show.individuals= show.individuals, lookup.table= lookup.table)
  
  muts <- rename_IDs$muts
  cnvs <- rename_IDs$cnvs
  svs <- rename_IDs$svs
  lookup.table <- rename_IDs$lookup.table
  REQ.cols <- rename_IDs$required.cols.lookup
  
  ###############################################
  # Create FONT.obj for calling simple.ht  ====
  ###############################################
  font.obj <- list(fig.title.font= fig.title.font,
                   legend.title.font= legend.title.font,
                   cols.font= cols.font,
                   pct.font= pct.font,
                   rows.font= rows.font,
                   barplot.font= barplot.font,
                   legend.label.font= legend.label.font
                   )

  my.fonts= font.obj
  
  ############################################################
  # == Find a subset of Mutations that have >= min.freq variants
  ############################################################
  
  muts <- muts %>% group_by(GENE) %>% mutate(gene.freq= n())
  
  muts <- muts %>% filter(gene.freq>= min.freq) 
  
  muts <- as.data.frame(muts)
  
  ###############################################################
  # == Add MUTATIONS  ====
  ##############################################################
  
  source(file.path("./sub_function/sort_variants.R"))
  
  A <- sort_variants(muts, muts.order, group.label= "Substitusions/Indels", variants.class= "Mutations")
  
  data <- A[[1]][,c("TARGET_NAME","GENE","EFFECT")]
  
  gene.list <- A$gene.list
  
  ###############################################################
  # == Add Structural Variants  ====
  ##############################################################
  
  if (!is.null(svs)){
    
    source(file.path("./sub_function/sort_variants.R"))
    
    B <- sort_variants(svs, svs.order, group.label= "SVs", variants.class= "SVs")
    
    data <- rbind(data, B[[1]][,c("TARGET_NAME","GENE","EFFECT")])
    
    gene.list <- rbind(gene.list,B$gene.list)
    
  }
  
  ###############################################################
  # == Add CNVs  ====
  ##############################################################
  
  if (!is.null(cnvs)){
    
    source(file.path("./sub_function/sort_variants.R"))
    
    D <- sort_variants(cnvs, cnvs.order, group.label= "CNVs", variants.class= "CNVs")
    
    data <- rbind(data, D[[1]][,c("TARGET_NAME","GENE","EFFECT")])
    
    gene.list <- rbind(gene.list,D$gene.list)
  }
  
  ###############################################################
  # == Prepare the Heatmap rows and columns  ====
  ##############################################################
  
  # data <- unique(data) # this would eliminate multis for gene mutations of the same type in a sample
  
  # == GENES is the FINAL order of genes that will be the row names
  GENES <- data.frame(genes= gene.list$GENE, EFFECT=gene.list$LAB)
  
  # BE CAREFUL do not unique data : you will loose cases where a gene has multiple variants in the same patient 
  
  SAMPLES = as.data.frame(with(data, table(TARGET_NAME)),stringsAsFactors = FALSE)
  
  ###############################################################
  # == Adjust EFFECT to uniform texts  ====
  ##############################################################
  data$EFFECT <- tolower(data$EFFECT)
  
  data$EFFECT <- gsub("^missense$|^non_synonymous_codon$|^missense_codon$","missense", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^stop_gained$|^stop_gain$|^stop_lost$|^stop_retained_variant$","stop_gain", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^splice_site_variant$","splice_site_variant", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^initiator_codon_change$","initiator_codon_change", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^inframe_codon_loss$|^inframe_indel$|^inframe_deletion$|^inframe_codon_gain$|^inframe_insersion$|^inframe_variant$","inframe_indel", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^complex_change_in_transcript$|^complex$","complex",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^other_snvs$","other_snvs",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^frameshift_indel$|^frameshift_del$|^frameshift_variant$","frameshift_indel",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^amp$|^amplification$|^gain$|^CN-gain$","amp", ignore.case = TRUE,  data$EFFECT)
  data$EFFECT <- gsub("^del$|^deletion$|^loss$|^CN-del$", "del", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^LOH$","loh", ignore.case = TRUE,  data$EFFECT)
  data$EFFECT <- gsub("^unknown$","unknown",ignore.case = TRUE, data$EFFECT)

  data$EFFECT <- gsub("^inv$|^inversion$", "inv", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^tandem duplications$|^tandem_duplications$|^tandem dup$","tdup", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^fusion$|^fus$","fusion",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^translocation$|^trans$","trans",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^other_svs$","other_svs",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^n/e$|^inconclusive$|^n_e$|^n_a$","inconclusive",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^complex_karyotype$","complex_karyotype",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^normal_karyotype$|^normal$","normal",ignore.case = TRUE, data$EFFECT)
  
  ####################################
  # Test EFFECTs  ====
  ####################################
  
  valid.effects <- tolower(c("missense","stop_gain","splice_site_variant","initiator_codon_change",
                     "inframe_indel","complex","other_snvs","frameshift_indel","unknown", "fusion",
                     "AMP","DEL","LOH","INV","TDUP","FUS","TRANS","OTHER_SVs","complex_karyotype", 
                     "DER","ADD","Normal","inconclusive"
                     ))
  
  data$EFFECT <- tolower(data$EFFECT)
  
  invalid.effects <- setdiff(data$EFFECT, valid.effects)
  
  if (length(invalid.effects) >0) {
    stop(cat(paste("\nThese variant(s) EFFECTs are not valid: ", paste(invalid.effects, collapse = ", "))))
  } else {
    cat(paste0("\nAll EFFECTs are valid. Good to go...\n"))
  }
  
  ####################################
  # Load colors  ====
  ####################################
  
  cat(paste0("\nLoading default oncopring colors...\n"))
  
  source(file.path("./sub_function/heatmap_colors.R"))
  list <- heatmap_colors()
  
  ###############################################################
  # == Prepare M matrix of variants  ====
  ##############################################################
  
  if (rem.empty) {
    M.num.cols <- length(unique(SAMPLES$TARGET_NAME))
    M.col.names <- SAMPLES$TARGET_NAME
  } else {
    M.num.cols <- length(unique(lookup.table$TARGET_NAME))
    M.col.names <- unique(lookup.table$TARGET_NAME)
  }
  
  M <- as.data.frame(matrix(0, nrow = length(GENES$genes), ncol = M.num.cols))
  
  row.names(M) <- GENES$genes
  
  colnames(M) <- M.col.names
  
  gene.order <- gene.list$GENE
  
  ###############################################################
  # == Add the Event.Type in the Matrix ====
  ##############################################################
  
  cat(paste0("\nGenerating the matrix of mutations (M)...\n"))
  
  M[M==0] = ""
  
  events <- factor(unique(data$EFFECT), levels=c("unknown","other_snvs","missense","splice_site_variant","initiator_codon_change",
                                                 "complex","complex_karyotype","stop_gain","inframe_indel",
                                                 "frameshift_indel","amp","del","loh","inv","fusion","trans","tdup","add","der",
                                                 "other_svs","inconclusive","normal"))
  events <- events[order(events)]
  
  events <- as.character(events)
  
  for (i in 1: length(events)){
    
    temp <- subset(data, EFFECT==events[i])
    
      for (j in 1:nrow(temp)) {
          
          M[temp$GENE[j], temp$TARGET_NAME[j]] <- paste0(M[temp$GENE[j], temp$TARGET_NAME[j]], unique(temp$EFFECT),";", collapse = "")
      }
  }
  
gogo=1
  ###############################################################
  # == Define "alter_fun" =====
  ##############################################################

  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.52, "mm"), h-unit(0.52, "mm"), gp = gpar(fill = "#f0f0f0", col = NA)) # alpha=0.5
    },
    unknown = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$mut.colors[["unknown"]][1], col = NA))
    },
    other_snvs = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$mut.colors[["other_snvs"]][1], col = NA))
    },
    missense = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$mut.colors[["missense"]][1]  , col = NA))
    },
    splice_site_variant = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$mut.colors[["splice_site_variant"]][1], col = NA))
    },
    initiator_codon_change = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$mut.colors[["initiator_codon_change"]][1], col = NA))
    },
    complex = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$mut.colors[["complex"]][1], col = NA))
    },
    stop_gain = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h-unit(0.3, "mm"), gp = gpar(fill = list$mut.colors[["stop_gain"]][1], col = NA))
    },
    inframe_indel = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h*0.33, gp = gpar(fill = list$mut.colors[["inframe_indel"]][1]  , col = NA))
    },
    frameshift_indel = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h*0.33, gp = gpar(fill = list$mut.colors[["frameshift_indel"]][1],  col = NA))
    },
    complex_karyotype = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$cyto.colors[["complex_karyotype"]][1], col = NA))
    },
    amp = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$cyto.colors[["amp"]][1], col = NA))
    },
    del = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$cyto.colors[["del"]][1], col = NA))
    },
    loh = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$cyto.colors[["loh"]][1], col = NA))
    },
    inv = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$cyto.colors[["inv"]][1], col = NA))
    },
    fusion = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h*0.33, gp = gpar(fill = list$cyto.colors[["fusion"]][1],  col = NA))
    },
    trans =function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h*0.33, gp = gpar(fill = list$cyto.colors[["trans"]][1],  col = NA))
    },
    other_svs =function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h*0.33, gp = gpar(fill = list$cyto.colors[["other_svs"]][1],  col = NA))
    },
    tdup =function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h*0.33, gp = gpar(fill = list$cyto.colors[["tdup"]][1],  col = NA))
    },
    inconclusive=function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h-unit(0.3, "mm"), gp = gpar(fill = list$cyto.colors[["inconclusive"]][1],  col = NA))
    },
    der = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$cyto.colors[["der"]][1],  col = NA))
    },
    add = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$cyto.colors[["add"]][1],  col = NA))
    },
    multi_hit = function(x, y, w, h) {
      grid.points(x, y, pch = 21, size = unit(multis.dot.size, "cm"), gp = gpar(col = "black", fill= "#FAEFD1"))
    },
    normal = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$cyto.colors[["normal"]][1]  , col = NA))
    }
  )
  
  EFFECT <- list(variants = c("missense","stop_gain","frameshift_indel",
                              "inframe_indel","splice_site_variant",
                              "initiator_codon_change",
                              "complex",
                              "unknown",  "complex_karyotype",
                              "amp",     
                              "del",  "loh",  "inv",      
                              "fusion",   
                              "trans", "other_svs","tdup",
                              "add","der",
                              "other_snvs",
                              "multi_hit","inconclusive","normal"), 
                 
                 labels= c("Missense","Stop-gain","Frameshift indel",
                           "Inframe indel","Splicing variant",
                           "Initiator_codon change",
                           "Complex", 
                           "Unknown", "Complex karyotype",  
                           "Amplification",
                           "Deletion","cnLOH", "Inversion",
                           "Fusion",
                           "Translocation","Other SVs","Tandem duplication",
                           "Add.","Der.",
                           "Other mutations",
                           "Multiple variants","Inconclusive","Normal"))
  
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
      stop("\n***An error occured in merging dataframe of variants with LOOKUP.TABLE, and you have INDIVIDUAL.IDs= NA.\n This can occur if the key TARGET_NAME in MUTs and LOOKUP.TABLE are inconsistent!\n")
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
    
    resp.col <- list$response.colors[names(list$response.colors) %in% unique(lookup.table$RESPONSE)]
    
    list.my.cols$RESPONSE <- resp.col
    
    show.annot.legend <- c(show.annot.legend, "TRUE")
    
  } 
  
  #============================
  # if showing NEW.BANNER ====   
  #============================
  
  if (show.another.banner){
    
    cat(paste0("\nPrepare new annotation bar...\n"))
    
    uniq.grps <- unique(lookup.table[,toupper(banner.name)])
    
    n1 <- length(uniq.grps)
    
    if (n1 <=4){
      # new.banner.col <- list$nice.cols.A[1:n1]
      # new.banner.col <- c("#C2E0D5","#ECCB80","#BE6698","#D18AB2")[1:n1]
      new.banner.col <- c("#E15759","#76B7B2","#F28E2B","#FF9DA7")[1:n1]
    } else {
      new.banner.col <-  distinctColorPalette(n1)
    } 
    
    names(new.banner.col) <- uniq.grps
    list.my.cols$new.banner.col <- new.banner.col
    
    show.annot.legend <- c(show.annot.legend, "TRUE")
    
    rm(n1)
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
  
  #=================================
  # if showing PATIENTS   ====
  #=================================  
  
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
 
  #========================================================
  # Prepare bottom heatmap bars (HeatmapAnnotation obj) ====
  #========================================================
  ###### This enforces to have at least ONE RIBBON + SURV.dots
  ###### update it in future if necessary
  #################################################################################################
  #### Define Top Annotation ====
  #################################################################################################
  
  if (show.survival) {
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
    
    
  }  else {
    cat(paste0("\nStart HeatMapAnnotation (line 570)...\n"))
    
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

  
  ###############################################################
  # == Set title/params and figure name ==== 
  ##############################################################
  cat(paste0("\nSet row/col orders...\n"))
  
  if (show.title){
    my.title <- paste0(title.str," \n# Alterations= ", nrow(data),"; # Genes with >= ", min.freq ," mutation(s)= ",length(unique(muts$GENE)),"; # Samples =", ncol(M))
  } else {
    my.title <- NULL    
  }
  
  if (is.null(patients.order)){
    column_order = NULL
  } else {
    column_order= as.character(patients.order)
  }
  
  if (is.null(muts.order) & is.null(cnvs.order) & is.null(svs.order)){
    row_order = NULL
  } else {
    row_order= c(muts.order,cnvs.order,svs.order)
  }
  
  ###############################################################
  # == Create a legend for multis if show.multis= TRUE ====
  ##############################################################
  if (show.multis){  
    
    ht.list = list(Legend(labels =  c(">1 Variant"),
                          labels_gp = gpar(fontsize = legend.label.font),
                          type = "points",
                          pch = 21,
                          size = unit(1, "cm"),
                          legend_gp = gpar(col = "black", fill= "#FAEFD1", lwd= 1, fontsize = legend.label.font),
                          background = NULL, grid_height = unit(1, "cm"),
                          grid_width = unit(1, "cm")))
  } else {
    ht.list = NULL
  }
  
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
  # == Generate Simple.ONCOPRINT ==== 
  ##############################################################
  cat(paste0("\nGenerating simple oncoprint...\n"))
  
  num.my.lgd.rows <- num.rows.heatmap.lgd 
  
  if (is.null(save.name)){
    saveFile <- file.path(savePath,paste0("Heatmap_TEMP_minFreq_",min.freq,".jpg"))
  } else {
    saveFile <- file.path(savePath,paste0("Heatmap_TEMP_minFreq_",min.freq,"_",save.name,".jpg"))
  }
  
  source(file.path("./sub_function/draw_basic_oncoprint.R"))
  
  simple.ht <- draw_basic_oncoprint(M, EFFECT, alter_fun, 
                                    
                                    saveFile= saveFile,
                                    
                                    list.colors= list, 
                                    
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
  
  
  if ((show.another.banner) | (show.response) | (show.individuals) ) { # col = list.my.cols
    
    cat(paste0("\nPrepare bottom annotation ...\n"))
    
    names(list.my.cols)[names(list.my.cols) == "new.banner.col"] <- toupper(banner.name)
    
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

    colnames(df)[colnames(df) == 'INDIVIDUAL.ID'] <- 'Patient.ID'
    names(list.my.cols)[names(list.my.cols) == 'INDIVIDUAL.ID'] <- 'Patient.ID'
     
    h2 = HeatmapAnnotation(df = df , name= "TEST", #df = data.frame(PATIENTS = pts), col= list(PATIENTS = col.assign), 
                           col = list.my.cols,
                           na_col = "grey",
                           simple_anno_size = unit(robbon.size, "cm"), # size of the ribbon
                           annotation_height =c(20,20), # this controls the height of the response/etc annotation that is added to the columns. However, in order to use mutiple features (e.g., response/celltype/etc) you have to use c(20,20,..) otherwise this generates error
                           gap = unit(c(5,5), "mm"), # this controls the gap between multiple annotation heatbars (for example, the space btw response and patient.id bars)
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
  } else {
    h2 = NULL # for example, you do not have any added bottom annotation but still like to see multis
  }
  
  # samples.order.mod <- colnames(simple.ht@matrix)
  
  samples.order.mod <- new.column_order
  
  ###############################################################
  ###############################################################
  # == Generate COMPLEX.ONCOPRINT ==== 
  ##############################################################
  ###############################################################
  ############################################################################################################################################
  # was trying to define legend as a list of list to allow separating CNVs and Mutations, but doesn't work. To be improved...
  ############################################################################################################################################
  # mutations.EFFECT.LAB <- c("Missense","Stop-gain","Frameshift indel","Inframe indel","Splicing variant","Complex Mutation","Unknown", "Other SNVs","Multiple variants")
  # cnvs.EFFECT.LAB <- c("CN Gain","CN Loss","CN-LOH","Deletion")
  # svs.EFFECT.LAB <- c("Inversion","Fusion","Translocation","Other SVs","Other SNVs","Inconclusive")
  # cyto.EEFECT.LAB <- c("Inversion","Fusion","Translocation","Multiple variants","Inconclusive")
  # if (!is.null(muts)){
  # 
  #   gg_list = list(
  #     mutations= list(title = "Mutations", at = EFFECT$variants[EFFECT$labels %in% mutations.EFFECT.LAB],
  #                               labels = EFFECT$labels[EFFECT$labels %in% mutations.EFFECT.LAB],
  #                               heatmap_legend_list= ht.list,
  #                               title_gp = gpar(fontsize = legend.title.font, fontface="bold"),
  #                               title_position = mut.legend.title.side,
  #                               # title_position= "topleft",
  #                               labels_gp = gpar(fontsize = legend.label.font),
  #                               grid_height= unit(1, "cm"),
  #                               nrow=num.rows.heatmap.lgd,
  #                               grid_width= unit(1, "cm"),
  #                               legend_height = unit(20, "cm")),
  #     cnvs= list(title = "CNVs", at = EFFECT$variants[EFFECT$labels %in% cnvs.EFFECT.LAB],
  #                labels = EFFECT$labels[EFFECT$labels %in% cnvs.EFFECT.LAB],
  #                heatmap_legend_list= ht.list,
  #                title_gp = gpar(fontsize = legend.title.font, fontface="bold"),
  #                title_position = mut.legend.title.side,
  #                # title_position= "topleft",
  #                labels_gp = gpar(fontsize = legend.label.font),
  #                grid_height= unit(1, "cm"),
  #                nrow=num.rows.heatmap.lgd,
  #                grid_width= unit(1, "cm"),
  #                legend_height = unit(20, "cm")))
  # }
  #########################################################################################################
  ###############################################################
  # == Generate COMPLEX.ONCOPRINT ==== 
  ##############################################################
  split <- rep(1:(ncol(M)/2), each = 2)
  
  cat(paste0("\nGenerate Final COMPLEX oncoprint ...\n"))
  
  ht <- oncoPrint(M, get_type = function(x) strsplit(x, ";")[[1]],
                  
                  alter_fun = alter_fun, col = append(list$mut.colors, list$cyto.colors),
                  
                  #axis_gp = gpar(fontsize = 8, fontface="bold"), # obsolete param
                  
                  column_order = samples.order.mod,
                  
                  row_order = row_order, #control the order of genes (rows)

                  remove_empty_columns = rem.empty,
                  
                  show_column_names = show.sample.names,
                  
                  column_split = split, 
                  
                  column_gap = unit(10, "mm"),
                  
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
                  column_names_max_height= unit(15,"cm") , # adjust this to control the name of samples (col names)
                  
                  # === Percent ====
                  pct_gp=gpar(fontsize = pct.font, fontface = "bold", col="black"), # specific control over percentage info on the left (add col="blue" to change colors)
                  row_names_gp = gpar(fontsize = rows.font, fontface="bold"), # gene-names and percent (if not prc_gp is defined above)
                  row_title_gp = gpar(fontsize =rows.font, col="blue",fontface = "bold"),
                  
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
  if (is.null(save.name)){
    saveFile <- file.path(savePath,paste0("Heatmap_minFreq_",min.freq,"_Complex_Heatmap.jpg"))
    
  } else {
    saveFile <- file.path(savePath,paste0("Heatmap_minFreq_",min.freq,"_Complex_Heatmap_",save.name,".jpg"))
  }
  
  jpeg(saveFile, width=w, height=h, pointsize =14, res = 100)
  
  if (heatmap.legend.side== annot.legend.side){
    draw(ht, split= LABS,  merge_legend = TRUE,  heatmap_legend_side = heatmap.legend.side, annotation_legend_side = heatmap.legend.side, annotation_legend_list = lgd_list,
               heatmap_legend_list = ht.list)
  } else {
    draw(ht, split= LABS,  merge_legend = FALSE,  heatmap_legend_side = heatmap.legend.side, annotation_legend_side = annot.legend.side, annotation_legend_list = lgd_list,
                 heatmap_legend_list = ht.list)
  }

  dev.off()
  
  # final.sample_order <- colnames(M)[column_order(ht)]
  
  
  ###############################################################################
  # === if automatic clustering is done, you can use the codes below 
  # ==== to decipher the exact order of clustered samples (add these to the calling code)
  # ==============================================================================
  col.list <- column_order(ht)
  htnames <- names(column_order(ht))
  col.orders <- col.list[[htnames[2]]]
  sample_order <- colnames(M)[col.orders]
  # =============================================================================
  
  # == ideas for future dev. 
  # draw(ht, padding = unit(c(40, 40), "mm")) 
  
  # == ideas for future dev. 
  # decorate_annotation("RESPONSE", {grid.text("value", unit(-2, "mm"), just = "right")})
  ###############################################################################
  
  cat(paste("\n\nThe file is saved at",saveFile,"\n"))
  
  return(list(ht.obj = ht,
              col.list= col.list,
              htnames= htnames,
              col.orders= col.orders,
              sample_order= sample_order
              # onco.samples= final.sample_order
              ))

}