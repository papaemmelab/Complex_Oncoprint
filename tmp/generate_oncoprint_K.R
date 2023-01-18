generate_oncoprint_K <-  function(muts= muts, cnvs= NULL, svs= NULL ,  # define variants DFs
                                
                                cnvs.order= NULL, svs.order= NULL, muts.order= NULL, patients.order= NULL,   # define pre-set orders
                                
                                show.response= FALSE, show.another.banner=FALSE, banner.name= NULL, show.individuals= FALSE, show.individuals.legend= FALSE, # ******* user-defined ORDERs
                                
                                lookup.table= NULL, #******* pass lookup.table 
                                
                                show.sample.names = TRUE, show.border= FALSE, show.multis= FALSE, rem.empty= T, # ******* what params to show in legend?
                                
                                heatmap.legend.side= "right", mut.legend.title.side= "topleft", num.rows.heatmap.lgd= NULL, #******* HEATMAP.legend 
                                
                                annot.legend.side= "bottom", annot.title.side= "leftcenter", num.rows.annot.lgd= NULL,  #******* ANNOT.legend 
                                
                                min.freq= 1, show.title= TRUE, title.str= NULL, save.path= save.path, #******* title and save path
                                
                                cols.font= 20, rows.font= 18, pct.font= 16, legend.label.font= 10, legend.title.font= 14, fig.title.font= 18,  barplot.font= 10,  #****FONTs: row.groupname.font is the same as rows.font
                                 
                                right.w= 10, top.w= 10 , w=3000, h=2200,  #**** Sizes of barplots and fig 
                                
                                axis.side= "left"){

  
  library(randomcoloR)
  
  column_order = patients.order
  
  #==================================================================================
  #   Written by Noushin Farnoud, Jul 2018. Last Update Dec Sep 2019.  ====
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
  
  ###############################################################
  # == Test Required cols and contents  ====
  ##############################################################
  
  
  source(file.path("./sub_function/test_required_fields.R"))

  #rename_IDs <- test_required_fields(muts= muts,  svs=svs, cnvs=cnvs, show.another.banner= show.another.banner, banner.name= banner.name, show.response= show.response, show.individuals= show.individuals, lookup.table= lookup.table)
  
  # muts <- rename_IDs$muts
  # cnvs <- rename_IDs$cnvs
  # svs <- rename_IDs$svs
  # lookup.table <- rename_IDs$lookup.table
  # REQ.cols <- rename_IDs$required.cols.lookup
  
  REQ.cols <- c("TARGET_NAME","IDH")
  
  ####################################
  # Load colors  ====
  ####################################
  
  source(file.path("./sub_function/heatmap_colors_K.R"))
  list <- heatmap_colors_K()
  
  ############################################################
  # == Find a subset of Mutations that have >= min.freq variants
  ############################################################
  
  mut.freq <- ddply(muts, c("GENE"), "nrow",.drop = TRUE)
  
  muts <- merge(muts, mut.freq, by="GENE")
  
  muts <- subset(muts, nrow >= min.freq)
  
  ###############################################################
  # == Add MUTATIONS  ====
  ##############################################################
  
  source(file.path("./sub_function/sort_variants.R"))
  
  A <- sort_variants(muts, muts.order, gene.list = NULL, group.label= "Substitutions/Indels")
  
  data <- A[[1]][,c("TARGET_NAME","GENE","EFFECT")]
  
  gene.list <- A$gene.list
  
  ###############################################################
  # == Add Structural Variants  ====
  ##############################################################
  
  if (!is.null(svs)){
    
    source(file.path("./sub_function/sort_variants.R"))
    
    B <- sort_variants(svs, svs.order, gene.list,group.label="SVs")
    
    data <- rbind(data, B[[1]][,c("TARGET_NAME","GENE","EFFECT")])
    
    gene.list <- B$gene.list
    
  }
  
  ###############################################################
  # == Add CNVs  ====
  ##############################################################
  
  if (!is.null(cnvs)){
    
    source(file.path("./sub_function/sort_variants.R"))
    
    D <- sort_variants(cnvs, cnvs.order, gene.list, group.label= "CNVs")
    
    data <- rbind(data, D[[1]][,c("TARGET_NAME","GENE","EFFECT")])
    
    gene.list <- D$gene.list
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
  
  M[M==0] = ""
  
  ###############################################################
  # == Adjust EFFECT to uniform texts  ====
  ##############################################################
  
  data$EFFECT <- gsub("^frameshift_indel$|^frameshift_indel$|^frameshift_variant$","frameshift_indel",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^inframe_codon_loss$|^inframe_indel$|^inframe_deletion$|^inframe_del$|^inframe_codon_gain$|^inframe_ins$|^inframe_insersion$|^inframe_variant$","inframe_indel", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^complex_change_in_transcript$|^complex$","complex",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^stop_gained$|^stopgain_SNV$|^stopgain$|^stop_gain$","stop_gain", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^missense$|^non_synonymous_codon$|^nonsynonymous_SNV$","missense", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^splice_site_variant$|^initiator_codon_change$","splice_site_variant", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^unknown$","unknown",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^other_snvs$|^other$","other_snvs",ignore.case = TRUE, data$EFFECT)
  
  data$EFFECT <- gsub("^amp$|^amplification$|^gain$|^CN-gain$","AMP", ignore.case = TRUE,  data$EFFECT)
  data$EFFECT <- gsub("^del$|^deletion$|^loss$|^CN-del$", "DEL", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^LOH$","LOH", ignore.case = TRUE,  data$EFFECT)
  data$EFFECT <- gsub("^SNV$","SNV", ignore.case = TRUE,  data$EFFECT)
  data$EFFECT <- gsub("^CNA$","CNA", ignore.case = TRUE,  data$EFFECT)
  data$EFFECT <- gsub("^inv$|^inversion$", "INV", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^tandem duplications$|^tandem_duplications$|^tandem dup$","TDUP", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^fusion$|^fuse$|^fus$","FUS",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^translocation$|^trans$","TRANS",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^other_svs$","OTHER_SVs",ignore.case = TRUE, data$EFFECT)
  #data$EFFECT <- gsub("^NA$","",ignore.case = TRUE, data$EFFECT)
  ###############################################################
  # == Add the Event.Type in the Matrix ====
  
  
  ##############################################################
  
  M[M==0] = ""
  
  events <- factor(unique(data$EFFECT), levels=c("one","two","four","eight","none","SNV","CNA","unknown","other_snvs","missense","sub","indel","splice_site_variant","complex","stop_gain","inframe_indel","frameshift_indel","AMP","DEL","LOH","INV","FUS","TRANS","TDUP","other_svs"))
  events <- events[order(events)]
  
  for (i in 1: length(events)){
    
    temp <- subset(data, EFFECT==events[i])
    
      for (j in 1:nrow(temp)) {
        M[temp$GENE[j], temp$TARGET_NAME[j]] <- paste0(M[temp$GENE[j], temp$TARGET_NAME[j]], unique(temp$EFFECT),";", collapse = "")
      }
  }
  
  
  
  ###############################################################
  # == Identify multi-hit sample-gene pairs =====
  ##############################################################
  
  if (show.multis){
    multi.hits <- data %>% group_by(TARGET_NAME, GENE) %>% mutate(N= n()) %>% filter(N>1) %>% select(TARGET_NAME, GENE) %>% unique()

    multi.hits <- data.frame(multi.hits)

    for (k in 1: nrow(multi.hits)){
      M[as.character(multi.hits$GENE[k]), as.character(multi.hits$TARGET_NAME[k])] <- paste0(M[as.character(multi.hits$GENE[k]), as.character(multi.hits$TARGET_NAME[k])], "multi_hit",";", collapse = "")
    }
  }

  ###############################################################
  # == Define "alter_fun" =====
  ##############################################################
  
  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.52, "mm"), h-unit(0.52, "mm"), gp = gpar(fill = "#f0f0f0", col = NA)) # alpha=0.5
    },
    one = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$mut.colors[["one"]][1], col = NA))
    },
    two = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$mut.colors[["two"]][1], col = NA))
    },
    four = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$mut.colors[["four"]][1], col = NA))
    },
    eight = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$mut.colors[["eight"]][1], col = NA))
    },
    none = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$mut.colors[["none"]][1], col = NA))
    },
    unknown = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$mut.colors[["unknown"]][1], col = NA))
    },
    other_snvs = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$mut.colors[["other"]][1], col = NA))
    },
    missense = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$mut.colors[["missense"]][1], col = NA))
    },
    sub = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$mut.colors[["sub"]][1], col = NA))
    },
    indel = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$mut.colors[["indel"]][1], col = NA))
    },
    SNV = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$mut.colors[["SNV"]][1], col = NA))
    },
    splice_site_variant = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$mut.colors[["splice_site_variant"]][1], col = NA))
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
    AMP = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$cyto.colors[["AMP"]][1], col = NA))
    },
    DEL = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$cyto.colors[["DEL"]][1], col = NA))
    },
    LOH = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$cyto.colors[["LOH"]][1], col = NA))
    },
    CNA = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$cyto.colors[["CNA"]][1], col = NA))
    },
    INV = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list$cyto.colors[["INV"]][1], col = NA))
    },
    FUS = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h-unit(0.3, "mm"), gp = gpar(fill = list$cyto.colors[["FUS"]][1],  col = NA))
    },
    TRANS=function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h*0.33, gp = gpar(fill = list$cyto.colors[["TRANS"]][1],  col = NA))
    },
    other_svs=function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h*0.33, gp = gpar(fill = list$cyto.colors[["other"]][1],  col = NA))
    },
    multi_hit = function(x, y, w, h) {
      grid.points(x, y, pch = 21, size = unit(0.8, "cm"), gp = gpar(col = "black", fill= "#FAEFD1"))
    }
  )
  
  EFFECT <- list(variants = c("one","two","four","eight","none","SNV","CNA","sub","indel","missense","stop_gain","frameshift_indel","inframe_indel","splice_site_variant","complex","unknown",  "AMP",     "DEL",  "LOH",  "INV",      "FUS",   "TRANS",           "other_svs","other_snvs","OTHER_SVs","multi_hit"), 
                    labels= c("one","two","four","eight","none","SNV","CNA","sub","indel","Missense","Stop-gain","Frameshift indel","Inframe indel","Splicing variant","Complex",   "Unknown",  "CN Gain","CN Loss","LOH", "Inversion","Fusion","Translocation","Other SVs","Other SNVs","Other SVs","Multiple variants"))
  
  LABS <- factor(gene.list$LAB, levels=c("Substitutions/Indels","SVs","Cytogenetics","CNVs"))
  
  #################################
  # == Top-annotation (1)  ====
  #################################
  
  # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual' & brewer.pal.info$colorblind==TRUE,]
  # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  df <-  data.frame(TARGET_NAME= colnames(M))

  
  if (show.another.banner | show.response | show.individuals){
    
    df <- merge(df, lookup.table[,REQ.cols], by=c("TARGET_NAME"))
    
    }
 
  #################################################################################
  # == PREPARE HeatmapAnnotation Obj (Response, Sample.Source, Individual.ID) ==== 
  #################################################################################
 
  df$TARGET_NAME <- NULL

  #==================================================
  # == Create COLOR palletes for HeatmapAnnotation
  #==================================================
  
  list.my.cols <- list()
  
  show.annot.legend <- c()
  
  #============================
  # if showing RESPONSE ====  
  #============================
  
  if (show.response){
    
    resp.col <- list$response.colors[names(list$response.colors) %in% unique(lookup.table$RESPONSE)]
    
    list.my.cols$Response <- resp.col
    
    show.annot.legend <- c(show.annot.legend, "TRUE")
    
  } 
  
  #============================
  # if showing NEW.BANNER ====   
  #============================
  
  if (show.another.banner){
    
    uniq.grps <- unique(lookup.table[,toupper(banner.name)])
    
    n1 <- length(uniq.grps)
    
    if (n1 <9){
      new.banner.col <- list$nice.cols.A[1:n1]
    }
    
    if (n1>=9){
      new.banner.col <-  distinctColorPalette(n1)
    } 
    
    names(new.banner.col) <- uniq.grps
    list.my.cols$new.banner.col <- new.banner.col
    
    show.annot.legend <- c(show.annot.legend, "TRUE")
    
    rm(n1)
  } 
  
  #=================================
  # if showing PATIENTS   ====
  #=================================  
  
  if (show.individuals){
    
    n2 <- length(unique(lookup.table$INDIVIDUAL.ID))

    indiv.col <-  c("thistle4","olivedrab3","salmon") #distinctColorPalette(n2)

    names(indiv.col) <- unique(lookup.table$INDIVIDUAL.ID)
    list.my.cols$Group <- indiv.col
    
    show.annot.legend <- c(show.annot.legend, show.individuals.legend)
    
    rm(n2)
  } 
 
  #========================================================
  # Prepare bottom heatmap bars (HeatmapAnnotation obj) ====
  #========================================================
  
  if ((show.another.banner) | (show.response) | (show.individuals)){ # col = list.my.cols
    
    names(list.my.cols)[names(list.my.cols) == "new.banner.col"] <- banner.name
    
    colnames(df)[colnames(df) == "RESPONSE"] <- "Response"
    colnames(df)[colnames(df) == "INDIVIDUAL.ID"] <- "Group"
    colnames(df)[colnames(df) == toupper(banner.name)] <- banner.name
    
    # anno_oncoprint_barplot(type = NULL, which = c("column", "row"),
    #                        bar_width = 0.6, axis = TRUE,
    #                        axis_param = if(which == "column") default_axis_param("column") else list(side = "top", labels_rot = 0),
    #                        width = NULL, height = NULL, border = FALSE)
    
    #################################################################################################
    #### Define Annotaion ribbons for RESPONSE etc. ====
    #################################################################################################
    
    response.title.pos= annot.title.side
    
    h2 = HeatmapAnnotation(df = df , name= "TEST", #df = data.frame(PATIENTS = pts), col= list(PATIENTS = col.assign), 
                           # barplot1 = anno_oncoprint_barplot(), 
                           col = list.my.cols , 
                           #simple_anno_size = unit(2, "cm"), # size of the ribbon
                           annotation_height =c(100), # this controls the height of the response/etc annotation that is added to the columns. However, in order to use mutiple features (e.g., response/celltype/etc) you have to use c(20,20,..) otherwise this generates error
                           gap = unit(c(3,3), "mm"), # this controls the gap between multiple annotation heatbars (for example, the space btw response and patient.id bars)
                           show_annotation_name= T,
                           show_legend = as.logical(show.annot.legend),
                           #show_annotation_name= show.annot.legend, 
                           annotation_name_offset = unit(10, "mm"),
                           gp = gpar(col = "black"),
                           annotation_name_gp= gpar(fontsize = legend.title.font, fontface= "bold", col="black"),
                           annotation_legend_param = list(#title = legend.tit.df,
                                                          title_gp = gpar(fontsize = legend.title.font, fontface="bold"),
                                                          title_position = response.title.pos, 
                                                          labels_gp = gpar(fontsize = legend.label.font),
                                                          grid_height= unit(2, "cm"),
                                                          grid_width= unit(2, "cm"),
                                                          nrow= num.rows.annot.lgd,
                                                          legend_height = unit(2, "cm")
                                                          )
                           )
  } else {
    h2 = NULL
  }
  
  ###############################################################
  # == Set title/params and figure name ==== 
  ##############################################################
  
  mat_letters = matrix(sample(letters[1:4], 1242, replace = TRUE), 69)
  
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
  
  
  main.legend.pos <- mut.legend.title.side
  
  saveFile <- file.path(savePath,paste0("Heatmap_minFreq_",min.freq,".jpg"))
  
  jpeg(saveFile, width=w, height=h, pointsize =14, res = 100)
  
  ###############################################################
  # == Define heatmap ====
  ##############################################################

    ht <- oncoPrint(M, name="hmp",get_type = function(x) strsplit(x, ";")[[1]],
                    
                    alter_fun = alter_fun, col = append(list$mut.colors, list$cyto.colors),
                    
                    #axis_gp = gpar(fontsize = 8, fontface="bold"), # obsolete param

                    column_order = column_order,
                    row_order = row_order, #control the order of genes (rows)
                    
                    remove_empty_columns = rem.empty,
                    
                    show_column_names = show.sample.names,
                    
                    # === Gene barplots on the left ====
                    
                    bottom_annotation = h2, 
                    # top_annotation_height = unit(10, "cm"),
                    #top_annotation=NULL, uncomment to remove top barplot
                    #**top_annotation = HeatmapAnnotation(column_bar = anno_oncoprint_barplot(),annotation_height = unit(8, "cm")),
                    # row_barplot_width = unit(10, "cm"), ###to get rid of barplot change to 0
                    
                    #top_annotation = h1,
                    
                    
                    #top_annotation = HeatmapAnnotation(column_bar = anno_oncoprint_barplot(type= NULL,
                                                                                          # border= show.border, # do you want the top-barplot to have a border?
                                                                                          # axis_param = list(side = axis.side, 
                                                                                                             # side = "right", 
                                                                                                             #labels = c("zero", "half", "one"),
                                                                                                             # at = c(0, 0.5, 1), 
                                                                                                             # labels_rot = 45,
                                                                                                            # gp= gpar(fontsize= barplot.font, fontface="bold"))
                                                            #                               ),
                                                       #annotation_height = unit(top.w, "cm")
                                                       #),
                      
                   # right_annotation = rowAnnotation(row_bar = anno_oncoprint_barplot(type= NULL,
                                                                           #           border= show.border, 
                                                                          #            axis_param = list(side= "top", gp= gpar(fontsize= barplot.font, fontface="bold"))),
                                                  #   annotation_width= unit(right.w,"cm")),   ## controls the width of the row.barplots
                    
                    #show_row_barplot = TRUE, # obsolete param
                    #row_barplot_width = unit(right.w, "cm"), # obsolete param
                    
                    split= LABS,
                    
                    # ==========================================
                    # ==========================================
 
                    # === Title ====
                    column_title = my.title,
                    column_title_gp = gpar(fontsize = fig.title.font, fontface = "bold"), # title font-size
                    gap = unit(2, "mm"),
                    
                    # === Column/Sample names ====
                    column_names_gp = gpar(cex=1, col= "black", fontsize = cols.font, fontface="bold"), #default size = 18
                    column_names_max_height= unit(15,"cm") , # adjust this to control the name of samples (col names)
                    
                    # === Percent ====
                    pct_gp=gpar(fontsize = pct.font, fontface = "bold", col="black"), # specific control over percentage info on the left (add col="blue" to change colors)
                    row_names_gp = gpar(fontsize = rows.font, fontface="bold"), # gene-names and percent (if not prc_gp is defined above)
                    row_title_gp = gpar(fontsize =rows.font, col="blue",fontface = "bold"),
                      
                    # === Percent ====
                    # layer_fun = function(j, i, x, y, w, h, fill) {
                    #   ind_mat = restore_matrix(j, i, x, y)
                    #   ind = unique(c(ind_mat[2, ], ind_mat[, 3]))
                    #   grid.points(x[ind], y[ind], pch = 16, size = unit(4, "mm"))},
                   
                    # === Legend ====
                    
                    heatmap_legend_param = list(title = "Alterations", at = EFFECT$variants,
                                                labels = EFFECT$labels,
                                                title_gp = gpar(fontsize = legend.title.font, fontface="bold"),
                                                title_position = main.legend.pos, 
                                                labels_gp = gpar(fontsize = legend.label.font),
                                                grid_height= unit(2, "cm"),
                                                nrow=num.rows.heatmap.lgd,
                                                grid_width= unit(2, "cm"),
                                                legend_height = unit(2, "cm")))
  
               

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
                          background = NULL, grid_height = unit(2, "cm"),
                          grid_width = unit(2, "cm")))
  } else {
    ht.list = NULL
  }
  
  ###############################################################################
  # == Draw the oncoprint (merge legends only if there are on the same side) ====
  ###############################################################################
  
  if (heatmap.legend.side== annot.legend.side){
    draw(ht, split= LABS,  merge_legend = TRUE,  heatmap_legend_side = heatmap.legend.side, annotation_legend_side = heatmap.legend.side,
         heatmap_legend_list = ht.list)
  } else {
    draw(ht, split= LABS,  merge_legend = FALSE,  heatmap_legend_side = heatmap.legend.side, annotation_legend_side = annot.legend.side,
         heatmap_legend_list = ht.list)
    
    # draw(ht, split= LABS,  merge_legend = FALSE)
  }
  
  
 # decorate_heatmap_body("hmp", {
#    i = which(colnames(M) == "I-H-118729-T1")
#    x1 = 0.32 #i/ncol(M)
#    x2 = 0.326
#    x3= 0.596
#    grid.polyline(c(x2,x3,x2,x3), c(-1.97,-1.97,1,1), id=c(2,3,2,3), gp = gpar(lwd = 2, lty = 1))
#  })
  
  dev.off()

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

  cat(paste("\n\nThe file is saved at",saveFile,"\n"))
  
  return(list(ht.obj = ht))

}