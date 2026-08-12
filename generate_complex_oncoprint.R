generate_complex_oncoprint <-  function(muts= muts, cnvs= NULL, svs= NULL ,  # ******* define variants DFs [mut is required]
                                        
                                           cell.type.heatmap = NULL,
                                           
                                           cnvs.order= NULL, svs.order= NULL, muts.order= NULL, patients.order= NULL,   # ******* allows pre-defined orders
                                        
                                           sec.1.label = "MUT",
                                           sec.2.label = "CNVs",
                                           sec.3.label = "SVs", 
                                        
                                           highlight.events = NULL,
                                        
                                           surval.data= NULL, show.survival= FALSE, # currently under development
                                        
                                           show.blast= FALSE,
                                        
                                           show.response= FALSE, response.order= NULL, # ******* allows pre-defined orders
                                        
                                           show.another.banner=FALSE, banner.name= NULL, 

                                           banner.case.exceptions = NULL, # default exceptions set in prepare_BOTTOM_annotation : c("Patient.ID", "MRD_subtype", "CNV.WGS.CNVKIT.RHO", "RNA.EE", "Complex.Karyotype")
                                           
                                           custom.banner.colors = NULL, # named list of palettes for banners with no default in heatmap_colors(), e.g. list(TREATMENT_ARM= c("ArmA"="#1f77b4","ArmB"="#ff7f0e")); merged into list.ht.colors for this run only
                                           
                                           show.ALL= FALSE, ## added specifically for ALL prj. keep it as temp for other adaptations
                                        
                                           show.MPN= FALSE,
                                        
                                           show.purity= FALSE,
                                        
                                           show.individuals= FALSE, show.individuals.legend= FALSE,  
                                           
                                           lookup.table= NULL, # ******* pass lookup.table 
                                           
                                           show.sample.names = TRUE, show.border= FALSE, 
                                        
                                           show.multis= TRUE, # adds the dots to multi hits with multis.dot.size
                                           multi.col = "black",
                                        
                                           highlight.multis.cell = FALSE, # in addition to the dot also highlights the multis cells (currently set to pink, can change in define_ALTER_fun)
                                        
                                           rem.empty= TRUE, # ******* what params to show in legend?
                                           
                                           split.cols.by = NULL, 
                                        
                                           column_split= NULL,
                                           
                                           heatmap.legend.side= "right",
                                           
                                           mut.legend.title.side= "topleft", # ******* HINT: this can only be topleft/topcenter/ etc. otherwise error
                                           
                                           num.rows.heatmap.lgd= NULL, # ******* HEATMAP.legend 
                                           
                                           annot.legend.side= "bottom",  
                                           
                                           annot.title.side= "topleft", # *****annot.title.side param can only be topleft/ topcenter / leftcenter / lefttop / leftcenter-rot/ lefttop-rot
                                           
                                           num.rows.annot.lgd= NULL,  # ******* ANNOT.legend 
                                           
                                           min.freq= 1,
                                           
                                           include.these.events= NULL,
                                        
                                           add.default.subtitle= TRUE, 
                                        
                                           title.str= NULL, 
                                        
                                           save.path= NULL, # ******* title and save path
                                           
                                           save.name= NULL,
                                           
                                           cols.font= 25, rows.font= 25, pct.font= 20, 
                                           legend.label.font= 20, legend.title.font= 25, 
                                           fig.title.font= 28,  barplot.font= 25,  
                                           
                                           multis.dot.size = 0.8, #****FONTs: row.groupname.font is the same as rows.font
                                           
                                           right.w= 13, top.w= 8 , 
                                           
                                           ribbon.size= 1, # height (cm) of each bottom banner/annotation row
                                           
                                           w=50, h=50,  #**** Sizes of barplots and fig 
                                           
                                           axis.side= "left",
                                           top.annot.axis.side= "left",
                                           top.annotation_name_side = "left",
                                        
                                           banner.label.col= "#5b859e",
                                           legend.height = 20,
                                           na_col = "white"
                                        ){
  
  ## must be main branch
  
  
  library(randomcoloR)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(stringr)
  library(plyr)
  library(dplyr)
  
  graphics.off()
  
  #==================================================================================
  #' @title Generate a complex oncoprint
  #'
  #' @description A wrapper around the ComplexHeatmap package (\code{oncoPrint}/\code{Heatmap})
  #'   that facilitates building oncoprints without hand-writing the underlying ComplexHeatmap
  #'   calls each time. Plots MUTATIONs [required], CNVs [optional], and Structural Variants (SVs)
  #'   [optional] as a single oncoprint, with optional survival, response, subtype, and custom
  #'   bottom-banner annotations.
  #' @author Written by Noushin Farnoud (rahnaman@mskcc.org), Jul 2018. Last Update Aug 2026.
  #' @seealso example_Heatmap, test_required_fields
  #'
  #' @references Thanks to Zuguang Gu for the ComplexHeatmap package this function wraps:
  #'   \url{https://github.com/jokergoo/ComplexHeatmap}. Gu Z (2016). "Complex heatmaps reveal
  #'   patterns and correlations in multidimensional genomic data." Bioinformatics.
  #'
  #' @section Input Variant Data:
  #' @param muts Data frame of mutations (required). Required cols: TARGET_NAME, EFFECT (e.g. "missense"), GENE. Extra cols allowed.
  #' @param cnvs Data frame of CNVs. Required cols: TARGET_NAME, EFFECT (e.g. "AMP"/"DEL"), VAR_ID (e.g. "del(2q)"). Default NULL.
  #' @param svs Data frame of SVs. Required cols: TARGET_NAME, EFFECT (e.g. "fusion"), VAR_ID (row name to display). Default NULL.
  #'
  #' @section Cell-type Panel:
  #' @param cell.type.heatmap Optional cell-type/CIBERSORT matrix; if supplied, an extra heatmap panel is appended via add_cibersort_panel.R. Default NULL.
  #'
  #' @section Row/Column Order:
  #' @param cnvs.order Row order of CNV genes/events; default NULL uses clustering.
  #' @param svs.order Row order of SV genes/events; default NULL uses clustering.
  #' @param muts.order Row order of mutated genes; default NULL uses clustering.
  #' @param patients.order Column order of samples/patients; default NULL uses clustering.
  #'
  #' @section Section Labels:
  #' @param sec.1.label Row-split label for the mutation section. Default "MUT".
  #' @param sec.2.label Row-split label for the CNV section. Default "CNVs".
  #' @param sec.3.label Row-split label for the SV section. Default "SVs".
  #'
  #' @section Highlighting & Filtering Events:
  #' @param highlight.events Gene/event names to always keep in CNVs/SVs regardless of `min.freq`; only protects MUTs if `include.these.events` is also set (see NOTE).
  #' @param include.these.events Gene/event names to always keep in MUTs/CNVs/SVs regardless of `min.freq`.
  #' @param min.freq Only show GENEs with >= min.freq mutations. Default 1.
  #'
  #' @section Survival (under development):
  #' @param surval.data Optional survival data frame; enables the survival legend when combined with `show.survival`. Default NULL.
  #' @param show.survival Turn on the Dead/Alive survival legend. Default FALSE.
  #'
  #' @section Blast / MPN Display:
  #' @param show.blast Add a Blast (<10/>10 or <20/>=20) point legend to the top annotation. Default FALSE.
  #' @param show.MPN Turn on MPN-specific display options (blast points, complex karyotype, etc). Default FALSE.
  #'
  #' @section Response Annotation:
  #' @param show.response Add a bottom banner colored by `lookup.table$RESPONSE`. Default FALSE.
  #' @param response.order Order of levels shown in the response legend. Default NULL.
  #'
  #' @section Bottom Banner Annotations:
  #' @param show.another.banner Turn on to add extra bottom banner(s) named in `banner.name`. Default FALSE.
  #' @param banner.name Column name(s) in `lookup.table` to show as bottom banner(s), e.g. c("Gender","Final_subtype"). Default NULL.
  #' @param banner.case.exceptions Banner/column names to keep exactly as-is instead of Title Case (default set in prepare_BOTTOM_annotation.R). Default NULL.
  #' @param custom.banner.colors Named list of palettes for banners with no default in heatmap_colors(), e.g. list(TREATMENT_ARM = c("ArmA"="#1f77b4")); merged into the default colors for this run only. Default NULL.
  #'
  #' @section Project-specific Modes:
  #' @param show.ALL Turn on the ALL-project-specific banner set (subtype/gender/purity/EE) built in add_ALL_banners.R. Default FALSE.
  #' @param show.purity Add a purity annotation to the top barplot. Default FALSE.
  #'
  #' @section Patient/Individual Annotation:
  #' @param show.individuals Add a bottom banner grouping samples that belong to the same patient (useful for timeline datasets). Default FALSE.
  #' @param show.individuals.legend Show a legend for the patient banner (only used when `show.individuals` is TRUE). Default FALSE.
  #' @param lookup.table Sample-level metadata table (e.g. TARGET_NAME/RESPONSE/GENDER/...); required whenever any annotation bar is turned on. Default NULL.
  #'
  #' @section Display Features:
  #' @param show.sample.names Show sample names as column names. Default TRUE.
  #' @param show.border Draw a box around the frequency barplots. Default FALSE.
  #' @param show.multis If TRUE, a dot is drawn on gene-sample pairs with >1 variant. NOTE: this currently affects clustering. Default TRUE.
  #' @param multi.col Color of the multi-hit dot drawn when `show.multis` is TRUE. Default "black".
  #' @param highlight.multis.cell If TRUE, also highlights the fill of multi-hit cells (in addition to the dot); color set in define_ALTER_fun.R. Default FALSE.
  #' @param rem.empty Remove samples (columns) with no variants from the oncoprint. Default TRUE.
  #'
  #' @section Column Splitting:
  #' @param split.cols.by Column name in `lookup.table` used to compute an additional sample split/order factor (see prepare_COMPLEX_aes.R). Default NULL.
  #' @param column_split Pre-computed split vector/factor passed directly to `oncoPrint(column_split=)`. Default NULL.
  #'
  #' @section Main Heatmap Legend:
  #' @param heatmap.legend.side Side where the main mutation legend is displayed. Default "right".
  #' @param mut.legend.title.side Position of the mutation legend title; must be topleft/topcenter/etc. Default "topleft".
  #' @param num.rows.heatmap.lgd Number of rows for the mutation legend. Default NULL.
  #'
  #' @section Annotation Legend:
  #' @param annot.legend.side Side where the optional annotation bar legend(s) are displayed. Default "bottom".
  #' @param annot.title.side Position of the annotation legend title(s); one of topleft/topcenter/leftcenter/lefttop/leftcenter-rot/lefttop-rot. Default "topleft".
  #' @param num.rows.annot.lgd Number of rows for the annotation legend(s). Default NULL.
  #'
  #' @section Title & Save Path:
  #' @param title.str Your custom title text (or NULL for none). Always used as the base/prefix of the figure's `column_title` - this is NOT overridden by `add.default.subtitle`.
  #' @param add.default.subtitle If TRUE, appends an auto-generated stats summary (total alterations, gene/sample counts, and `min.freq` if `show.min.freq=TRUE`) after `title.str`. If FALSE, the title is exactly `title.str` with nothing appended. Default TRUE.
  #' @param save.path Output directory for the saved oncoprint. Defaults to the current working directory if NULL.
  #' @param save.name Base file name for the saved oncoprint (passed to initialize_data.R). Default NULL.
  #'
  #' @section Font Sizes:
  #' @param cols.font Sample name (column) font size. Default 25.
  #' @param rows.font Gene/CNV/SV name (row) font size. Default 25.
  #' @param pct.font Font size of the left-side percentage frequency labels. Default 20.
  #' @param legend.label.font Font size of legend labels (heatmap + annotation legends). Default 20.
  #' @param legend.title.font Font size of legend titles. Default 25.
  #' @param fig.title.font Font size of the figure title. Default 28.
  #' @param barplot.font Font size of the top/right frequency barplot axis ticks. Default 25.
  #' @param multis.dot.size Size of the multi-hit dot (same scale used for `rows.font`). Default 0.8.
  #'
  #' @section Figure & Barplot Size:
  #' @param right.w Width (cm) of the right gene-frequency barplot area. Default 13.
  #' @param top.w Height (cm) of the top sample-frequency barplot area. Default 8.
  #' @param ribbon.size Height (cm) of each bottom banner/annotation row. Default 1.
  #' @param w Width (px, at res=300) of the saved figure. Default 50.
  #' @param h Height (px, at res=300) of the saved figure. Default 50.
  #'
  #' @section Top Annotation Axis:
  #' @param axis.side Side of the main oncoprint frequency-axis ticks. Default "left".
  #' @param top.annot.axis.side Side of the top-annotation barplot axis ticks. Default "left".
  #' @param top.annotation_name_side Side where top-annotation names are drawn. Default "left".
  #'
  #' @section Bottom Banner Appearance:
  #' @param banner.label.col Text color for the bottom banner name labels. Default "#5b859e".
  #' @param legend.height Height (cm) of the annotation legend(s). Default 20.
  #' @param na_col Fill color for missing/NA values in the bottom banners. Default "white".
  #'
  #' @return A list with:
  #'   \describe{
  #'     \item{ht.obj}{The drawn ComplexHeatmap/oncoPrint object (or cibersort-combined object).}
  #'     \item{onco.samples}{Final column (sample) order used in the plot.}
  #'     \item{onco.genes}{Final row (gene/event) order used in the plot.}
  #'     \item{Fig.Path}{Path to the saved oncoprint PDF.}
  #'     \item{mut_matrix}{The alteration matrix (M) used to build the oncoprint.}
  #'     \item{draw.specs}{List of draw() params (LABS, legend sides/lists) to re-draw outside this function.}
  #'   }
  #==================================================================================
  
  # browser()
  
  suppressMessages(library("argparse", quietly = TRUE))
  
  if(!is.data.frame(muts)) {muts= as.data.frame(muts)}
  if (!is.null(cnvs) & !is.data.frame(cnvs)) {cnvs = as.data.frame(cnvs)}
  if (!is.null(svs) & !is.data.frame(svs)) {svs = as.data.frame(svs)}
  
  if (is.null(save.path)){
    save.path <- getwd()
    message(paste0("\n ***** NOTE: You did not pass 'save.path' param when calling the function. The default path used to save the generated oncoprints is --> ", save.path,"\n\n"))
  }
  
  dir.create(file.path(save.path,"TEMP"), showWarnings=FALSE)
  
  my.params = as.list(match.call(expand.dots=FALSE))
  
  RD <- met.brewer("Redon", type = "discrete")
  
  ###############################################################
  # == Test Required cols and contents  ====
  ##############################################################
  
  source(file.path("./sub_function/test_required_fields.R"))
  
  rename_IDs <- test_required_fields(muts= muts,  svs=svs, cnvs=cnvs, show.another.banner= show.another.banner, banner.name= banner.name, show.response= show.response, 
                                     split.cols.by= split.cols.by, show.individuals= show.individuals, lookup.table= lookup.table, annot.title.side= annot.title.side,
                                     show.ALL= show.ALL)
  
  muts <- rename_IDs$muts
  cnvs <- rename_IDs$cnvs
  svs <- rename_IDs$svs
  lookup.table <- rename_IDs$lookup.table
  REQ.cols <- rename_IDs$required.cols.lookup
  
  # cols.font <- as.numeric(cols.font)
  # browser()
  
  ############################################################
  # == Find a subset of Mutations that have >= min.freq variants
  ############################################################

  highlight.genes <- setdiff(highlight.events, c(svs$GENE, cnvs$GENE))
  
  include.these.genes <- setdiff(c(include.these.events, highlight.genes), c(svs$GENE, cnvs$GENE))

  if (!is.null(include.these.events)) {
    muts <- muts %>% group_by(GENE) %>% mutate(gene.freq= n()) %>% filter(gene.freq>= min.freq | GENE %in% include.these.genes) 
  } else {
    muts <- muts %>% group_by(GENE) %>% mutate(gene.freq= n()) %>% filter(gene.freq>= min.freq) 
  }
  
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
    
    highlight.svs <- setdiff(highlight.events, c(muts$GENE, cnvs$GENE))
    
    include.these.svs <- setdiff(c(include.these.events, highlight.svs), c(muts$GENE, cnvs$GENE))
    
    svs.test <- svs %>% group_by(GENE) %>% mutate(gene.freq= n()) %>% filter(gene.freq >= min.freq | GENE %in% include.these.svs) 
  
    if (nrow(svs.test)==0){
      cat(paste0("\n >>>>>>>  There are no SVs with min.freq you specified; So, plotting svs with at least 1 hit instead..."))
      svs <- svs %>% group_by(GENE) %>% mutate(gene.freq= n()) %>% filter(gene.freq>= 1) 
    } else {svs <- svs.test}
  
  rm(svs.test)
  
  }
  
  # CNVs ---
  #=====================
  if (!(is.null(cnvs))){
    
    highlight.cnvs <- setdiff(highlight.events, c(muts$GENE, svs$GENE))
    
    include.these.cnvs <- setdiff(c(include.these.events, highlight.cnvs), c(muts$GENE, svs$GENE))
    
    cnvs.test <- cnvs %>% group_by(GENE) %>% mutate(gene.freq= n()) %>% filter(gene.freq >= min.freq | GENE %in% include.these.cnvs) 
  
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
  Init.List <- initialize_data(data, muts, cnvs, svs, muts.order, cnvs.order, svs.order, min.freq,
                               sec.1.label=  sec.1.label , sec.2.label= sec.2.label, sec.3.label= sec.3.label , 
                               lookup.table, save.path, save.name= save.name)
  
  saveFile.1 <-  Init.List$saveFile.1
  saveFile.2 <-  Init.List$saveFile.2
  saveFile.2 <- gsub("png", "pdf", saveFile.2)
  
  my.fonts <-  Init.List$font.obj
  
  data <- Init.List$data
  muts <-  Init.List$muts
  
  gene.list= Init.List$gene.list
  
  ####################################
  # Load colors  ====
  ####################################
  
  # browser()
  
  cat(paste0("\nLoading default oncopring colors...\n"))
  
  source(file.path("./sub_function/heatmap_colors.R"))
  list.ht.colors <- heatmap_colors()
  
  if (!is.null(custom.banner.colors)) {
    list.ht.colors <- modifyList(list.ht.colors, custom.banner.colors)
  }
  
  source(file.path("./sub_function/match_default_colors.R"))
  
  # test colors if not shown properly
  #----------------------------------------------
  # source("./sub_function/color_alpha_test.R")
  # color_alpha_test(list.my.cols$response.colors)
  
  ###############################################################
  # == Adjust EFFECT to uniform texts  ====
  ##############################################################
  
  valid.effects <- tolower(rename_IDs$valid.effects)
  
  source(file.path("./sub_function/make_uniform_EFFECT_values.R"))
  data <- make_uniform_EFFECT_values(data) 
  
  ###############################################################
  # == Prepare the Heatmap rows and columns  ====
  ##############################################################
  
  # data <- unique(data) # this would eliminate multis for gene mutations of the same type in a sample
  
  # == GENES is the FINAL order of genes that will be the row names
  GENES <- data.frame(genes= gene.list$GENE, EFFECT=gene.list$LAB)
  
  # BE CAREFUL do not unique data : you will loose cases where a gene has multiple variants in the same patient 
  
  # SAMPLES = as.data.frame(with(data, table(TARGET_NAME)),stringsAsFactors = FALSE)
  
  if (!is.null(lookup.table)) {
    SAMPLES <- data.frame(TARGET_NAME = unique(as.character(lookup.table$TARGET_NAME)), stringsAsFactors = FALSE)
    
    missing.from.lookup <- setdiff(unique(as.character(data$TARGET_NAME)), SAMPLES$TARGET_NAME)
    
    if (length(missing.from.lookup) > 0) {
      stop(paste0("\n*** These TARGET_NAME(s) have variants in your data but are missing from lookup.table: '",
                  paste(missing.from.lookup, collapse = "', '"), "'. Fix lookup.table before proceeding.\n"))
    }
  } else {
    SAMPLES <- as.data.frame(with(data, table(TARGET_NAME)),stringsAsFactors = FALSE)
  }
  
  ###############################################################
  # == Prepare M and populate matrix of variants  ====
  ##############################################################
  
  source(file.path("./sub_function/prepare_fill_M.R"))
  M.List <- prepare_fill_M(data, SAMPLES$TARGET_NAME, GENES$genes, remove.empty.cols = rem.empty, show.multis = show.multis)
  
  M <- M.List$M
  gene.order <- M.List$gene.order
  events <- M.List$events
  
  ###############################################################
  # == Define "alter_fun" =====
  ##############################################################
  
  cat(paste0("\nLoading Default ALTER func...\n"))
  
  source(file.path("./sub_function/define_ALTER_fun.R"))
  alter_fun <- define_ALTER_fun(list.ht.colors, multis.dot.size, multi.col= multi.col, pink.multi = highlight.multis.cell)
  
  ###############################################################
  # == Define Labels for MUT/CNV/... segments  =====
  ##############################################################
  
  EFFECT.all <- list(variants = c("biallelic", "missense","stop_gain","frameshift_indel", "frameshift",
                                  "inframe_indel","inframe","splice_site_variant", "splicing",
                                  "initiator_codon_change",
                                  
                                  "complex", "complex_karyotype", "truncating",
                                  "unknown", 
                                  "amp", "gain",    
                                  "del", "loss",
                                  "loh",  "cnloh",
                                  "inv",   "INV",   
                                  "iso", "ISO",
                                  "fusion", "TRA",  
                                  "trans", "other_svs","tdup","dup","rearr",
                                  "add","der",
                                  "other_snvs",
                                  "other_cnvs",
                                  "unavailable","normal","karyotypic_abnormal"), 
                     
                     labels= c("biallelic","missense","stop_gain","frameshift_indel", "frameshift",
                               "Inframe indel","Inframe","splicing variant","splicing",
                               "Initiator_codon change",
     
                               "complex", "Complex karyotype", "truncating",
                               "Unknown",   
                               "Amplification", "GAIN",    
                               "Deletion", "LOSS",
                               "cnLOH", "cnLOH",
                               "Inversion", "INV",
                               "ISO","ISO",
                               "FUS", "TRA",
                               "TRA","Other SVs","Tandem duplication", "Duplication","Rearrangement",
                               "Add.","Der.",
                               "Other mutations",
                               "Other CN alterations",
                               "Unavailable","Normal","Karyotypic abnormal"))
  
  EFFECT <- list (variants = EFFECT.all[[1]][EFFECT.all[[1]] %in% data$EFFECT],
                  labels = EFFECT.all[[2]][EFFECT.all[[1]] %in% data$EFFECT]
  )
  

  LABS <- factor(gene.list$LAB, levels=c(sec.1.label, sec.2.label, sec.3.label))
  
  #################################
  # == Top-annotation (1)  ====
  #################################
  
  cat(paste0("\nPrepare Top Annotation...\n"))
  
  # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual' & brewer.pal.info$colorblind==TRUE,]
  # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  df <-  data.frame(TARGET_NAME= colnames(M)) #### IMPORTANT :: here we make sure the order of df = colnames of M (this guarantees the correct order of annotation)
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
    
    banner.name = unique(c("RESPONSE", banner.name))
    
    show.annot.legend <- c(show.annot.legend, "TRUE")
    
  } 
  
  #============================
  # if showing NEW.BANNER ====   
  #============================
  # browser()
  
  if (show.another.banner){
    
    if (show.ALL) {
      # browser()
      source(file.path("./sub_function/add_ALL_banners.R"))
      list.ALL.banners <- add_ALL_banners(list.my.cols, show.annot.legend, list.ht.colors, lookup.table)
      
      list.my.cols <- list.ALL.banners$updated.list.my.cols
      show.annot.legend <- list.ALL.banners$updated.show.annot.legend
      
      # banner.name <- c(banner.name, c("CNV.WGS.CNVKIT.RHO","RNA.EE","CNV.WGS.ACE.RHO"))

    } else {
      
      source(file.path("./sub_function/add_new_banner.R"))
      BannerList <- add_new_banner(banner.name, lookup.table, list.my.cols, show.annot.legend)
      list.my.cols= BannerList$list.my.cols
      new.banner.col = BannerList$new.banner.col
      show.annot.legend= BannerList$show.annot.legend
    }

  }
  
  banner.name = toupper(banner.name)
  
  # Case/punctuation/suffix-insensitive match of banner names to their default palette,
  # plus alias lookup for banners that share another palette's name (see banner_palette_aliases)
  for (banner_name in banner.name) {
    pal <- resolve_banner_palette(banner_name, list.ht.colors)
    if (!is.null(pal)) {
      list.my.cols[[banner_name]] <- pal
    }
  }
  
  ####################################
  # == Prepare the survival data ====
  ####################################
  
  if (show.ALL){
    
    if (!"CNV.WGS.CNVKIT.RHO" %in% colnames(df)) {
      df <- merge(df, lookup.table %>% dplyr::select(TARGET_NAME, CNV.WGS.CNVKIT.RHO, EXPRESSION.PROFILING.EFFICIENCY), by=c("TARGET_NAME"), all.x= TRUE)
      
    }
    
    df <- df %>% dplyr::mutate(CNV.WGS.CNVKIT.RHO= ifelse(CNV.WGS.CNVKIT.RHO=="#N/A", NA, as.numeric(CNV.WGS.CNVKIT.RHO)*100),
                                             RNA.EE= ifelse(EXPRESSION.PROFILING.EFFICIENCY=="#N/A", NA, as.numeric(EXPRESSION.PROFILING.EFFICIENCY)*100)
                                             )
    
    # surv.df$status.col= ifelse(surv.df$Death.status=="1","red","blue")
    # 
    # surv.df$pch= ifelse(surv.df$Death.status=="1",13,16)
    
    # show.annot.legend <- c(show.annot.legend, "TRUE")
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
    
    indiv.col <-  distinctColorPalette(n2)
    
    names(indiv.col) <- unique(lookup.table$INDIVIDUAL.ID)
    list.my.cols$INDIVIDUAL.ID <- indiv.col
    
    banner.name = unique(c(banner.name, "INDIVIDUAL.ID"))
    
    show.annot.legend <- c(show.annot.legend, "TRUE")
    
    rm(n2)
  } 
  
  ################################################################
  # Sort lookup table with the same order of M cols (TARGET_NAMEs)
  ################################################################

  if (!is.null(lookup.table)){
    rownames(lookup.table) <- lookup.table$TARGET_NAME
    lookup.table <- lookup.table[colnames(M),]
  }
  
  ################################################################
  #### Define Top Annotation ====
  ################################################################
  ###### This enforces to have at least ONE RIBBON + SURV.dots
  ###### update it in future if necessary
  #========================================================
  # browser()
  
  source(file.path("./sub_function/prepare_TOP_annotation.R"))

  h1 <- prepare_TOP_annotation(list.colors, show.border, 
                               barplot.font, legend.title.font, 
                               top.w, 
                               
                               show.purity= show.purity, 
                               purity.df= df, 
                               
                               show.blast= show.blast, 
                               show.MPN= show.MPN, 
                               show.ALL= show.ALL,
                               
                               lookup.table= lookup.table,
                               axis.side= axis.side,
                               top.annot.axis.side= top.annot.axis.side,
                               top.annotation_name_side = top.annotation_name_side,
                               banner.label.col= banner.label.col)
  
  # col_fun = colorRamp2(c(0, 50, 100), c("blue", "white", "red"))
  # ha = HeatmapAnnotation(foo = purity.df %>% pull(CNV.WGS.CNVKIT.RHO), col = list(foo = col_fun))
  
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
  
  # ###############################################################
  # # == Create a legend for BLAST if show.MPN= TRUE ====  <<<< must be done here but make sure u use the same thresh
  # ##############################################################
  
  blast_col <- c("<20"=DM[6], "≥20"=FH[1])
  # blast_col <- add.alpha(blast_col, alpha = 0.6)
  
  # browser()
  
  if (show.blast|show.MPN) {
    lgd_blast <- Legend(
      labels = c("<20", "≥20"),
      labels_gp = gpar(fontsize = legend.label.font),
      title_gp  = gpar(fontsize = legend.title.font, fontface = "bold"),
      title = "Blast",
      title_position = "topleft",
      type = "points",
      pch  = c(16, 16),
      size = unit(0.5, "cm"),
      background = "white",   # soft white
      nrow = 2,
      grid_width  = unit(1, "cm"),
      grid_height = unit(1, "cm"),
      legend_gp = gpar(
        col = c(blast_col[1], blast_col[2]),
        fontsize = legend.label.font
      )
    )
    
    lgd_list <- c(lgd_list, list(lgd_blast))
  }
  
  ###############################################################
  # == Create a legend for multis if show.multis= TRUE ==== MUST BE AFTER h1
  ##############################################################
  
  if (show.multis){  
    
    ht.list = list(Legend(labels =  c(">1 variant"),
                          labels_gp = gpar(fontsize = legend.label.font),
                          type = "points",
                          pch = 16, #21
                          size = unit(0.5, "cm"),
                          legend_gp = gpar(col = multi.col, fill= multi.col, lwd= 1, fontsize = legend.label.font), #FAEFD1
                          background = NULL, grid_height = unit(1, "cm"),
                          grid_width = unit(1, "cm")))
  } else {
    ht.list = NULL
  }
  
  ###############################################################
  # == Set title/params and figure name ==== 
  ##############################################################
  
  cat(paste0("\nSet row/col orders...\n"))
  
  cat(paste("\n **** You chose min.freq to filter all events = ", min.freq))

  # browser()
  
  if (add.default.subtitle){
      # my.title <- paste0(title.str," \n# Alterations= ", nrow(data),"; # Genes with >= ", min.freq ," mutations = ",length(unique(muts$GENE)),"; # Samples =", ncol(M))
      my.title <- paste0(title.str," \nTotal Alterations: ", nrow(data),", Unique Events: ", length(unique(muts$GENE)),", Total Samples: ", ncol(M),", # Min.Freq >= ",  min.freq)
    # } else{
    #   my.title <- paste0(title.str," \n# Alterations= ", nrow(data),"; # Top Genes = ",length(unique(muts$GENE)),"; # Samples=", ncol(M))
    # }
  } else {
    my.title <- title.str    
  }
  
  ###############################################################
  # == Set Col/sample order  ==== 
  ##############################################################
  # browser()
  
  if (is.null(patients.order)){
    column_order = NULL
  } else {
    column_order= as.character(patients.order)
  }
  
  ###############################################################
  # == Set Row (muts/cnvs/svs) order  ==== 
  ##############################################################
  
  # Set default orders from dataframes if specific orders are not provided
  muts.order.new <- if (is.null(muts.order)) as.character(unique(muts$GENE)) else muts.order
  cnvs.order.new <- if (is.null(cnvs.order)) as.character(unique(cnvs$GENE)) else cnvs.order
  svs.order.new <- if (is.null(svs.order)) as.character(unique(svs$GENE)) else svs.order
  
  # Create row_order based on provided and default orders
  row_order <- if (is.null(muts.order) && is.null(svs.order) && is.null(cnvs.order)) {
    NULL
  } else {
    row_order <- c(muts.order.new, cnvs.order.new, svs.order.new)
  }
  
  ###############################################################
  # == Generate Simple.ONCOPRINT ==== 
  ##############################################################

  # browser()
  
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
  ##=================================================================================
  ## *** IMPORTANT: Get the sample.order of simple.ht to sort the annotation, UNLESS
  ##                the user has strict patient order in input
  ##=================================================================================
  
  # browser()
  
  if (is.null(patients.order)){
    new.column_order <- colnames(M)[column_order(simple.ht)] #this is the order of the simple oncoprint with basic clustering
    
  } else {
    new.column_order <- patients.order
  }
  
  M <- M[, new.column_order]
  row.names(df) = df$TARGET_NAME
  df <- df[new.column_order,]
  # my.temp.column_order <- colnames(simple.ht@matrix)
  
  ###############################################################
  # == Prepare Heatmap Annotation/Aesthetics ==== 
  ##############################################################
  # browser()
  
  source(file.path("./sub_function/prepare_COMPLEX_aes.R"))
  
  complex.Annot <- prepare_COMPLEX_aes(data= data, M= M, highlight.events=highlight.events, 
                                       df= df, 
                                       list.my.cols= list.my.cols, 
                                      show.multis= show.multis, 
                                      show.another.banner= show.another.banner, 
                                      show.response = show.response, 
                                      show.individuals= show.individuals, 
                                      show.individuals.legend= show.individuals.legend,
                                      legend.title.font= legend.title.font, 
                                      legend.label.font= legend.label.font, 
                                      annot.title.side= annot.title.side, 
                                      num.rows.annot.lgd= num.rows.annot.lgd, 
                                      show.annot.legend= show.annot.legend, 
                                      ribbon.size= ribbon.size, 
                                      banner.name= banner.name, 
                                      rows.font= rows.font,
                                      split.cols.by= split.cols.by,
                                      show.ALL= show.ALL,
                                      show.MPN = show.MPN,
                                      banner.label.col= banner.label.col,
                                      legend.height = legend.height,
                                      banner.case.exceptions= banner.case.exceptions,
                                      na_col = na_col)
    
  #################################################################################################
  #################################################################################################
  #### Start Complex plot. ====
  #################################################################################################
  #################################################################################################
  # This works if multi hits are added as a label after initial set up however! I noticed it may not be good, so I moved it to M creation. You can change later if you want
  # if (show.multis){
  #   cat(paste0("\nStart multi.hit Oncoprint preparation...\n"))
  #   
  #   multi.hits <- data %>% dplyr::group_by(TARGET_NAME, GENE) %>% dplyr::mutate(N= n()) %>% dplyr::filter(N>1) %>% dplyr::select(TARGET_NAME, GENE) %>% unique()
  #   
  #   multi.hits <- data.frame(multi.hits)
  #   
  #   if (nrow(multi.hits) > 0) {
  #     for (k in seq_len(nrow(multi.hits))) {
  #       gene <- as.character(multi.hits$GENE[k])
  #       sample <- as.character(multi.hits$TARGET_NAME[k])
  #       val <- M[gene, sample]
  #       
  #       if (!grepl("biallelic", val)) {
  #         M[gene, sample] <- if (nchar(val) > 0) paste0(val, ";multi_hit") else "multi_hit"
  #       }
  #     }
  #   }
  # }
  
  # samples.order.mod <- colnames(simple.ht@matrix)
  
  # browser()
  
  ###############################################################
  # == If you have defined column_split ==== 
  ##############################################################
  
  if (!is.null(column_split)){
    column_split <- factor(column_split, levels = unique(column_split))
  }
  ###############################################################
  
  cat(paste0("\nGenerate Final COMPLEX oncoprint ...\n"))
  
  ht <- oncoPrint(M, get_type = function(x) strsplit(x, ";")[[1]],
                  
                  name= "oncoplot",
                  
                  # cluster_columns= col_hclust,
                  
                  alter_fun = alter_fun, col = append(list.ht.colors$mut.colors, list.ht.colors$cyto.colors),
                  
                  #axis_gp = gpar(fontsize = 8, fontface="bold"), # obsolete param
                  
                  column_order = new.column_order,
                  
                  row_order = row_order, #control the order of genes (rows)
                  
                  row_split = LABS,
                  
                  remove_empty_columns = rem.empty,
                  
                  show_column_names = show.sample.names,

                  column_gap = unit(5, "mm"),
                  
                  # === Gene barplots on the left ====
                  
                  column_split= column_split, # this supposed to add a vertical gap between columns based on a selected characteristic of samples (SPLIT col in lookup-table)
                  
                  # row_gap = unit(5, "mm"),
                  
                  bottom_annotation= complex.Annot$BotAnnot,
                  top_annotation = h1,
                  
                  left_annotation= complex.Annot$rowAnno,
                  
                  right_annotation = rowAnnotation(row_bar = anno_oncoprint_barplot(type= NULL,
                                                                                    border= show.border,
                                                                                    axis_param = list(side= "top",
                                                                                                      gp= gpar(fontsize= barplot.font, fontface="bold"))),
                                                   annotation_width= unit(right.w,"cm")),   ## controls the width of the row.barplots
                  
                  
                  #show_row_barplot = TRUE, # obsolete param
                  #row_barplot_width = unit(right.w, "cm"), # obsolete param
                  
                  row_names_side = "left", 
                  show_row_names = FALSE,
                  
                  pct_side = "right", # pct_digits = 2,
                  
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
                  
                  # ===  Percent/Rows/Genes ====
                  pct_gp=gpar(fontsize = pct.font, fontface = "bold", col="black"), # specific control over percentage info on the left (add col="blue" to change colors)
                  row_names_gp = gpar(fontsize = rows.font, fontface="bold"), # gene-names and percent (if not prc_gp is defined above)
                  row_title_gp = gpar(fontsize =rows.font+3, col=RD[3],fontface = "bold"), #blue
                  
                  # === Legend ====
                  # heatmap_legend_param = gg_list # list of list does not work here!
                  
                  heatmap_legend_param = list(title = "Alterations", at = EFFECT$variants,
                                              labels = EFFECT$labels,
                                              heatmap_legend_list= ht.list,
                                              title_gp = gpar(fontsize = legend.title.font, fontface="bold"),
                                              title_position = mut.legend.title.side,
                                              # title_position= "topleft",
                                              labels_gp = gpar(fontsize = legend.label.font),
                                              grid_height= unit(0.5, "cm"), # size of the mutation legend color-boxes
                                              nrow=num.rows.heatmap.lgd,
                                              grid_width= unit(0.5, "cm"),
                                              legend_height = unit(10, "cm"))
  ) 
  
  #======================================================
  # Draw  ====
  #======================================================
  
  # ht.2 <- ht  +  Heatmap(matrix(rnorm(nrow(M)*10), ncol = 10), name = "expr", width = unit(4, "cm"))
  # ht.2
  # draw(ht_list, row_split = sample(c("a", "b"), nrow(mat), replace = TRUE))
  
  ############################################
  # Cibersort plot ----
  ############################################
  
  source(file.path("./sub_function/prepare_TOP_annotation.R"))
  
  if (!is.null(cell.type.heatmap)){
    
    source(file.path("./sub_function/add_cibersort_panel.R"))
    
    hh <- add_cibersort_panel(ht, M, cell.type.heatmap, 
                              legend.title.font= legend.title.font,  legend.label.font= legend.label.font, rows.font= rows.font, 
                              annot.title.side= annot.title.side, show.sample.names= show.sample.names)
    } 
  else {hh <- ht} #bottom.heatmap.list$heat.1 %v%
    
  ###############################
  # Generate plots ----
  ###############################
  highlight_genes <- c("TP53", "SRSF2", "IDH2")
  
  # browser()
  
  cairo_pdf(saveFile.2, width = w / 2, height = h / 2)
  
  # png(saveFile.2, units="in", width = w / 2, height = h / 2, res = 300)
  
  # saveFile.2 <- gsub("png", "pdf", saveFile.2)
  # 
  # pdf(saveFile.2, width = w / 2, height = h / 2, useDingbats = FALSE)
  
  if (heatmap.legend.side== annot.legend.side){
    suppressMessages(draw(hh, split= LABS,  merge_legend = TRUE,  
                          heatmap_legend_side = heatmap.legend.side, 
                          annotation_legend_side = heatmap.legend.side, 
                          annotation_legend_list = lgd_list,
         heatmap_legend_list = ht.list))
    
    # This actually works ==== 
    # decorate_annotation("FINAL_SUBTYPE", {
    #   grid.text("Gender", x = unit(-2, "mm"), just = "right")
    # })
    # decorate_annotation("DNA_SUBTYPE", {
    #   grid.text("DNA", x = unit(-2, "mm"), just = "right")
    # })
    
  } else {  # ALL is this
    suppressMessages(draw(hh, split= LABS,  merge_legend = FALSE,  
                          heatmap_legend_side = heatmap.legend.side, 
                          annotation_legend_side = annot.legend.side, 
                          # annotation_legend_list = lgd_list,
                          annotation_legend_list = NULL,
                          # heatmap_legend_list = ht.list
                          heatmap_legend_list= c(ht.list, lgd_list)
                          )
                     )
  }
  
  # # Add horizontal lines at the gene positions
  # decorate_heatmap_body("oncoplot", {  # <--- match the name used in `oncoPrint()`
  #   gene_idx <- which(rownames(M) %in% c("TP53", "SRSF2", "IDH2"))
  #   ordered_idx <- row_order(ht)[["oncoplot"]][gene_idx]
  #   for (i in ordered_idx) {
  #     grid.lines(x = c(0, 1), y = unit(i, "native"), gp = gpar(col = "red", lwd = 2))
  #   }
  # })
  
  dev.off()
  
  cat(paste("\n *** Final oncoprint saved at: ",saveFile.2))
  
  # htShiny(hh, width1 = 1000)
  
  ############################################
  # Return specs so u can draw outside ----
  ############################################
  
  draw.specs <- list(LABS= LABS,
                     heatmap_legend_side= heatmap.legend.side,
                     annotation_legend_side= annot.legend.side,
                     annotation_legend_list= lgd_list,
                     heatmap_legend_list= ht.list)
  
  ############################################
  # Report final sample and gene/alt order ----
  ############################################
  # browser()
  if (!is.null(split.cols.by)){
    cat(paste("\n*** NOTE ***You can not get the final ordered list of samples (column_order) if you have chosen to split the columns by RESPONSE.\n 
              You can still get the list if you re-run the function and set split.by.response= FASLE. \n---> Future dev."))
    final.sample_order = NULL
    
    final.sample_order <- colnames(M)[column_order(ht)]
    final.row_order <- rownames(M)
    
  } else {
    
    co <- column_order(ht)
    
    if (is.list(co)) {
      co <- unlist(co, use.names = FALSE)
    }
    
    final.sample_order <- colnames(M)[co]
    final.row_order <- rownames(M)
    

  }

  #########################################################################################
  # === if automatic clustering is done, you can use the codes below 
  # ==== to decipher the exact order of clustered samples (add these to the calling code)
  # =======================================================================================
  # col.list <- column_order(ht)
  # htnames <- names(column_order(ht))
  # col.orders <- col.list[[htnames[2]]]
  # sample_order <- colnames(M)[col.orders]
  # =======================================================================================
  
  # == ideas for future dev. 
  # draw(ht, padding = unit(c(40, 40), "mm")) 
  
  # == ideas for future dev. 
  # decorate_annotation("RESPONSE", {grid.text("value", unit(-2, "mm"), just = "right")})
  #########################################################################################
  
  cat(paste("\n\nThe file is saved at",saveFile.2,"\n"))
  
  return(list(ht.obj = hh, 
              onco.samples= final.sample_order,
              onco.genes= final.row_order,
              Fig.Path = saveFile.2,
              mut_matrix= M,
              # rendered_ht= rendered_ht,
              draw.specs = draw.specs))
  
}