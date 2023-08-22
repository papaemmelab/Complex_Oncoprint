test_required_fields  <- function(muts= muts, cnvs= NULL, svs= NULL, annot.title.side= annot.title.side,
                                  patients.order= NULL, muts.order= NULL, cnvs.order= NULL, svs.order= NULL,
                                  show.another.banner=FALSE, banner.name= NULL, 
                                  show.response= FALSE, split.cols.by = NULL, show.individuals= FALSE, lookup.table= NULL)
  {

  #==================================================================================
  #   Written by Noushin Farnoud, Jul 2018, Last Update Jan 2023  ----
  #----------------------------------------------------------------------------------
  #   This function tests the fields of MUTATIONS/CNVs/SVs dataframes that are
  #   mandatory field is missing.
  # 
  #   See also example_Heatmap, heatmap_colors, generate_oncoprint.
  #==================================================================================

  cat(paste("\nStarting testing required fields ...\n"))
  
  mut.VT.options <- c("non_synonymous_codon","missense",
                      
                      "missense_codon", # noush added  
                      
                      "stop_gained","stop_lost","stop_retained_variant",  #stop_lost can be from both indels and subs     
                      
                      "stop_gain", 
                      
                      "missense_mutation",
                      
                      "in_frame_del", "in_frame_ins",
                      
                      "nonsense_mutation",
                      
                      "frame_shift_del",
                      
                      "splice_site",
                      
                      "frame_shift_ins",
                      
                      "nonstop_mutation",
                      
                      "translation_start_site",
                      
                      # "nonsense",  # noush added
                      
                      "splice_site_variant","initiator_codon_change",
                      
                      "inframe_codon_loss", "inframe_codon_gain", "inframe_variant", 
                      
                      "inframe_indel","inframe_deletion", "inframe_insersion", # noush added   
                      
                      "complex_change_in_transcript", 
                      
                      "complex", # noush added
                      
                      "frameshift_variant",
                      
                      "frameshift_indel", "frameshift_del", # noush added  
                      
                      "frameshift_insertion",
                      
                      "other_snvs", # noush added 
                      
                      # == Noush Added for CNVs and SVs
                      
                      "amp","amplification","gain","CN-gain",
                      
                      "del","deletion","loss","CN-del",
                      
                      "LOH",
                      
                      "unknown",
                      
                      "inv","inversion",
                      
                      "tandem duplications","tandem_duplications","tandem dup", "tdup",
                      
                      "fusion","fus",
                      
                      "translocation","trans", "complex_karyotype", "normal", "normal_karyotype",
                      
                      "other_svs","other_cnvs","rearr","rearrangements","rearrangement",
                      
                      "N/E", "inconclusive","unavailable","N_E","N/A",
                      "karyotypic_abnormal"
  )
                      
  # stop_lost, splice_site_variant, stop_retained_variant, initiator_codon_change can be from both indels or subs in VEP
  
  required.cols <- c("TARGET_NAME","EFFECT")
  
  #=====================================================
  # *** HINT: first change all col names to uppercase 
  #=====================================================
  if (!is.null(lookup.table)){
    lookup.table <- data.frame(lookup.table)
  }
  #=================================================
  # == check field of MUTs and Mutations.VT ----
  #=================================================
  colnames(muts) <- toupper(colnames(muts))
  
  if ("EFFECT" %in% toupper(colnames(muts))){
        muts$EFFECT <- tolower(muts$EFFECT)
  } 
  
  if (!all(c(required.cols ,"GENE") %in% colnames(muts))) {
    stop(paste("\n****Required cols are missing from input MUTATION dataframe =>\t", setdiff(c(required.cols ,"GENE"), colnames(muts))),"\n")
  } else {
    cat(paste("\n**** test_required_fields TEST (1): successfull!\tAll required columns exist in MUTATION dataframe...\n"))
  }
  
  if (!all(tolower(muts$EFFECT) %in% tolower(mut.VT.options))){
    cat(paste0("\n\nERROR --- Valid options for EFFECT are : \n", paste(mut.VT.options,collapse  = ", ")))
    stop(paste("\n****Not valid entry for mutation EFFECT  => ", setdiff(muts$EFFECT, mut.VT.options),
               "\n\nCheck available mutation types in README, or change the troublesome EFFECT to 'other_SNVs' or 'other_SVs' to continue."))
  } else {
    cat(paste("\n**** test_required_fields TEST (2): successfull!\tAll mutations EFFECTs are consistent with available options (for coloring)...\n"))
  }
  
  #===============================
  # == Check INPUT params  ----
  #===============================
  
  if (show.another.banner | show.response | show.individuals){
    if (is.null(lookup.table)){
    stop(paste("\n****Required INPUT for added annotation is missing =>\t","lookup.table" ,"\n"))
    }}
  
  #===============================
  # == Check INPUT params ----
  #===============================

  if (!is.null(patients.order) & !is.vector(patients.order)){stop("input patients.order must be a simple vector.\n")}
  if (!is.null(muts.order) & !is.vector(muts.order)){stop("input muts.order must be a simple vector.\n")}
  if (!is.null(cnvs.order) & !is.vector(cnvs.order)){stop("input cnvs.order must be a simple vector.\n")}
  if (!is.null(svs.order) & !is.vector(svs.order)){stop("input svs.order must be a simple vector.\n")}
  
  #==============================================================================
  # == ComplexHeatmap.OncoPrint func accepted mut.legend.title.side values ----
  #==============================================================================
  
  if (!any((annot.title.side) %in% c("topleft", "topcenter", "leftcenter", "lefttop", "leftcenter-rot", "lefttop-rot"))){
    stop("annot.title.side param can only be topleft/ topcenter / leftcenter / lefttop / leftcenter-rot/ lefttop-rot")
  }
  
  
  #===============================
  # == check field of SVs  ----
  #===============================
  if (!is.null(svs)){
    colnames(svs) <- toupper(colnames(svs))
    
    if (!all(c(required.cols, "VAR_ID") %in% colnames(svs))) {
      stop(paste("\n****Required cols are missing from input SV dataframe  =>\t", setdiff(c(required.cols, "VAR_ID"), colnames(svs))),"\n") 
    } else {
      colnames(svs)[colnames(svs) == 'VAR_ID'] <- 'GENE' 
      cat(paste("\n**** test_required_fields TEST (3): successfull!\tAll required columns exist in SV dataframe...\n"))
    }}
  
  #===============================
  # == check field of CNVs  ----
  #===============================
  
  if (!is.null(cnvs)){
    
    colnames(cnvs) <- toupper(colnames(cnvs))
    
    if (!all(c(required.cols, "VAR_ID") %in% colnames(cnvs))) {
      stop(paste("\n****Required cols are missing from input CNV dataframe =>\t", setdiff(c(required.cols, "VAR_ID"), colnames(cnvs))),"\n") 
    } else {
      colnames(cnvs)[colnames(cnvs) == 'VAR_ID'] <- 'GENE'
      cat(paste("\n**** test_required_fields TEST (4): successfull!\tAll required columns exist in CNV dataframe...\n"))
    }}
  

  #=========================================================
  # == check field of SAMPLE.SOURCE (i.e., cell.type)  ----
  #=========================================================
  required.cols.lookup <- c("TARGET_NAME")
  
  if (show.another.banner){
    if (is.null(banner.name)) {
      stop(paste("\n****Required input is msissing: when you request to draw a new banner (other than response), you must specify which field in the lookup.table deonotes this info! \n\n"))
    } else {
      # required.cols.lookup <- c(required.cols.lookup, "SAMPLE.SOURCE")
      required.cols.lookup <- c(required.cols.lookup, toupper(banner.name))
    }
  }

  if (show.response){
    required.cols.lookup <-  c(required.cols.lookup,"RESPONSE")
  }
  
  if (!is.null(split.cols.by)){
     required.cols.lookup <-  unique(c(required.cols.lookup, toupper(split.cols.by)))
  }
 
  if (show.individuals){
    required.cols.lookup <-  c(required.cols.lookup,"INDIVIDUAL.ID")
  }
  
  if ((show.another.banner | show.response | show.individuals ) | (!is.null(split.cols.by))) {
    
    colnames(lookup.table) <- toupper(colnames(lookup.table))
    
    if (!all(required.cols.lookup %in% colnames(lookup.table))) {
      
      stop(paste("\n****Required cols for added annotations are missing from input lookup.table (please note lookup table MUST have the exact columns)) =>\t", 
                 setdiff(required.cols.lookup, colnames(lookup.table))),"\n") 
      
    } else {
      
      cat(paste("\n**** test_required_fields (5): successfull!\tAll required columns exist in lookup.table... \n"))
    }}

  cat(paste("\n**** Finished... SUCCESSFUL test of all input data fields (test_required_fields). \n"))
  
  #=========================================================

  return(list(muts= muts, svs= svs, cnvs= cnvs, 
              lookup.table=lookup.table, 
              required.cols.lookup= unique(required.cols.lookup), 
              valid.effects= mut.VT.options))
  
}