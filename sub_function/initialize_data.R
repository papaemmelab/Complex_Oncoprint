initialize_data <- function(data, muts, cnvs, svs, muts.order, cnvs.order, svs.order, min.freq, lookup.table, sec.1.label= sec.1.label, sec.2.label=NULL, sec.3.label=NULL, REQ.cols, save.path, params){
  
  ###############################################
  save.name = eval(params$save.name)
  

  if (is.null(save.name)){
    saveFile.1 <- file.path(save.path,"TEMP",paste0("Heatmap_TEMP_minFreq_",min.freq,".png"))
    saveFile.2 <- file.path(save.path,paste0("Oncoprint_minFreq_",min.freq,".png"))
    
  } else {
    saveFile.1 <- file.path(save.path,"TEMP",paste0("Heatmap_TEMP_minFreq_",save.name,"_",min.freq,".png"))
    saveFile.2 <- file.path(save.path,paste0("Oncoprint_minFreq_",save.name,"_",min.freq,".png"))
    
  }
  
  ###############################################
  # Create FONT.obj for calling simple.ht  ====
  ###############################################
  font.obj <- list(fig.title.font= eval(substitute(fig.title.font), parent.frame()),
                   legend.title.font= eval(substitute(legend.title.font), parent.frame()),
                   cols.font= eval(substitute(cols.font), parent.frame()),
                   pct.font= eval(substitute(pct.font), parent.frame()),
                   rows.font= eval(substitute(rows.font), parent.frame()),
                   barplot.font= eval(substitute(barplot.font), parent.frame()),
                   legend.label.font= eval(substitute(legend.label.font), parent.frame())
  )
  

  ###############################################################
  # == Add MUTATIONS  ====
  ##############################################################
  
  source(file.path("./sub_function/sort_variants.R"))
  
  A <- sort_variants(muts, muts.order,group.label= sec.1.label, variants.class= "MUTs")
  
  data <- A[[1]][,c("TARGET_NAME","GENE","EFFECT")]
  
  gene.list <- A$gene.list
  
  ###############################################################
  # == Add CNVs  ====
  ##############################################################
  
  if (!is.null(cnvs)){
    
    if (nrow(cnvs)==0){
      stop("\n Size of cnvs = 0. This may be due to the fact that you currently do not have any CNV with the min.freq your specified.\nby default min.freq will be applied to mutations as well as CNVs/SVs. Try re-defining svs-specific 'min.freq.svs' param...\n")
    }
    
    source(file.path("./sub_function/sort_variants.R"))
    
    D <- sort_variants(cnvs, cnvs.order, group.label= sec.2.label , variants.class= "CNVs") 
    
    data <- rbind(data, D[[1]][,c("TARGET_NAME","GENE","EFFECT")])
    
    gene.list <- rbind(gene.list, D$gene.list)
  }
  
  ###############################################################
  # == Add Structural Variants  ====
  ##############################################################
  
  if (!is.null(svs)){
    
    if (nrow(svs)==0){
      stop("\n Size of svs = 0. This may be due to the fact that you currently do not have any svs with the min.freq your specified. Try re-defining the 'min.freq.svs' param...\n")
    }
    
    source(file.path("./sub_function/sort_variants.R"))
    
    B <- sort_variants(svs, svs.order, group.label= sec.3.label , variants.class= "SVs") 
    
    data <- rbind(data, B[[1]][,c("TARGET_NAME","GENE","EFFECT")])
    
    gene.list <- rbind(gene.list, B$gene.list)
    
  }
  
  ###############################################################

  
  return(list(saveFile.1= saveFile.1,
              saveFile.2= saveFile.2,
              font.obj= font.obj,
              
              data= data,
              muts= muts,
              
              gene.list= gene.list

              ))
}
