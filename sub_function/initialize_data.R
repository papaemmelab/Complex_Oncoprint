initialize_data <- function(data, muts, cnvs, svs, muts.order, cnvs.order, svs.order, lookup.table, REQ.cols, save.path, params){
  
  ###############################################
  save.name = params$save.name
  
  min.freq <- params$min.freq
  

  if (is.null(save.name)){
    saveFile.1 <- file.path(save.path,"TEMP",paste0("Heatmap_TEMP_minFreq_",min.freq,".jpg"))
    saveFile.2 <- file.path(save.path,paste0("Heatmap_minFreq_",min.freq,".jpg"))
    
  } else {
    saveFile.1 <- file.path(save.path,"TEMP",paste0("Heatmap_TEMP_minFreq_",min.freq,"_",save.name,".jpg"))
    saveFile.2 <- file.path(save.path,paste0("Heatmap_minFreq_",min.freq,"_",save.name,".jpg"))
    
  }
  
  ###############################################
  # Create FONT.obj for calling simple.ht  ====
  ###############################################
  font.obj <- list(fig.title.font= params$fig.title.font,
                   legend.title.font= params$legend.title.font,
                   cols.font= params$cols.font,
                   pct.font= params$pct.font,
                   rows.font= params$rows.font,
                   barplot.font= params$barplot.font,
                   legend.label.font= params$legend.label.font
  )
  

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
    
    gene.list <- rbind(gene.list, B$gene.list)
    
  }
  
  ###############################################################
  # == Add CNVs  ====
  ##############################################################
  
  if (!is.null(cnvs)){
    
    source(file.path("./sub_function/sort_variants.R"))
    
    D <- sort_variants(cnvs, cnvs.order, group.label= "Cytogenetics", variants.class= "CNVs")
    
    data <- rbind(data, D[[1]][,c("TARGET_NAME","GENE","EFFECT")])
    
    gene.list <- rbind(gene.list, D$gene.list)
  }
  

  
  return(list(saveFile.1= saveFile.1,
              saveFile.2= saveFile.2,
              font.obj= font.obj,
              
              data= data,
              muts= muts,
              
              gene.list= gene.list

              ))
}
