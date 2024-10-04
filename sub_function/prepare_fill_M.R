prepare_fill_M <- function(data, SAMPLES, GENES, lookup.table, rem.empty, gene.list ){
  
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
  # == Add the Event.Type in the Matrix ====
  ##############################################################
  
  cat(paste0("\nGenerating the matrix of mutations (M)...\n"))
  
  # M[M==0] = ""
  
  events <- factor(unique(data$EFFECT), levels=c("unknown","other_snvs","missense","splice_site_variant","splicing","initiator_codon_change",
                                                 "complex","complex_karyotype","biallelic","multi_hit","stop_gain","truncating","inframe_indel", "inframe",
                                                 "frameshift_indel","frameshift","amp","cngain","del","cnloss","loh","cnloh","inv","rearr","fusion","trans","tra","tdup","dup","add","der",
                                                 "other_svs","other_cnvs","other","unavailable","normal","karyotypic_abnormal"))
  events <- events[order(events)]
  
  events <- as.character(events)
  
  data$GENE <- as.character(data$GENE)
  
  for (i in 1: length(events)){
    
    temp <- subset(data, EFFECT==events[i])
    
    for (j in 1:nrow(temp)) {
      
      M[temp$GENE[j], temp$TARGET_NAME[j]] <- paste0(M[temp$GENE[j], temp$TARGET_NAME[j]], unique(temp$EFFECT),";", collapse = "")
    }
  }
 
  return(list(M= M,
              gene.order= gene.order,
              events= events)) 
}
