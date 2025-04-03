prepare_fill_M <- function(long_df, SAMPLES, GENES, remove.empty.cols = TRUE){ # lookup.table DEV for remove.empty.cols
  
  ###############################################################
  # == Prepare M matrix of variants  ====
  ##############################################################
  
  gene.order <- GENES
  
  long_df <- long_df %>% dplyr::filter(GENE %in% GENES & TARGET_NAME %in% SAMPLES)
  
  event_levels <- c("unknown", "other_snvs", "missense", "splice_site_variant", "splicing", "initiator_codon_change",
                    "complex", "complex_karyotype", "biallelic", "multi_hit", "stop_gain", "truncating",
                    "inframe_indel", "inframe", "frameshift_indel", "frameshift", "amp", "gain", "del", "loss",
                    "loh", "cnloh", "inv", "rearr", "fusion", "fus", "trans", "tra", "tdup", "dup", "add", "der",
                    "other_svs", "other_cnvs", "other", "unavailable", "normal", "karyotypic_abnormal")
  
  # Collapse EFFECTs by gene and sample with proper ordering
  
  long_df$EFFECT <- factor(long_df$EFFECT,  levels = event_levels)

  long_df <- long_df %>% arrange(EFFECT) %>% ungroup()
  
  grouped_df <- long_df %>%
    group_by(TARGET_NAME, GENE) %>%
    dplyr::summarise(EFFECT = paste(as.character(EFFECT), collapse = ";"), .groups = "drop") %>% ungroup()
  
  # Pivot to wide format
  mat <- grouped_df %>%
    tidyr::pivot_wider(names_from = TARGET_NAME, values_from = EFFECT, values_fill = "") %>% ungroup() %>% as.data.frame()
  
  mat[mat==0] = ""
  
  rownames(mat) <- mat$GENE
  
  mat$GENE = NULL
  
  events <- factor(unique(long_df$EFFECT), levels=event_levels)
  events <- as.character(events[order(events)])

  
  # M <- as.data.frame(matrix(0, nrow = length(GENES$genes), ncol = M.num.cols))
  # 
  # row.names(M) <- GENES$genes
  # 
  # colnames(M) <- M.col.names
  # 
  # gene.order <- gene.list$GENE
  
  # M[M==0] = ""
  
  
  
  ###############################################################
  # == Add the Event.Type in the Matrix ====
  ##############################################################
  
  cat(paste0("\nGenerating the matrix of mutations (M)...\n"))
  
  # # M[M==0] = ""
  # 
  # events <- factor(unique(data$EFFECT), levels=c("unknown","other_snvs","missense","splice_site_variant","splicing","initiator_codon_change",
  #                                                "complex","complex_karyotype","biallelic","multi_hit","stop_gain","truncating","inframe_indel", "inframe",
  #                                                "frameshift_indel","frameshift","amp","gain","del","loss","loh","cnloh","inv","rearr","fusion","trans","tra","tdup","dup","add","der",
  #                                                "other_svs","other_cnvs","other","unavailable","normal","karyotypic_abnormal"))
  # events <- events[order(events)]
  # 
  # events <- as.character(events)
  # 
  # data$GENE <- as.character(data$GENE)
  # 
  # for (i in 1: length(events)){
  #   
  #   temp <- subset(data, EFFECT==events[i])
  #   
  #   for (j in 1:nrow(temp)) {
  #     
  #     M[temp$GENE[j], temp$TARGET_NAME[j]] <- paste0(M[temp$GENE[j], temp$TARGET_NAME[j]], unique(temp$EFFECT),";", collapse = "")
  #   }
  # }
 
  return(list(M= mat,
              gene.order= gene.order,
              events= events)) 
}
