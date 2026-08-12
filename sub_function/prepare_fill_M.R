prepare_fill_M <- function(long_df, SAMPLES, GENES, remove.empty.cols = TRUE, show.multis= FALSE){ # lookup.table DEV for remove.empty.cols
  
  ###############################################################
  # == Prepare M matrix of variants  ====
  ##############################################################
  # browser()
  
  gene.order <- GENES
  
  long_df <- long_df %>% dplyr::filter(GENE %in% GENES & TARGET_NAME %in% SAMPLES)
  
  event_levels <- c("unknown", "other_snvs", "missense", "splice_site_variant", "splicing", "initiator_codon_change",
                    "complex", "complex_karyotype", "biallelic", "multi_hit", "stop_gain", "truncating",
                    "inframe_indel", "inframe", "frameshift_indel", "frameshift", "amp", "gain", "del", "loss",
                    "loh", "cnloh", "inv", "rearr", "fusion", "fus", "trans", "tra", "tdup", "dup", "add", "der",
                    "other_svs", "other_cnvs", "other", "unavailable", "normal", "karyotypic_abnormal","iso")
  
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
  
  # pivot_wider only creates a column for samples that had >=1 matching mutation, so
  # samples with 0 selected mutations are silently missing here (remove.empty.cols was
  # never actually consulted) - add them back as all-empty columns, then honor the flag
  missing.samples <- setdiff(SAMPLES, colnames(mat))
  
  if (length(missing.samples) > 0) {
    mat[missing.samples] <- ""
  }
  
  mat <- mat[, intersect(SAMPLES, colnames(mat)), drop = FALSE]
  
  if (remove.empty.cols) {
    mat <- mat[, colSums(mat != "") > 0, drop = FALSE]
  }
  
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
  
  if (show.multis){
    
    cat(paste0("\nStart multi.hit Oncoprint preparation...\n"))
    
    multi.hits <- long_df %>% dplyr::group_by(TARGET_NAME, GENE) %>% dplyr::mutate(N= n()) %>% dplyr::filter(N>1) %>% dplyr::select(TARGET_NAME, GENE) %>% unique()
    
    multi.hits <- data.frame(multi.hits)
    
    if (nrow(multi.hits) > 0) {
      for (k in seq_len(nrow(multi.hits))) {
        gene <- as.character(multi.hits$GENE[k])
        sample <- as.character(multi.hits$TARGET_NAME[k])
        val <- mat[gene, sample]
        
        if (!grepl("biallelic", val)) {
          mat[gene, sample] <- if (nchar(val) > 0) paste0(val, ";multi_hit") else "multi_hit"
        }
      }
    }
  }
 
  mat = mat[gene.order, ]
  
  return(list(M= as.matrix(mat),
              gene.order= gene.order,
              events= events)) 
}
