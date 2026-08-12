match_numerics_to_oncoprint_layout<- function(oncoprint.mat, lookup.table = NULL, lookup.val.col= NULL){
  
  # Extract the column name as a symbol
  val_col_sym <- rlang::ensym(lookup.val.col)
  
  # Convert matrix to long format
  mat_df <- as.data.frame(oncoprint.mat) %>%
    rownames_to_column("GENE") %>%
    pivot_longer(-GENE, names_to = "Individual.System.ID", values_to = "event")
  
  # Join survival time info
  mat_df <- mat_df %>%
    left_join(lookup.table %>% dplyr::select(Individual.System.ID, !!val_col_sym), by = "Individual.System.ID") %>%
    mutate(.lookup_val = !!val_col_sym) %>%
    mutate(value = ifelse(event != "", .lookup_val, ""))
  
  # Reshape back to wide matrix
  output_mat <- suppressWarnings(mat_df %>%
    dplyr::select(GENE, Individual.System.ID, value) %>% mutate(value=as.numeric(value)) %>%
    pivot_wider(names_from = Individual.System.ID, values_from = value) %>%
    column_to_rownames("GENE") %>%
    as.matrix())
  
  return(output_mat)
}