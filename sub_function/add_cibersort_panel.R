add_cibersort_panel <- function(ht, M, added.heatmap, legend.title.font= 10, legend.label.font= 10, rows.font= 10, annot.title.side= "top", show.sample.names= FALSE){
  
  # browser()
  
  source(file.path("./sub_function/create_lineage_heatmaps.R"))
  
  # == what cols are missing from RNA deconv that are in oncoplot M
  
  missing_cols <- setdiff(colnames(M), colnames(added.heatmap$cell_state_zscores))
  
  # == Add NA to deconv for missing
  
  cell.state.mod <- added.heatmap$cell_state_zscores %>% bind_cols(data.frame(matrix(NA, nrow = nrow(added.heatmap$cell_state_zscores), ncol = length(missing_cols), 
                                                                                     dimnames = list(NULL, missing_cols))))
  rownames(cell.state.mod) <- rownames(added.heatmap$cell_state_zscores)
  
  dev.index.mod <- added.heatmap$B_dev_index %>% bind_cols(data.frame(matrix(NA, nrow = nrow(added.heatmap$B_dev_index), ncol = length(missing_cols), 
                                                                             dimnames = list(NULL, missing_cols))))
  
  ling.index.mod <- added.heatmap$lineage_index %>% bind_cols(data.frame(matrix(NA, nrow = nrow(added.heatmap$lineage_index), ncol = length(missing_cols), 
                                                                                dimnames = list(NULL, missing_cols))))
  # == Now sort according to the M sample order
  
  added.heatmap$cell.state <- cell.state.mod[,colnames(M)]
  added.heatmap$dev.index  <- dev.index.mod[,colnames(M)]
  added.heatmap$ling.index <- ling.index.mod[,colnames(M)]
  
  # == Run adding extra heatmaps
  
  cibersort.heatmaps <- create_lineage_heatmaps (heat.2.df= added.heatmap, 
                                                 # heat.2.name= c("Cell-population z-score","Scaled Index"),
                                                 legend.title.font = legend.title.font,
                                                 legend.label.font = legend.label.font,
                                                 annot.title.side = annot.title.side,
                                                 rows.fs= rows.font,
                                                 cols.fs= 10, 
                                                 show.sample.names = show.sample.names )
  
  hh <- ht %v%  cibersort.heatmaps$heat.1 %v% cibersort.heatmaps$heat.2 %v% cibersort.heatmaps$heat.3 #cibersort.heatmaps$heat.1 %v%
  
  return(hh)
}
