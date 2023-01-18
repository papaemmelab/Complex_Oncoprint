add_new_banner <- function(banner.name, lookup.table, list.my.cols, show.annot.legend) {
  
  cat(paste0("\nPrepare new annotation bar...\n"))
  
  # uniq.grps <- unique(lookup.table[,toupper(banner.name[1])])
  
  uniq.grps <- unique(as.vector(as.matrix(lookup.table[,toupper(banner.name)])))
  
  if (is.factor(uniq.grps)) {
    uniq.grps = levels(uniq.grps)
  }
  
  n1 <- length(uniq.grps)
  
  if (n1 <=6){
    # new.banner.col <- list$nice.cols.A[1:n1]
    # new.banner.col <- c("#C2E0D5","#ECCB80","#BE6698","#D18AB2")[1:n1]
    new.banner.col <- c("#E15759","#76B7B2","#F28E2B","#FF9DA7","#D18AB2","#808080")[1:n1]
  } else {
    new.banner.col <-  distinctColorPalette(n1)
  } 
  
  names(new.banner.col) <- uniq.grps
  list.my.cols$new.banner.col <- new.banner.col
  
  show.annot.legend <- c(show.annot.legend, "TRUE")
  
  rm(n1)
  
  num.features = length(banner.name)
  
  #----------------------------------------------------
  # Assign universal colors to most common features
  #----------------------------------------------------
  
  for (m in c("N/A","NA")) new.banner.col[[m]] <- "#444444" 
  
  for (m in c("MRD -"," 0: CR, MRD Negative","CR/CRi/MLFS MRD- (Conv.)","MRD-","CR/CRi/MLFS MRD-")) new.banner.col[[m]] <- "#80B1D3"
  
  for (m in c("CR/CRi/MLFS MRD- (after R1)")) new.banner.col[[m]] <-"#2166ac" 
  
  for (m in c("MRD+","CR/CRi/MLFS MRD+")) new.banner.col[[m]] <-"#F39C12"
  
  for (m in c("CR","C/R","Complete Response","Complete response")) new.banner.col[[m]] <- "#80B1D3"
  
  
  for (m in c(" 1: Complete Remission")) new.banner.col[[m]] <- "#80B1D3"
  
  for (m in c(" 2: CR w/Incomplete count recovery (CRi)")) new.banner.col[[m]] <- "#8c96c6"
  
  for (m in c(" 3: MFLS (Morphologic Leukemia-Free State)")) new.banner.col[[m]] <- "#d4b9da"
  
  for (m in c("PR","P/R","Partial Response"," 4: Partial Remission")) new.banner.col[[m]] <-"#F39C12"
  
  for (m in c(" 5: Primary Refractory")) new.banner.col[[m]] <- "#f768a1"
  
  for (m in c(" 5: Progressive Disease","Persistent Disease"," 4: Persistent Disease")) new.banner.col[[m]] <-"#DA2310"
  
  for (m in c("10: Hematologic Relapse (After CR)")) new.banner.col[[m]] <- "#800026"
  
  for (m in c("NR","N/R","No Response","No response")) new.banner.col[[m]] <- "#DA2310" 
  
  
  
  for (m in c("no","No","NO")) new.banner.col[[m]] <- "#C1DAD6"
  for (m in c("yes","Yes","YES")) new.banner.col[[m]] <- "#FC968B"
  
  #----------------------------------------------------
  # Assign universal colors to most common platforms
  #----------------------------------------------------
  
  for (m in c("ThunderBolts MSK-Heme")) new.banner.col[[m]] <- "#444444" 
  for (m in c("Raindance")) new.banner.col[[m]] <- "#80B1D3"
  for (m in c("Raindance-and-RAINDANCE THUNDERSTORM")) new.banner.col[[m]] <-"#2166ac"
  for (m in c("IMPACT for Hematology-and-ThunderBolts MSK-Heme")) new.banner.col[[m]] <-"#F39C12"
  for (m in c("IMPACT for Hematology")) new.banner.col[[m]] <- "#DA2310"
  # for (m in c("MRD -","MRD-")) new.banner.col[[m]] <- "#df65b0"
  # for (m in c("no","No","NO")) new.banner.col[[m]] <- "#C1DAD6"
  # for (m in c("yes","Yes","YES")) new.banner.col[[m]] <- "#FC968B"
  
  #----------------------------------------------------
  
  TEMP <- replicate(num.features, new.banner.col, simplify = FALSE)
  
  names(TEMP) = toupper(banner.name)
  
  list.my.cols = c( list.my.cols, TEMP)
  
  return(list(list.my.cols= list.my.cols,
              new.banner.col= new.banner.col,
              show.annot.legend= show.annot.legend
              
              ))
  
}