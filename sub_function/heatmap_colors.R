heatmap_colors <-  function() {
  
  library(prettyGraphs)
  
  #==================================================================================
  #   Written by Noushin Farnoud, Jul 2018, Latest Update Jan 2024
  #----------------------------------------------------------------------------------
  #   The aim of this function is to use a default color for mutations and indels in 
  #   oncoprint. 
  #   You can change this to accomodate further VARIANT TYPEs, but for now refer to 
  #   function generate_oncoprint.R to see the impact of variant types in your data.
  # 
  #   Contact Noushin Farnoud (rahnaman@mskcc.org) if you faced any error.
  # 
  #    See also example_Heatmap, test_required_fields, generate_oncoprint.
  #==================================================================================

  mut.colors <- list("missense"= "#3182bd",
                     "stop_gain"= "#000000", 
                     "frameshift_indel"= "#de2d26", 
                     "initiator_codon_change"= "#8460db", 
                     "splice_site_variant"= "#fdb863", ###"#598238",
                     "inframe_indel"= "#9ecae1",
                     "unknown"= "#bdbdbd",
                     "other_snvs"= "#ffeda0",
                     "complex"= "#FF7400FF")
  
  mut.colors <- unlist(mut.colors)
  
  mut.colors <- add.alpha(mut.colors, alpha = .9)
  
  ###############################
  # == Define CNV/Cytogenetics colors
  ###############################  
  
  cyto.colors <- list(
    # "amp"= "#3288bd", # AMP
    # "del"= '#d53e4f', # DEL
    # "inv"= '#feb24c', # INV
    # "normal"= '#c7eae5'  # Normal_Karyotype
    
    "amp"= "#3C33FF",        # AMP
    "del"= "#FB3F28",        # DEL
    "loh"= "#fdae61",        # LOH
    "complex_karyotype"= "#6a3d9a", # complex_karyotype
    "inv"= "#c51b8a",       # INV
    "complex"= '#FF7400FF', # complex
    "trans"= "#66c2a5",     # TRANS
    "rearr"= "darkgreen",   # REARRANGEMENT
    "tdup"= "#cab2d6",      # TDUPs
    "dup"= "#31a354",       # DUP
    "other_svs"= "#c51b8a", # OTHER_SVs
    "other_cnvs"= "#fa9fb5",# OTHER_CNVs
    "fusion"= "#fdbf6f",    # FUSION
    "N/A"= "grey",           # N/E or inconclusive
    "add"= '#fdae61',       # ADD
    "der"= '#fee08b',       # DER
    "normal"= '#a8ddb5'     # Normal_Karyotype
  )
  
  cyto.colors <- unlist(cyto.colors)
  
  cyto.colors <- add.alpha(cyto.colors, alpha = .7)

  ###############################
  # == ELN colors (special case)
  ###############################  
  
  eln.molecular.response.colors <- list("CR"= "#80b1d3","PR"= "#b2df8a","NR"="#fb8072","N/E"="#d9d9d9")
  response.colors <- unlist(eln.molecular.response.colors)
  
  #############################################
  # == pathology.report colors (special case)
  #############################################
  
  path.colors <- c("#80b1d3","#fb8072","#d9d9d9")
  names(path.colors) <- c("HR","NR","N/E")
  
  ###############################
  # == Response colors
  ###############################  
  
  response.colors <- list("persistent"="#016c59", 
                       "partial response"= "#1c9099", 
                       "non-responder"= "#67a9cf", 
                       "stable disease"= "#bdc9e1", 
                       "responder"= "#df65b0", 
                       "CR"= "#80B1D3", 
                       "CR-i"= "#80B1D3",
                       "PR"= "#F39C12", 
                       "NR"= "#DA2310", 
                       "N/A"= "#bdbdbd", 
                       "N/E"= "#bdbdbd"  
                       )
  
  response.colors <- unlist(response.colors)
  
  ###############################
  # == therapy colors
  ###############################  
  
  therapy.colors <- list("No resistance/intolerance"= "#80b1d3", 
                      "No.Resistance"= "#80b1d3", 
                      "Intolerant"= "#fee090", 
                      "Resistant"= "#fb8072", 
                      "MRD Negative"= "#df65b0", 
                      "Inevaluable"= "#d9d9d9", 
                      "#N/E"="#d9d9d9") #N/E
  
  therapy.colors <- unlist(therapy.colors)
  
  ###############################
  # == Nice.cols
  ############################### 
  # this was initially used for "J Grinfeld et al. Classification" ribbon in 157 manuscript and I liked the combo (TP53 mutation etc); the names must be adjusted based on the ribbon features
  
  nice.cols.A <- c("#C1DAD6", #TP53 mutation"
                  "#6D929B", #"Chromatin/Spliceosome/RAS mutation
                  "#CCFFCC", #"CALR mutation"
                  "#FFCF79", #"MPL mutation"
                  "#B7AFA3", #"homozygous JAK2/NFE2 mutation"
                  "#E8D0A9", #"heterozygous JAK2"
                  "#CCCCCC", #"Other drivers"
                  "#666666"#"No drivers"
  )
  
  ###############################
  # == ALL cols
  ############################### 
  
  library(colorspace)
  library(unikn)
  
  # hcl_palettes(palette= "Zissou1", plot = TRUE)
  
  
  hcl.pals(type = "divergingx")
  A <- hcl.colors(9, palette = "Zissou1", alpha = 0.7, rev = FALSE, fixup = TRUE)
  seecol(A)
  
  hcl.pals(type = "diverging")
  B <- hcl.colors(9, palette = "Reds", alpha = 1, rev = FALSE, fixup = TRUE)
  # seecol(B)
  
  C <- hcl.colors(9, palette = "Grays", alpha = 0.8, rev = FALSE, fixup = TRUE)
  # seecol(C)
  
  D <- hcl.colors(9, palette = "Roma", alpha = 0.6, rev = FALSE, fixup = TRUE)
  
  
  G <- hcl.colors(9, palette = "Reds", alpha = 1, rev = FALSE, fixup = TRUE)
  # seecol(G)
  
  K <- hcl.colors(5, palette = "Blue-Yellow 3", rev = FALSE, fixup = TRUE)
  # seecol(K)
  
  # Noushin's prev cols
  
  # ALL.SUBTYPE <- c("#636363", # No-RNA
  #                  C[5], # "No-DNA"
  #                      C[2],  # excluded
  #                      "dark grey",  #Unclassified
  # 
  #                      D[7],  # KMT2A Group
  #                      D[7],  # KMT2A
  #                  
  #                      D[7],  # Gene fusions  # for simplified graph
  #                      D[8],  # Ph
  #                      D[8],  # BCR-ABL1
  # 
  #                      A[6],  # Hypodiploid
  #                      A[6],  # Low hypodiploid
  # 
  #                      A[1],  # CEBP
  #                      A[4],  # Near haploid
  #                      G[5],  # PAX5alt
  #                      K[2],  # TCF3-PBX1
  # 
  #                      D[1],   # DUX4
  #                      D[9],   # High hyperdiploid   ####
  #                      D[9],   # Hyperdiploid
  #                      D[9],   # Ploidy Subtype
  # 
  # 
  #                      D[6],  # Ph-like  ####
  #                      D[6],  # BCR-ABL1-like  ####
  # 
  #                      G[2],  # ZNF384  ####
  #                      G[2],  # ZNF384 Group ####
  #                   "yellow", # "ZNF384-like"
  #                      G[6],  # BCL2/MYC  ####
  #                      A[2],  # MEF2D  ####
  #                      B[1],  # iAMP21  ####
  #                      D[4],  # Other  ####
  #                      D[2],  #"CDX2/UBTF
  #                      D[2],  # CDX2_UBTF
  #                      G[7],  # PAX5 P80R"
  #                      K[5]
  # )
  
  subtype_colors <- list(
    "No-RNA" = "#e5e4e2",
    "No-DNA" = "#a9a9a9",
    "Excluded" = "#C2C2C2",
    "Unclassified" = "#000000",
    "KMT2A Group" = "#b25333",
    "KMT2A" = "#b25333",
    "Gene fusions" = "#C2C2C2",
    "Ph" = "#762a83",
    "BCR-ABL1" = "#762a83",
    "Hypodiploid" = "#8073ac",
    "Low hypodiploid" = "#8073ac",
    "CEBP" = "#e08214",
    "Near haploid" = "#dfc27d",
    "PAX5alt" = "#5aae61",
    "TCF3-PBX1" = "#92c5de",
    "DUX4" = "#313695",
    "High hyperdiploid" = "#d6604d",
    "Hyperdiploid" = "#d6604d",
    "Ploidy subtype" = "#C2C2C2",
    "Ph-like" = "#c2a5cf",
    "BCR-ABL1-like" = "#c2a5cf",
    "ZNF384" = "#35978f",
    "ZNF384 Group" = "#35978f",
    "ZNF384-like" = "#542788",
    "BCL2/MYC" = "#01665e",
    "MEF2D" = "#8c510a",
    "iAMP21" = "#fdb863",
    "Other" = "#C2C2C2",
    "CDX2/UBTF" = "#bf812d",
    "CDX2_UBTF" = "#bf812d",
    "PAX5 P80R" = "#a6dba0"
  )
  
  # Convert the list to a named vector if needed
  ALL.SUBTYPE <- unlist(subtype_colors)
  
  
  # names(ALL.SUBTYPE) <- c("No-RNA","No-DNA","excluded" ,  "Unclassified", 
  #                              "KMT2A Group","KMT2A", "Gene fusions" ,"Ph"     ,   "BCR-ABL1", 
  #                              "Hypodiploid"    ,  "Low hypodiploid", "CEBP"     ,    "Near haploid","PAX5alt"  ,    "TCF3-PBX1"  ,  "DUX4", 
  #                              "Hyperdiploid", "High hyperdiploid", "Ploidy subtype",
  #                              "Ph-like", "BCR-ABL1-like", "ZNF384", "ZNF384 Group","ZNF384-like","BCL2/MYC", 
  #                              "MEF2D","iAMP21","Other","CDX2/UBTF","CDX2_UBTF","PAX5 P80R"
  #                         )
  
  ###############################
  # == GENDER cols
  ############################### 
  
  GENDER <- list("MALE"= K[5], 
              "FEMALE"= K[1], 
              "UNKNOWN"= K[3]  
  )
  GENDER <- unlist(GENDER)
  
  
  ############################################################## 

  return(list(mut.colors= mut.colors, 
              cyto.colors=cyto.colors,  
              response.colors= response.colors,
              nice.cols.A= nice.cols.A,
              therapy.colors=therapy.colors,
              eln.molecular.response.colors,
              path.colors= path.colors,
              ALL.SUBTYPE= ALL.SUBTYPE,
              GENDER= GENDER
              )) 
  
}