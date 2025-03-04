heatmap_colors <-  function() {
  
  library(prettyGraphs)
  library(colorspace)
  library(unikn)
  library(MetBrewer)
  
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
  # Select these from MET ----
  #==================================================================================
  
  MN <- met.brewer("Monet", type = "discrete")
  RD <- met.brewer("Redon", type = "discrete")
  FH <- met.brewer("Isfahan1", type = "discrete") 
  DM <- met.brewer("Demuth", type = "discrete") 
  HK1 <- met.brewer("Hokusai1", type = "discrete") 
  HK3 <- met.brewer("Hokusai3", type = "discrete") 
  DeR <- met.brewer("Derain", type = "discrete")
  TP <- met.brewer("Tiepolo", type = "discrete")
  LK <- met.brewer("Lakota", type = "discrete")
  #############################################
  # == Define Default Mutation colors ----
  ##############################################
  
  mut.colors.def <- list("missense"= "#3182bd",
                         "stop_gain"= "#000000", 
                         "frameshift_indel"= "#de2d26", 
                         "initiator_codon_change"= "#8460db", 
                         "splice_site_variant"= "#fdb863", ###"#598238",
                         "inframe_indel"= "#9ecae1",
                         "unknown"= "#bdbdbd",
                         "other_snvs"= "#ffeda0",
                         "complex"= "#FF7400FF")
  
  mut.colors.def <- unlist(mut.colors.def)
  
  mut.colors.def <- add.alpha(mut.colors.def, alpha = .9)
  
  #############################################
  # == Define Mutation colors ----
  ##############################################
  mut.colors <-  list( "biallelic"= DM[1],
                       "Biallelic"= DM[1],
                       
                       "complex_change_in_transcript"= RD[5],
                       "complex"= RD[5],
                       
                       # "complex_change_in_transcript"= MN[3],
                       # "complex"= MN[3],

                       "missense" = MN[3],
                       
                       "initiator_codon_change"= RD[4],
                       "stop_retained_variant"= RD[2],
                       "stop_lost"= RD[7],
                       "extended_intronic_splice_region_variant"= RD[9],
                       
                       
                       "stop_gained"= MN[7],
                       "stop_gain"= MN[7],
                       "stop_retained_variant" =  RD[10],
                       
                       "splice_site_variant"= RD[3],
                       "splicing"=RD[3],
                       
                       "frameshift"= DM[2],
                       "truncating" = DM[1],
                       
                       "inframe_codon_gain"= MN[4],
                       "inframe_codon_loss"= FH[5],
                       "inframe_variant"= MN[9],
                       "inframe"= MN[4],
                       
                       "multi_hit" = RD[9],
                       # "multi_hit" ="black",
                       
                       "multi_muts"= RD[10],
                       "multi_muts"= RD[10],
                       
                       "loss" = RD[12],
                       "Loss" = RD[12],
                       "LOSS" = RD[12],
                       
                       
                       "gain" = RD[1],
                       "Gain" = RD[1],
                       "GAIN" = RD[1],
                       
                       "loh" = RD[11],
                       "cnLOH"=  RD[11],
                       "cnloh"=  RD[11],
                       
                       "trans"= RD[3],
                       "TRA"= RD[3],
                       
                       "INV"= FH[5],
                       "inv"= FH[5],
                       
                       
                       "fusion"= RD[4],
                       "FUS"= RD[4],
                       
                       "SNV" = RD[4],
                       "INDEL"= RD[2]
  )
  
  # mut.colors <- list("missense"= "#3182bd",
  #                    "non_synonymous_codon"= "#3182bd",
  #                    
  #                    "stop_gain"= "#000000", 
  #                    "stop_gained" = "#000000", 
  #                    
  #                    # "biallelic"= "#000000", 
  # 
  #                    "frameshift_indel"= "#de2d26", 
  #                    "frameshift_variant"= "#de2d26", 
  #                    
  #                    "initiator_codon_change"= "#8460db", 
  #                    
  #                    "splice_site_variant"= "#fdb863", ###"#598238",
  #                    
  #                    
  #                    "inframe_indel"= "#9ecae1",
  #                    "inframe_variant" = "#9ecae1",
  #                    "inframe_codon_gain"= HK3[6],
  #                    "inframe_codon_loss" = MN[1],
  #                    
  #                    "unknown"= "#bdbdbd",
  #                    "other_snvs"= "#ffeda0",
  #                    "complex"= MN[4],
  #                    "complex_change_in_transcript"= MN[4],
  #                    "stop_retained_variant" =  HK3[7],
  #                    
  #                    "extended_intronic_splice_region_variant"= RD[9],
  #                    "biallelic"= RD[8],
  #                    
  #                    "complex_change_in_transcript"= HK1[4],
  #                    "complex"= RD[7],
  #                    
  #                    
  #                    
  #                    "Missense" = MN[3],
  #                    
  #                    "initiator_codon_change"= RD[4],
  #                    "stop_retained_variant"= RD[2],
  #                    "stop_lost"= RD[8],
  #                    "extended_intronic_splice_region_variant"= RD[9],
  #                    
  #                    
  #                    "stop_gained"= MN[7],
  #                    "stop_retained_variant" =  RD[10],
  #                    
  #                    "splice_site_variant"= RD[3],
  #                    "splicing"=RD[3],
  #                    
  #                    "frameshift"= DM[2],
  #                    
  #                    "inframe_codon_gain"= MN[4],
  #                    "inframe_codon_loss"= FH[5],
  #                    # "inframe_variant"= MN[9],
  #                    "inframe"= MN[4],
  #                    
  #                    "truncating" = DM[1],
  #                    
  #                    "SNV" = RD[1],
  #                    "INDEL"= RD[2],
  #                    
  #                    "Multi_hit_MUT" = HK1[3],
  #                    "multi_hit" = RD[9],
  #                    
  #                    "cnloh" = RD[11],
  #                    
  #                    "TRA"= RD[3],
  #                    "INV"= RD[4])
  
  mut.colors <- unlist(mut.colors)
  
  # Define specific names to exclude from alpha modification
  specific_names <- c("biallelic", "Biallelic")
  
  # Apply alpha only to colors not in `specific_names`
  mut.colors[!names(mut.colors) %in% specific_names] <- add.alpha(mut.colors[!names(mut.colors) %in% specific_names], 0.6)
  
  # mut.colors <- add.alpha(mut.colors, alpha = 0.6)
  
  #############################################
  # == Define CNV/Cytogenetics colors ----
  #############################################
  
  cyto.colors <- list(
    # "amp"= "#3288bd", # AMP
    # "del"= '#d53e4f', # DEL
    # "inv"= '#feb24c', # INV
    # "normal"= '#c7eae5'  # Normal_Karyotype
    
    "amp"= "#3C33FF",        # AMP
    "del"= RD[8],        # DEL #"#FB3F28",
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
    "normal"= '#a8ddb5',     # Normal_Karyotype
    # below events are added by Jesus' reco
    'CNLOH' = '#4d9221',
    'GAIN' = MN[8],     # gain for CN-state=3
    # 'LOSS' = RD[8],
    'LOSS' = "#313695",
    'deep_LOSS' = RD[8],
    "cnloh" = RD[11],
    
    "tra"= RD[3],
    "fus"= RD[3],
    "inv"= RD[4],
    
    "biallelic"= RD[8]
    
  )
  
  cyto.colors <- unlist(cyto.colors)
  
  cyto.colors <- add.alpha(cyto.colors, alpha = .7)
  
  ###########################################
  # == Complex Karyotype colors (Publication case) ----
  ###########################################
  
  complex.colors <- list("complex"= RD[10],"not complex"= RD[11],"not available"="grey")
  complex.colors <- unlist(complex.colors)
  
  ###########################################
  # == ELN colors (Publication case) ----
  ###########################################
  
  eln.molecular.response.colors <- list("CR"= "#80b1d3","PR"= "#b2df8a","NR"="#fb8072","N/E"="#d9d9d9")
  response.colors <- unlist(eln.molecular.response.colors)
  
  ####################################################
  # == pathology.report colors (special case) ----
  ####################################################
  
  path.colors <- c("#80b1d3","#fb8072","#d9d9d9")
  names(path.colors) <- c("HR","NR","N/E")
  
  ###############################
  # == Response colors ----
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
  # == therapy colors ----
  ###############################  
  
  therapy.colors <- list("No resistance/intolerance"= "#80b1d3", 
                         "No.Resistance"= "#80b1d3", 
                         "Intolerant"= "#fee090", 
                         "Resistant"= "#fb8072", 
                         "MRD Negative"= "#df65b0", 
                         "Inevaluable"= "#d9d9d9", 
                         "#N/E"="#d9d9d9") #N/E
  
  therapy.colors <- unlist(therapy.colors)
  
  ########################################
  # == Nice.cols for Publication ----
  ########################################
  
  # this was initially used for "J Grinfeld et al. Classification" ribbon in 157 manuscript and I liked the combo (TP53 mutation etc); 
  # the names must be adjusted based on the ribbon features
  
  nice.cols.A <- c("#C1DAD6", #TP53 mutation"
                   "#6D929B", #"Chromatin/Spliceosome/RAS mutation
                   "#CCFFCC", #"CALR mutation"
                   "#FFCF79", #"MPL mutation"
                   "#B7AFA3", #"homozygous JAK2/NFE2 mutation"
                   "#E8D0A9", #"heterozygous JAK2"
                   "#CCCCCC", #"Other drivers"
                   "#666666"#"No drivers"
  )
  
  ######################################
  # == Some other coloring options ---
  ######################################
  
  hcl.pals(type = "divergingx")
  A <- hcl.colors(9, palette = "Zissou1", alpha = 0.7, rev = FALSE, fixup = TRUE)
  #seecol(A) #### <<< displays colors
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
  
  ###############################
  # == ALL cols
  ############################### 
  
  # subtype_colors <- list(
  #   "No-RNA" = "#e5e4e2",
  #   "No-DNA" = "#a9a9a9",
  #   "Excluded" = "#C2C2C2",
  #   "Unclassified" = "#000000",
  #   "KMT2A Group" = "#b25333", 
  #   "KMT2A" = "#b25333",
  #   "Gene fusions" = "#C2C2C2",
  #   "Ph" = "#762a83",
  #   "BCR-ABL1" = "#762a83",
  #   "Hypodiploid" = "#8073ac",
  #   "Low hypodiploid" = "#8073ac",
  #   "CEBP" = "#e08214",
  #   "Near haploid" = "#dfc27d",
  #   "PAX5alt" = "#5aae61",
  #   "TCF3-PBX1" = "#92c5de",
  #   "DUX4" = "#313695",
  #   "High hyperdiploid" = "#d6604d",
  #   "Hyperdiploid" = "#d6604d",
  #   "Ploidy subtype" = "#C2C2C2",
  #   "Ph-like" = "#c2a5cf",
  #   "BCR-ABL1-like" = "#c2a5cf",
  #   "ZNF384" = "#35978f",
  #   "ZNF384 Group" = "#35978f",
  #   "ZNF384-like" = "#542788",
  #   "BCL2/MYC" = "#01665e",
  #   "MEF2D" = "#8c510a",
  #   "iAMP21" = "#fdb863",
  #   "Other" = "#C2C2C2",
  #   "CDX2/UBTF" = "#bf812d",
  #   "CDX2_UBTF" = "#bf812d",
  #   "PAX5 P80R" = "#a6dba0"
  # )
  
  ###############################
  # NEW Add Met colors -------
  ###############################
  
  subtype_colors <- list(
    "No-RNA" = "#FFFFFF",
    "No-DNA" = "#FFFFFF",
    "Excluded" = "#000000",
    "excluded" = "#000000",
    "Not Available" = "#FFFFFF",
    "No-WGS"=  "#FFFFFF",
    "No-RNA"=  "#FFFFFF",
    "Normal"= MN[9],
    "Other" = "#C2C2C2",
    "other"= "#C2C2C2",
    "Singleton" = DM[10],
    "Gene fusions" = "#C2C2C2",
    
    #===========================
    
    "BCL2/MYC" = RD[9], #DeR[6],
    "BCR-ABL1" = MN[2],
    "BCR-ABL1-like" = MN[3],
    
    "CDX2/UBTF" = RD[12],
    "CDX2_UBTF" = RD[12],
    "CEBP" = RD[1],
    
    "DUX4" = RD[2],
    "ETV6-RUNX1"= HK3[1],
    
    "iAMP21" = "#fdb863",
    "IKZF1 N159Y"= FH[5],
    
    "KMT2A Group" = HK1[3], #nice.cols.A[6],
    
    "Hypodiploid" = MN[4],
    "Low hypodiploid" = MN[4], #RD[8],
    
    "High hyperdiploid" = HK1[1], #RD[8], #MN[4],
    "Hyperdiploid" = HK1[1], #RD[8], #MN[4],
    
    "MEF2D" = FH[1],
    "Near haploid" = MN[9],
    
    "PAX5 P80R" = RD[11],
    "PAX5alt" = RD[7],
    
    "Ph" = MN[2],
    "Ph-like" = MN[3],
    
    "TCF3-PBX1" = "#92c5de",
    
    "ZNF384" = RD[9],
    "ZNF384 Group" = DeR[6], # RD[9],
    "ZNF384-like" = RD[10],
    
    "Unclassified" = DM[9],
    "Low Quality" = DM[7],
    "failed"= DM[10],
    "Failed"= DM[10],
    "FAILED"= DM[10],
    "NoData"= "#FFFFFF", # white
    "ABL-class"= HK3[1],
    "CRLF2" = HK3[2],
    "P2RY8-CRLF2"= HK3[3],
    "Complex karyotype" = DM[4]
    
    #===========================
    
  )
  
  # Convert the list to a named vector if needed
  ALL.SUBTYPE <- unlist(subtype_colors)
  ALL.SUBTYPE <- ALL.SUBTYPE[sort(names(ALL.SUBTYPE))]
  
  ALL.SUBTYPE <- unlist(ALL.SUBTYPE)
  
  ALL.SUBTYPE <- add.alpha(ALL.SUBTYPE, alpha = .9)
  
  ###############################
  # == GENDER cols ----
  ############################### 
  
  GENDER <- list("MALE"= K[5], 
                 "FEMALE"= K[1], 
                 "UNKNOWN"= K[3]  
  )
  GENDER <- unlist(GENDER)
  
  ###############################
  ALL.mut.cols <- list("non_synonymous_codon"= MN[3],
                       "frameshift_variant" = MN[1],
                       "stop_gained"= MN[7],
                       
                       "complex_change_in_transcript"= HK1[4],
                       "stop_retained_variant" =  RD[10],
                       
                       "inframe_variant"= MN[9],
                       "splice_site_variant"= RD[3],
                       
                       "inframe_codon_gain"= MN[4],
                       "inframe_codon_loss"= FH[5],
                       
                       "initiator_codon_change"= RD[4],
                       "stop_retained_variant"= RD[2],
                       "stop_lost"= RD[8],
                       "extended_intronic_splice_region_variant"= RD[9]
  )
  
  mut.colors.def <- add.alpha(mut.colors.def, alpha = .6)
  
  ###############################
  # == Purity color
  ###############################
  
  purity_colors <- c(
    "<20" = "#3E4A89",      # Dark blue for <20
    "20-40" = "#586BA4",    # Medium-dark blue for 20-40
    "40-60" = "#7D93B2",    # Medium blue for 40-60
    "60-80" = "#A9BEDB",    # Light blue for 60-80
    ">80" = "#D4E4F7",    # Very light blue for 80-100
    "No-WGS"= "white",
    "No-RNA"= "white"
  )
  
  #====================================================
  
  return(list(mut.colors= mut.colors, 
              default.mut.colors= mut.colors.def,
              #cyto.colors=cyto.colors,  
              response.colors= response.colors,
              nice.cols.A= nice.cols.A,
              therapy.colors=therapy.colors,
              eln.molecular.response.colors,
              path.colors= path.colors,
              ALL.SUBTYPE= ALL.SUBTYPE,
              ALL.mut.colors= ALL.mut.cols,
              GENDER= GENDER,
              PURITY= purity_colors,
              EE= purity_colors,
              COMPLEX.KARYOTYPE= complex.colors
  )) 
  
}