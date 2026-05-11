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
  # browser()
  
  MN <- met.brewer("Monet", type = "discrete")
  RD <- met.brewer("Redon", type = "discrete")
  FH <- met.brewer("Isfahan1", type = "discrete") 
  DM <- met.brewer("Demuth", type = "discrete") 
  HK1 <- met.brewer("Hokusai1", type = "discrete") 
  HK3 <- met.brewer("Hokusai3", type = "discrete") 
  DeR <- met.brewer("Derain", type = "discrete")
  TP <- met.brewer("Tiepolo", type = "discrete")
  LK <- met.brewer("Lakota", type = "discrete")
  CAS1 <- met.brewer("Cassatt1", type = "discrete")
  CAS2 <- met.brewer("Cassatt2", type = "discrete")
  
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
                       
                       "complex_change_in_transcript"= HK1[3],
                       "complex"= HK1[3],
                       
                       # "complex_change_in_transcript"= MN[3],
                       # "complex"= MN[3],

                       # "missense" = MN[2],
                       
                       "missense" = HK1[6],
                       
                       "initiator_codon_change"= RD[4],
                       
                       "stop_lost"= MN[1],
                       
                       "stop_retained_variant"= RD[2], #RD[2],
                       "stop_gained"=  RD[2], #MN[7],
                       "stop_gain"=  RD[2],

                       "splice_site_variant"= HK1[5], #RD[3],
                       "splicing"=HK1[5],
                       "extended_intronic_splice_region_variant"= HK1[3],
                       
                       "frameshift"= DM[2],
                       "truncating" = DM[1],
                       
                       "inframe_codon_gain"= DeR[4],
                       "inframe_codon_loss"= DeR[5],
                       "inframe_variant"= DeR[9],
                       "inframe"= DeR[4],
                       
                       "multi_hit" = HK1[4],
                       "multi_hits" =  HK1[4],
                       "multi_muts"=  HK1[4],

                       "loss" = RD[12],
                       "Loss" = RD[12],
                       "LOSS" = RD[12],
                       
                       
                       "gain" = MN[9],  # TP[5],
                       "Gain" = MN[9],
                       "GAIN" = MN[9],
                       
                       # "loh" = RD[11], # BALL
                       # "cnLOH"=  RD[11],
                       # "cnloh"=  RD[11],
                       
                       "loh" = "#91987A", #FH[2], #TP[5], # BALL FH[1],
                       "cnLOH"= "#91987A",# FH[2],
                       "cnloh"= "#91987A",# FH[2],
                       
                       "trans"= LK[1], #RD[8],
                       "TRA"= LK[1],
                       
                       "INV"= RD[9],
                       "inv"= RD[9],
                       
                       
                       "fusion"= RD[4],
                       "FUS"= RD[4],
                       
                       "SNV" = RD[4],
                       "INDEL"= RD[2],
                       
                       "iso"= "red", # CAS2[10],
                       "ISO" ="red", # CAS2[10],
                       
                       "other"= FH[7]
  )
  
  # browser()
  
  mut.colors <- unlist(mut.colors)
  
  # Define specific names to exclude from alpha modification
  #=============================================================================
  specific_names <- c("biallelic", "cnloh", #"gain",
                      "stop_gained", "stop_retained_variant","stop_gain",
                      "splice_site_variant", "splicing",
                      "multi_hit","multi_hits", "multi_muts",
                      "trans","trans","inv")

  # Apply alpha only to colors not in `specific_names`
  mut.colors[!tolower(names(mut.colors)) %in% specific_names] <- add.alpha(mut.colors[!(tolower(names(mut.colors)) %in% specific_names)], 0.8)
  
  # specific_names <- c("gain", "loss", "cnloh")
  # # only mut these (MPN)
  # # Apply alpha only to colors not in `specific_names`
  # mut.colors[tolower(names(mut.colors)) %in% specific_names] <- add.alpha(mut.colors[(tolower(names(mut.colors)) %in% specific_names)], 0.6)
  
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
    # "cnloh" = RD[11],
    "cnloh" = RD[3],
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
  
  complex.colors <- list("complex"= MN[1],"not complex"= MN[3],"not available"="ghostwhite")
  # complex.colors <- add.alpha(complex.colors, alpha = 0.6)
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
  # browser()
  response.colors <- list("persistent"="#016c59", 
                          "partial response"= "#1c9099", 
                          "non-responder"= "#67a9cf", 
                          "stable disease"= "#bdc9e1", 
                          "Stable disease" = RD[11], 
                          "responder"= "#df65b0", 
                          "CR"= MN[2], 
                          "CR-i"= HK1[1],
                          "CR/CRi" = MN[3],
                          "PR/stable disease"= RD[12], 
                          "PR"= "#b2df8a", 
                          "NR"= "#DA2310", 
                          "N/A"=  DM[8], 
                          "N/E"=  DM[8],
                          "NA" =  DM[8]
  )
  
  response.colors <- unlist(response.colors, use.names = TRUE)
  response.colors <- prettyGraphs::add.alpha(response.colors, 0.6)
  
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
  # == GENDER cols ----
  ############################### 
  
  GENDER <- list("MALE"= RD[9], 
                 "FEMALE"= MN[9], 
                 "UNKNOWN"= DM[6],
                 "NA"= DM[6]
  )
  GENDER <- unlist(GENDER)
  
  
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
  
  ###############################
  # == Complex.Karyotype color
  ###############################
  
  # Complex.Karyotype <- c(
  #   "complex"      = "firebrick",
  #   "not complex"  = MN[3],
  #   "not available"= "lightgray"
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
    "Singleton" = CAS1[4],
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
    "ZNF384 Group" = RD[9] , #DeR[6]
    "ZNF384-like" = RD[10],
    
    "Unclassified" = DM[7],
    "Low Quality" = DM[7],
    "failed"= DM[10],
    "Failed"= DM[10],
    "FAILED"= DM[10],
    "NoData"= "#FFFFFF", # white
    "ABL-class"= HK3[1],
    "CRLF2" = HK3[2],
    "P2RY8-CRLF2"= HK3[3],
    "Complex karyotype" = DM[4],
    "IGH-ID4" = DeR[2],
    "ZEB2" = DeR[3]
  )
  
  # Convert the list to a named vector if needed
  ALL.SUBTYPE <- unlist(subtype_colors)
  ALL.SUBTYPE <- ALL.SUBTYPE[sort(names(ALL.SUBTYPE))]
  
  ALL.SUBTYPE <- unlist(ALL.SUBTYPE)
  
  ALL.SUBTYPE <- add.alpha(ALL.SUBTYPE, alpha = .6)
  
  ####################################
  # == ALL CNV.WGS.CNVS.CALLS cols ----
  ####################################
  
  CNV.WGS.CNVS.CALLS <- list("#N/A"= "white", 
                             "CNVs"= MN[2], 
                             "No CNVs"= MN[7],
                             "oncoVAF/No CNVs"= MN[9])
  
  CNV.WGS.CNVS.CALLS <- unlist(CNV.WGS.CNVS.CALLS)
  
  ####################################
  ## ALL George's global clusters  ----
  ####################################
  
  # cl2.colors <- c(RD[1:12],FH[2:3]) # DM[8]
  # names(cl2.colors) <- paste0("cl.",1:14)
  # cl2.colors <- add.alpha(cl2.colors, alpha = 0.65)
  
  ####################################
  ## ALL George's global clusters  ----
  ####################################
  cl2.colors <- c(
    "cl.1"  = "red",  # red
    "cl.2"  = "#377EB8",  # blue
    "cl.3"  = "#4DAF4A",  # green
    "cl.5"  = "#FF7F00",  # orange
    "cl.6"  = "#984EA3",  # purple
    "cl.7"  = "#A65628",  # brown
    "cl.8"  = "#F781BF",  # pink
    "cl.9"  = "#999999",  # gray
    "cl.10" = "#00CED1",  # dark turquoise
    "cl.14" = "#FFD700"   # gold
  )
  
  ###############################
  # == Disease cols ----
  ############################### 
  
  OK  <- met.brewer("OKeeffe1", type = "discrete")
  
  DISEASE <- list("MF"= RD[9], 
                 "PV"= RD[5], 
                 "AML"= MN[6],
                 "NA"= DM[6],
                 
                 "Chronic.Phase" = TP[4], 
                 "Chronic Phase" = TP[4], 
                 "Chronic" = TP[4], 
                 "Chronic MPN"= TP[4],
                 
                 "MPN AP/BP" = TP[1],
                 
                 "Blast.Phase"= RD[2]
                 

                 # "MPN AP/BP" = MN[4],
                 # "Chronic" = MN[8], 
                 # "Chronic MPN" = MN[8]
  )
  
  DISEASE <- unlist(DISEASE)
  DISEASE <- add.alpha(DISEASE, alpha = 0.8)
  
  
  DISEASE.PROGRESSION <- list("ET --> AML"= OK[2], 
                              "ET --> MF --> AML"= OK[4])
  DISEASE.PROGRESSION <- unlist(DISEASE.PROGRESSION)
  # DISEASE.PROGRESSION <- add.alpha(DISEASE.PROGRESSION, alpha = 0.8)
  ###############################
  # == Disease cols ----
  ############################### 
  
  INITIAL.MPN.DX <- list("MF"= "#F4A460", 
                  "PV"= "#9DB07A",
                  "ET"="#C0392B",
                  "NA"= DM[6]
  )
  
  INITIAL.MPN.DX <- unlist(INITIAL.MPN.DX)
  ###############################
  # == CNACS colros ----
  ############################### 
  
  COPY.NUMBER.STATUS <- list("QC pass"= FH[4], 
                     "QC fail" = FH[8],
                     "QC borderline"= FH[3],
                     "VERY NOISY" = "red",
                     "NOISY"= "green"
                     
  )
  COPY.NUMBER.STATUS <- unlist(COPY.NUMBER.STATUS)
  COPY.NUMBER.STATUS <- add.alpha(COPY.NUMBER.STATUS, alpha = 0.8)
  
  #====================================================
  
  return(list(mut.colors= mut.colors, 
              default.mut.colors= mut.colors.def,
              #cyto.colors=cyto.colors,  
              response.colors= response.colors,
              eln.molecular.response.colors = eln.molecular.response.colors,
              COPY.NUMBER.STATUS = COPY.NUMBER.STATUS,
              nice.cols.A= nice.cols.A,
              therapy.colors=therapy.colors,
              path.colors= path.colors,
              
              ALL.SUBTYPE= ALL.SUBTYPE, ### < ALL-specific
              global.clusters = cl2.colors,  ### < ALL-specific
              CNV.WGS.CNVS.CALLS = CNV.WGS.CNVS.CALLS,
              
              GENDER= GENDER,
              DISEASE=DISEASE,
              INITIAL.MPN.DX= INITIAL.MPN.DX,
              PURITY= purity_colors,
              EE= purity_colors,
              COMPLEX.KARYOTYPE= complex.colors,
              
              DISEASE.PROGRESSION= DISEASE.PROGRESSION
              # Complex.Karyotype = Complex.Karyotype
  )) 
  
}