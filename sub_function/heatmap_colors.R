heatmap_colors <-  function() {
  
  library(prettyGraphs)
  
  #==================================================================================
  #   Written by Noushin Farnoud, Jul 2018, Latest Update June 2020
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

  mut.colors <- c("#3182bd","#000000", "#de2d26", "#8460db","#598238","#9ecae1","#bdbdbd","#ffeda0","#FF7400FF")
  
  mut.colors <- add.alpha(mut.colors, alpha = .9)
  
  names(mut.colors) <- c("missense","stop_gain","frameshift_indel","initiator_codon_change",
                         "splice_site_variant","inframe_indel","unknown","other_snvs","complex")
  
  ###############################
  # == Define CNV/Cytogenetics colors
  ###############################  
  
  cyto.colors <- c(
    # "#3288bd", # AMP
    "#3C33FF", #AMP
    # '#d53e4f', # DEL
    "#FB3F28", #DEL
    "#f46d43", # LOH
    "#6a3d9a", # complex_karyotype
    # '#feb24c', # INV
    "#c51b8a", # INV
    '#FF7400FF', #complex
    "#66c2a5", # TRANS
    "darkgreen", # REARRANGEMENT
    '#cab2d6', # TDUPs
    '#c51b8a', # OTHER_SVs
    '#fa9fb5', # OTHER_CNVs
    '#fdbf6f', # FUSION
    'grey',    # N/E or inconclusive
    '#fdae61', # ADD
    '#fee08b', # DER
    # '#c7eae5'  # Normal_Karyotype
    '#a8ddb5'  # Normal_Karyotype
  )
  
  # cyto.colors <- add.alpha(cyto.colors, alpha = .7)
  
  names(cyto.colors) <- c("amp",
                          "del",
                          "loh",
                          "complex_karyotype",
                          "inv",
                          "complex",
                          "trans",
                          "rearr",
                          "tdup",
                          "other_svs",
                          "other_cnvs",
                          "fusion",
                          "NA", # N/E
                          "add",
                          "der",
                          "normal") 

  ###############################
  # == ELN colors (special case)
  ###############################  
  
  eln.molecular.response.colors <- c("#80b1d3","#b2df8a","#fb8072","#d9d9d9")
  names(eln.molecular.response.colors) <- c("CR","PR","NR","N/E")
  
  #############################################
  # == pathology.report colors (special case)
  #############################################
  
  path.colors <- c("#80b1d3","#fb8072","#d9d9d9")
  names(path.colors) <- c("HR","NR","N/E")
  
  ###############################
  # == Response colors
  ###############################  
  
  response.colors <- c("#016c59", # persistent
                       "#1c9099", # partial response
                       "#67a9cf", # non-responder
                       "#bdc9e1", # stable disease
                       "#df65b0", # responder
                       "#80B1D3", # CR
                       "#80B1D3", # CR-i
                       "#F39C12", # PR
                       "#DA2310", # NR
                       "#bdbdbd", # NA
                       "#bdbdbd"  # N/E
                       )
  
  names(response.colors) <- c("persistent","partial response","non-responder","stable disease","responder","CR","CR-i","PR","NR","N/A","N/E")
  
  ###############################
  # == therapy colors
  ###############################  
  
  therapy.colors <- c("#80b1d3", #No resistance/intolerance
                       "#80b1d3", #No.Resistance
                       "#fee090", #Intolerant
                       "#fb8072", #Resistant
                       "#df65b0", #MRD Negative
                       "#d9d9d9", #Inevaluable
                       "#d9d9d9") #N/E
  
  names(therapy.colors) <- c("No resistance/intolerance","No.Resistance","Intolerant","Resistant","MRD Negative","Inevaluable")
  
   ###############################################################
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
  
  
  #cyto.colors <- add.alpha(cnv.colors, alpha = .8)
  
  return(list(mut.colors= mut.colors, 
              cyto.colors=cyto.colors,  
              response.colors= response.colors,
              nice.cols.A= nice.cols.A,
              therapy.colors=therapy.colors,
              eln.molecular.response.colors,
              path.colors= path.colors
              )) 
  
}