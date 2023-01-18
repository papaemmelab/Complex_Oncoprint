heatmap_colors <-  function() {
  
  library(prettyGraphs)
  #==================================================================================
  #   Written by Noushin Farnoud, Jul 2018
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
  # 3182bd
  
  # mut.colors <- c("#74add1", #missense
  #                 "#000000", #stop_gain
  #                 "#d73027", #frameshift_indel
  #                 "#fdae61", #splice_site_variant
  #                 "#a6d96a", #inframe_indel
  #                 "#bdbdbd", #unknown
  #                 "#636363", #other
  #                 "#636363", #other_snvs
  #                 "#FF7400FF") #complex
  
  mut.colors <- c("#3182bd","#000000", "#de2d26", "#8460db","#598238","#9ecae1","#bdbdbd","#636363","#FF7400FF")
  
  mut.colors <- add.alpha(mut.colors, alpha = .7)
  
  names(mut.colors) <- c("missense","stop_gain","frameshift_indel","frameshift_insersion","splice_site_variant","inframe_indel","unknown","other","complex")
  
  # names(mut.colors) <- c("missense","stop_gain","frameshift_indel","splice_site_variant","inframe_indel","unknown","other","other_snvs","complex")
  
  
  ###############################
  # == Flow.MRD.response.Colors
  ###############################  
  
  #show_col(hue_pal()(3))
  
  # Flow.MRD.response.Colors <- c("#d4b9da", "#980043", "#dd1c77","#df65b0","#f1eef6")
  # names(Flow.MRD.response.Colors) <- c("Not Performed","Persistent Disease","MRD Positive","MRD Negative","Inevaluable")
  
  ###############################
  # == ELN colors
  ###############################  
  
  eln.molecular.response.colors <- c("#80b1d3","#b2df8a","#fb8072","#d9d9d9")
  names(eln.molecular.response.colors) <- c("CR","PR","NR","N/E")
  
  ###############################
  # == path.colors colors
  ############################### 
  
  path.colors <- c("#80b1d3","#fb8072","#d9d9d9")
  names(path.colors) <- c("HR","NR","N/E")
  
  ###############################
  # == Response colors
  ###############################  
  
  response.colors <- c("#016c59", "#1c9099", "#67a9cf","#bdc9e1","#df65b0",
                       "#80B1D3", #CR
                       "#80B1D3", #CR-i
                       "#F39C12", #PR
                       "#DA2310", #NR
                       "#bdbdbd", #NA
                       "#bdbdbd" #N/E
                       )
  
  names(response.colors) <- c("persistent","partial response","non-responder","stable disease","responder","CR","CR-i","PR","NR","N/A","N/E")
  
  # response.colors <- c("#016c59", #persistent
  #                      "#1c9099", #partial response
  #                      "#67a9cf", #non-responder
  #                      "#bdc9e1", #stable disease
  #                      "#df65b0", #responder
  #                      "#80b1d3", #CR
  #                      "darkorange", #PR
  #                      "#fb8072", #NR
  #                      "#ccebc5", #NR–CR
  #                      "#fccde5", #NR–PR
  #                      "#fb8072", #NR–NR
  #                      "#d9d9d9") #N/E
  # 
  # names(response.colors) <- c("persistent","partial response","non-responder","stable disease","responder","CR","PR","NR","NR–CR","NR–PR","NR–NR","N/E")
  
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
  
  
  
  ###############################
  # == Define CNV colors
  ###############################  
  
  cyto.colors <- c(
    "#3288bd", #AMP
    '#d53e4f', #DEL
    "#f46d43", #LOH
    "#6a3d9a", #complex_karyotype
    "#66c2a5", #inconclusive
    '#D3D3D3', #INV
    '#6a3d9a', #COMPX
    "#66c2a5", # TRANS
    '#cab2d6', #TDUPs
    '#ff7f00', #OTHER_SVs
    '#fdbf6f', #FUSION
    'grey', #"Karyotype Failure"
    '#fb9a99', #11q23
    '#33a02c', #t(9;11)
    '#b2df8a', #CBF inv(16)
    '#1f78b4', #inv(3)
    '#a6cee3', #+8/add(8)
    '#80b1d3', #add(18)
    '#fdae61', # ADD
    '#fee08b', # DER
    '#fee08b', # del(17p)
    "#66c2a5", #-7
    '#c7eae5' # Normal_Karyotype
  )
  
  # cyto.colors <- add.alpha(cyto.colors, alpha = .7)
  
  names(cyto.colors) <- c("AMP",
                          "DEL",
                          "LOH",
                          "complex_karyotype",
                          "inconclusive",
                          "INV",
                          "COMPX",
                          "TRANS",
                          "TDUP",
                          "OTHER_SVs",
                          "FUSION",
                          "Karyotype Failure",
                          "11q23",
                          "t(9;11)",
                          "CBF inv(16)", 
                          "inv(3)",
                          "+8/add(8)",
                          "add(18)",
                          "ADD",
                          "DER",
                          "del(17p)",
                          "-7",
                          "Normal_Karyotype") 
 
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