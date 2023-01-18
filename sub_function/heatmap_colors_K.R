heatmap_colors_K <-  function() {
  
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
  
  mut.colors <- c("lightcyan3",
                  "thistle4",
                  "lightsalmon",
                  "slategray",
                  "cornflowerblue",
                  "#3182bd", #missense
                  "#000000", #stop_gain
                  "#de2d26", #frameshift_indel
                  "#7CAE00", #splice_site_variant
                  "#80B1D3", #inframe_indel
                  "#bdbdbd", #unknown
                  "#636363", #other
                  "#636363", #other_snvs
                  "#FF7400FF") #complex
  
  #mut.colors <- add.alpha(mut.colors, alpha = .8)
  
  names(mut.colors) <- c("one","two","four","eight","none","missense","stop_gain","frameshift_indel","splice_site_variant","inframe_indel","unknown","other","other_snvs","complex")
  
  
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
                       "#F39C12", #PR
                       "#DA2310", #NR
                       "#bdbdbd" #NA
                       )
  names(response.colors) <- c("persistent","partial response","non-responder","stable disease","responder","CR","PR","NR","N/A")
  
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
  
  cyto.colors <- c(#'#F8766D', #AMP
                  "lightsteelblue", #CNA
                   "lightsalmon", #AMP
                   'lightsteelblue', #DEL
                   "#C6CDF7", #LOH
                   '#6a3d9a',
                   '#cab2d6',
                   '#ff7f00',
                   '#fdbf6f',
                   '#fdbf6f',
                   'darkorange', #fusion
                   '#fb9a99',
                   '#33a02c',
                   '#b2df8a',
                   '#1f78b4',
                   '#a6cee3',
                   '#80b1d3',
                   '#b3de69',
                   '#fccde5',
                   '#bc80bd'
  )
  
  names(cyto.colors) <- c("CNA",
                          "AMP",
                          "DEL",
                          "LOH",
                          "INV",
                          "COMPX",
                          "TDUP",
                          "OTHER_SVs",
                          "other_svs",
                          "FUS",
                          "TRANS",
                          "Karyotype Failure",
                          "11q23",
                          "t(9;11)",
                          "CBF inv(16)", 
                          "inv(3)",
                          "+8/add(8)",
                          "del(17p)",
                          "-7") 
 
  # this was initially used for "J Grinfeld et al. Classification" ribbon in 157 manuscript and I liked the combo (TP53 mutation etc); the names must be adjusted based on the ribbon features
  
  nice.cols.A <- c("#C1DAD6", #TP53 mutation"
                  "#CCFFCC", #"CALR mutation"
                  "#FFCF79", #"MPL mutation"
                  "#B7AFA3", #"homozygous JAK2/NFE2 mutation"
                  "#E8D0A9", #"heterozygous JAK2"
                  "#CCCCCC", #"Other drivers"
                  "#6D929B", #"Chromatin/Spliceosome/RAS mutation
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