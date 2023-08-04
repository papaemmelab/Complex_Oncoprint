This is to create oncoprints of variants (mutations [required], copy number variants [optional], structural variants [optional]).
Look at the example file and documentation for an example of how to call the function. 

Simple example: 

    source("./generate_complex_oncoprint.R")
  
    ht <-  generate_complex_oncoprint(muts= MUT, 

                                    show.sample.names = TRUE, show.border= FALSE,

                                    min.freq= 1, show.title= TRUE, title.str=  "Example 1 - YOHOOOO oncoprint with NO ERROR!", 

                                    save.name= "Example_1.A",

                                    save.path= "./example_oncoprints/") 
 

Feb 2019: features added to show RESPONSE, and SAMPLE.SOURCE (for example you can choose B-cell, T-cell, or you can even pass on the disease type such as PV, ET, etc.). Also, for timepoint data, you have the option to show a heatbar of samples that belong to the same patient. 
Sep 2019: new features added to allow defining the second annotation ribbon as wish (e.g., RACE or DISEASE). Other visualization added and will be explained in an analysis session.
You can install/upgrade using:

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")


Current recognized terms for mutation/CNV/SV types (col EFFECT) are:

 [1] "amp"                          "amplification"                "CN-del"                       "CN-gain"                      "complex"                     
 [6] "complex_change_in_transcript" "complex_karyotype"            "del"                          "deletion"                     "frame_shift_del"             
[11] "frame_shift_ins"              "frameshift_del"               "frameshift_indel"             "frameshift_insertion"         "frameshift_variant"          
[16] "fus"                          "fusion"                       "gain"                         "in_frame_del"                 "in_frame_ins"                
[21] "inconclusive"                 "inframe_codon_gain"           "inframe_codon_loss"           "inframe_deletion"             "inframe_indel"               
[26] "inframe_insersion"            "inframe_variant"              "initiator_codon_change"       "inv"                          "inversion"                   
[31] "karyotypic_abnormal"          "LOH"                          "loss"                         "missense"                     "missense_codon"              
[36] "missense_mutation"            "N_E"                          "N/A"                          "N/E"                          "non_synonymous_codon"        
[41] "nonsense_mutation"            "nonstop_mutation"             "normal"                       "normal_karyotype"             "other_cnvs"                  
[46] "other_snvs"                   "other_svs"                    "rearr"                        "rearrangement"                "rearrangements"              
[51] "splice_site"                  "splice_site_variant"          "stop_gain"                    "stop_gained"                  "stop_lost"                   
[56] "stop_retained_variant"        "tandem dup"                   "tandem duplications"          "tandem_duplications"          "trans"                       
[61] "translation_start_site"       "translocation"                "unavailable"                  "unknown"  

IMPORTANT: 

Please contact rahnaman@mskcc.org if you had any mutation-type that wanted to have a specific representation for it in Oncoprint.

