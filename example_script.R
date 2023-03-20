#==================================================================================
#   Written by Noushin Farnoud, Jul 2018. Last modified July 2023
#----------------------------------------------------------------------------------
#   This is an example of how to run Oncoprint-Scripts for a set of mutations/cnvs/svs.
# 
#   Contact Noushin Farnoud (rahnaman@mskcc.org) if you faced any error.
#
#   Many thanks to Teja for proving the SV example data that is used in this example.
# 
#   See also generate_complex_oncoprint.R
#==================================================================================

rm(list=ls()) #clear all vars 

library(reshape)
library(gsheet)

library(readxl)

library(plyr)
library(dplyr)

#########################
# Initialize
#########################

setwd("/Users/rahnaman/Documents/Oncoprint/")

muts.file <- "./example_data/Mutations.xlsx" # 1st dataframe (Noushin's)

teja.url <- "https://docs.google.com/spreadsheets/d/1yZzCSA38mZbGyeCjlzh6xups2t0S6y_JY0jcTwDKp5A/edit#gid=1786134805" # 2nd dataframe (Teja's)

savePath <- "./example_oncoprints/"

min.freq <- 1 

#########################
# Load SUBs and INDELs 
#########################

MUT <- read_excel(muts.file) %>% filter(Noush_annot2 %in% c("ONCOGENIC","LIKELY"))

###########################################################################
# EXAMPLE 1.A: Generate Simple Oncoprint of Mutations ----
###########################################################################

cat(paste("\n>>>>> Staring Example 2. Show genes with >= 2 freq, Add Response and a new added ribbon (for Race). Control visualizations ... \n"))

source("./generate_complex_oncoprint.R")

ht <-  generate_complex_oncoprint(muts= MUT, 
                                  
                                  show.sample.names = TRUE, show.border= FALSE,
                                  
                                  min.freq= 1, show.title= TRUE, title.str=  "Example 1 - YOHOOOO I managed to get this without any error, but these default font sizes are NOT good.", 
                                  
                                  save.name= "Example_1.A", #******* NAME of the saved oncoprint fig
                                  
                                  save.path= "./example_oncoprints/", #******* DIR of the saved oncoprint fig
)

###########################################################################
# EXAMPLE 1.B: Run previ. Example with better fonts ----
###########################################################################

cat(paste("\nFinished Example 1. Now I like to improve the size of the fonts.\n"))

ht <-  generate_complex_oncoprint(muts= MUT, 
                                  
                                  show.sample.names = TRUE, show.border= FALSE,
                                  
                                  min.freq= 1, show.title= TRUE, title.str=  "Example 1.B - YOHOOOO - GOOD FONTs!", 
                                  
                                  save.name= "Example_1.B", #******* NAME/DIR of the saved oncoprint fig
                                  
                                  save.path= "./example_oncoprints/", #******* DIR of the saved oncoprint fig
                                  
                                  cols.font= 18, 
                                  
                                  rows.font= 18, 
                                  
                                  pct.font= 20, 
                                  
                                  legend.label.font= 20, 
                                  
                                  legend.title.font= 25, 
                                  
                                  fig.title.font= 25,  
                                  
                                  barplot.font= 20,  
                                  
                                  right.w= 8, top.w= 5 , w=2000, h=1600,  #**** Sizes of barplots and fig 
                                  
                                  axis.side= "left")


cat(paste("\nFinished Example 1. Now I like to define the order of the response in the legend.\n"))

##############################################################################################
# EXAMPLE 5: How to override default clustering with an input order of patients (cols)? ---- 
##############################################################################################

my.pts.order <- c("E-H-131710-T1-1-D1-1", "E-H-131712-T1-1-D1-1", "E-H-131721-T1-1-D1-1", "E-H-131713-T1-1-D1-1", "E-H-131717-T1-1-D1-1", "E-H-131719-T1-1-D1-1", "E-H-131720-T1-1-D1-1",
              "E-H-131714-T1-1-D1-1", "E-H-131711-T1-1-D1-1",  "E-H-131716-T1-1-D1-1","E-H-131715-T1-1-D1-1", "E-H-131718-T1-1-D1-1")

ht.order <-  generate_complex_oncoprint(muts= MUT, 
                                  
                                  show.sample.names = TRUE, patients.order = my.pts.order,
                                  
                                  min.freq= 1, show.title= TRUE, title.str=  "Example 6 user-specified patients.order", 
                                  
                                  save.name= "Example_5", 
                                  
                                  save.path= "./example_oncoprints/", 
                                  
                                  cols.font= 18, 
                                  
                                  rows.font= 18, 
                                  
                                  pct.font= 20, 
                                  
                                  legend.label.font= 20, 
                                  
                                  legend.title.font= 25, 
                                  
                                  fig.title.font= 25,  
                                  
                                  barplot.font= 20,  
                                  
                                  right.w= 8, top.w= 5 , w=2000, h=1600,  #**** Sizes of barplots and fig 
                                  
                                  axis.side= "left")

cat(paste("\nFinished USER-defined ORDER example.\n"))

###########################################################################
# Load and prepare CNVs : Required columns (TARGET_NAME, VAR_ID, VT)
###########################################################################
rm(list=ls()) #clear all vars 

savePath <- "./example_oncoprints/"

load("./example_data/Example2.RData")

############################################################
# EXAMPLE (2): Oncoprint with MUT + CNV + SV  ====
############################################################

cat(paste("\n>>>>> Staring Example 2 ... \n"))

source("./generate_complex_oncoprint.R")

ht <- generate_complex_oncoprint(muts= MUTs, svs= SVs, cnvs=cnvs, min.freq= 1, 
                                 save.name= "Example_2",
                                 save.path= save.path, title.str= "Example 2 - Mutuations + CNVs + SVs",
                                 heatmap.legend.side= "bottom", #**change the default location of mutation legend
                                 num.rows.heatmap.lgd= 2, 
                                 mut.legend.title.side = "leftcenter",
                                 pct.font = 18,
                                 legend.label.font = 20,
                                 legend.title.font = 18,
                                 barplot.font = 18,
                                 right.w = 5)


# this saved the jpg version of the heatmap at the specified location; 
# you can load the heatmap and save the pdf ver using ht$ht.obj

cat(paste("================================\n\n"))

###################################################
# EXAMPLE (3): Oncoprint with RESPONSE annotation
###################################################

cat(paste("\n>>>>> Staring Example 3. Show genes with >= 2 freq, Add Response and a new added ribbon (for Race). Control visualizations ... \n"))

source("./generate_complex_oncoprint.R")

ht <-  generate_complex_oncoprint(muts= MUTs, cnvs= cnvs, svs= SVs ,  # define variants DFs
                                  
                                  save.name= "Example_3",
                                        
                                  cnvs.order= NULL, svs.order= NULL, muts.order= NULL, patients.order= NULL,   # define pre-set orders (if needed to)
                                  
                                  show.response= TRUE, response.order= NULL,  ## ***** show RESPONSE 
                                  
                                  lookup.table= my.table, #******* REQUIRED: PASS lookup.table 
                                  
                                  show.sample.names = TRUE, show.border= FALSE, show.multis= FALSE, rem.empty= TRUE, # ******* params to show in legend
                                  
                                  heatmap.legend.side= "right", mut.legend.title.side= "topcenter", num.rows.heatmap.lgd= NULL, #******* HEATMAP.legend 
                                  
                                  annot.legend.side= "bottom", annot.title.side= "topleft", num.rows.annot.lgd= NULL,  #******* ANNOT.legend 
                                  
                                  min.freq= 2, show.title= TRUE, title.str=  "Notice the impact of setting show.multis to TRUE", 
                                  
                                  save.path= savePath, #******* title and save path
                                  
                                  cols.font= 22, rows.font= 22, pct.font= 22, legend.label.font= 20, 
                                  
                                  legend.title.font= 25, fig.title.font= 22,  barplot.font= 20,  
                                  
                                  multis.dot.size = 0.8, #****FONTs: (for noush => 'row.groupname.font' is the same as rows.font)
                                  
                                  right.w= 8, top.w= 5 , w=3200, h=2800,  #**** Sizes of barplots (right.w and top.w) and saved figure (w,h) 
                                  
                                  axis.side= "left")

cat(paste("\nFinished Example 3. Now I like to also visualize my multi-hit genes.\n"))
cat(paste("===========================================================================\n"))

#########################################################
# EXAMPLE (4): Example 3 + multi-hit Visualization ====
#########################################################

source("./generate_complex_oncoprint.R")

ht <-  generate_complex_oncoprint(muts= MUTs, cnvs= cnvs, svs= SVs ,  # define variants DFs
                                  
                                  cnvs.order= NULL, svs.order= NULL, muts.order= NULL, patients.order= NULL,   # define pre-set orders
                                  
                                  surval.data= NULL,
                                  
                                  show.response= TRUE, response.order= c("CR","PR","NR"), show.another.banner=TRUE, banner.name= "Race", show.individuals= FALSE, show.individuals.legend= FALSE, show.survival= FALSE, # ******* user-defined ORDERs
                                  
                                  lookup.table= my.table, #******* pass lookup.table 
                                  
                                  show.sample.names = FALSE, show.border= FALSE, show.multis= TRUE, rem.empty= TRUE, # ******* what params to show in legend?
                                  
                                  heatmap.legend.side= "right", mut.legend.title.side= "topcenter", num.rows.heatmap.lgd= NULL, #******* HEATMAP.legend 
                                  
                                  annot.legend.side= "bottom", annot.title.side= "topleft", num.rows.annot.lgd= NULL,  #******* ANNOT.legend 
                                  
                                  min.freq= 2, show.title= TRUE, title.str=  "Notice the order of RESPONSE in the legend", 
                                  
                                  #******* title and save path
                                  
                                  save.name = "Example_4",
                                  
                                  save.path= savePath, 
                                  
                                  cols.font= 22, rows.font= 22, pct.font= 22, legend.label.font= 20, legend.title.font= 25, fig.title.font= 22,  barplot.font= 10,  
                                  
                                  multis.dot.size = 0.8, #****FONTs: row.groupname.font is the same as rows.font
                                  
                                  right.w= 8, top.w= 5 , w=3200, h=2800,  #**** Sizes of barplots and fig 
                                  
                                  axis.side= "left")


cat(paste("============================================================================"))
cat(paste("\nEnd of Tutorial. Try playing with the params and observe the difference.\n"))
cat(paste("============================================================================\n"))




