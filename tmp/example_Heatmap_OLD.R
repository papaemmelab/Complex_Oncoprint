#==================================================================================
#   Written by Noushin Farnoud, Jul 2018. Last modified Sep 2019
#----------------------------------------------------------------------------------
#   This is an example of how to run Oncoprint-Scripts for a set of mutations/cnvs/svs.
# 
#   Contact Noushin Farnoud (rahnaman@mskcc.org) if you faced any error.
#
#   Many thanks to Teja for proving the mutation data that is used in this example.
# 
#   See also generate_oncoprint.R.
#==================================================================================

rm(list=ls()) #clear all vars 

library(reshape)
library(gsheet)

#########################
# Initialize
#########################

setwd("/Users/rahnaman/analysisTools/src/tools/R/Oncoprint")

teja.url <- "https://docs.google.com/spreadsheets/d/1yZzCSA38mZbGyeCjlzh6xups2t0S6y_JY0jcTwDKp5A/edit#gid=1786134805"

min.freq <- 1

#########################
# Load SUBs and INDELs 
#########################

A <- read.csv(construct_download_url(teja.url, format = "csv", sheetid = "376689594"), stringsAsFactors=FALSE)

#======================================================
# === choose LIKELY/ONCOGENIC subs for the oncoprint
#======================================================
SUBs <- subset(A, MANUAL_ANNOTATION %in% c("ONCOGENIC","LIKELY"))
rm(A)
#======================================================
# === choose LIKELY/ONCOGENIC indels for the oncoprint
#======================================================
B <- read.csv(construct_download_url(teja.url, format = "csv", sheetid = "407131731"), stringsAsFactors=FALSE)
INDELs <- subset(B, MANUAL_ANNOTATION %in% c("ONCOGENIC","LIKELY"))
rm(B)
#========================================================================
# === Merge indels and subs to prepare a unique df of mutation data MUTs
#========================================================================
my.cols <- intersect(colnames(INDELs), colnames(SUBs))
MUTs <- rbind(INDELs[,my.cols],SUBs[,my.cols])
MUTs[,52] <- NULL

rm(INDELs,SUBs)

###########################################################################
# Load and prepare CNVs : Required columns (TARGET_NAME, VAR_ID, VT)
###########################################################################

D <- read.csv(construct_download_url(teja.url, format = "csv", sheetid = "716153617"), stringsAsFactors=FALSE)

cnvs <- melt(D)

ix <- which(cnvs$value==0)
cnvs <- cnvs[-ix,]
cnvs$VT <- NA

ix.del <- grep("del",cnvs$variable)
cnvs$VT[ix.del] <- "DEL"

ix.amp <- grep("amp",cnvs$variable)
cnvs$VT[ix.amp] <- "AMP"

cnvs <- cnvs[-which(is.na(cnvs$VT)),]

rm(D)

colnames(cnvs)[colnames(cnvs) == 'Sample'] <- 'TARGET_NAME'
colnames(cnvs)[colnames(cnvs) == 'variable'] <- 'VAR_ID'
colnames(cnvs)[colnames(cnvs) == 'VT'] <- 'EFFECT'

##########################################################
# Load SVs : Required columns (TARGET_NAME, VAR_ID, VT)
##########################################################

C <- read.csv(construct_download_url(teja.url, format = "csv", sheetid = "356721653"), stringsAsFactors=FALSE)

SVs <- subset(C, gene!="_")

colnames(SVs)[colnames(SVs) == 'svtype'] <- 'EFFECT'
colnames(SVs)[colnames(SVs) == 'sample'] <- 'TARGET_NAME'
colnames(SVs)[colnames(SVs) == 'gene'] <- 'VAR_ID'

######################################
# == Add annotation ribbon info ====
######################################

my.table <- read.csv(construct_download_url(teja.url, format = "csv", sheetid = "773700969"), stringsAsFactors=FALSE)

#==========================
# == add made-up Response
#==========================
my.table$RESPONSE <- "NA"
my.table$RESPONSE[1:20] = "CR"
my.table$RESPONSE[21:30] = "PR"
my.table$RESPONSE[31:76] = "NR"

my.table$SAMPLE.SOURCE <- my.table$Race
my.table$TARGET_NAME <- my.table$SAMPLE.ID

#########################
# Run Oncoprint examples
#########################

savePath <- "./example_oncoprints/"

###################################################
# EXAMPLE (1): ====
# Successful run with minimal input
###################################################

cat(paste("\n>>>>> Staring Example 1 ... \n"))

source("./generate_complex_oncoprint.R")

ht <- generate_complex_oncoprint(muts= MUTs, svs= NULL, cnvs=cnvs, min.freq= 1, save.path= save.path, title.str= "Example 1 - MYTYPE p212")


# this saved the jpg version of the heatmap at the specified location; 
# you can load the heatmap and save the pdf ver using ht$ht.obj

cat(paste("\nFinished Example 1. \n"))
cat(paste("====================================================\n\n"))

###################################################
# EXAMPLE (2): ====
###################################################

cat(paste("\n>>>>> Staring Example 2. Show genes with >= 2 freq, Add Response and a new added ribbon (for Race). Control visualizations ... \n"))

source("./generate_complex_oncoprint.R")

ht <-  generate_complex_oncoprint(muts= MUTs, cnvs= cnvs, svs= SVs ,  # define variants DFs
                                        
                                  cnvs.order= NULL, svs.order= NULL, muts.order= NULL, patients.order= NULL,   # define pre-set orders
                                  
                                  surval.data= NULL,
                                  
                                  show.response= TRUE, response.order= NULL, show.another.banner=TRUE, banner.name= "Race", show.individuals= FALSE, show.individuals.legend= FALSE, show.survival= FALSE, # ******* user-defined ORDERs
                                  
                                  lookup.table= my.table, #******* pass lookup.table 
                                  
                                  show.sample.names = TRUE, show.border= FALSE, show.multis= TRUE, rem.empty= TRUE, # ******* what params to show in legend?
                                  
                                  heatmap.legend.side= "right", mut.legend.title.side= "topcenter", num.rows.heatmap.lgd= NULL, #******* HEATMAP.legend 
                                  
                                  annot.legend.side= "bottom", annot.title.side= "topleft", num.rows.annot.lgd= NULL,  #******* ANNOT.legend 
                                  
                                  min.freq= 2, show.title= TRUE, title.str=  "Notice the impact of setting show.multis to TRUE", 
                                  
                                  save.path= savePath, #******* title and save path
                                  
                                  cols.font= 22, rows.font= 22, pct.font= 22, legend.label.font= 20, legend.title.font= 25, fig.title.font= 22,  barplot.font= 10,  
                                  
                                  multis.dot.size = 0.8, #****FONTs: row.groupname.font is the same as rows.font
                                  
                                  right.w= 8, top.w= 5 , w=3200, h=2800,  #**** Sizes of barplots and fig 
                                  
                                  axis.side= "left")


cat(paste("\nFinished Example 3. Now I like to define the order of the response in the legend.\n"))
cat(paste("====================================================\n"))

source("./generate_complex_oncoprint.R")

ht <-  generate_complex_oncoprint(muts= MUTs, cnvs= cnvs, svs= SVs ,  # define variants DFs
                                  
                                  cnvs.order= NULL, svs.order= NULL, muts.order= NULL, patients.order= NULL,   # define pre-set orders
                                  
                                  surval.data= NULL,
                                  
                                  show.response= TRUE, response.order= c("CR","PR","NR"), show.another.banner=TRUE, banner.name= "Race", show.individuals= FALSE, show.individuals.legend= FALSE, show.survival= FALSE, # ******* user-defined ORDERs
                                  
                                  lookup.table= my.table, #******* pass lookup.table 
                                  
                                  show.sample.names = FALSE, show.border= FALSE, show.multis= TRUE, rem.empty= TRUE, # ******* what params to show in legend?
                                  
                                  heatmap.legend.side= "right", mut.legend.title.side= "topcenter", num.rows.heatmap.lgd= NULL, #******* HEATMAP.legend 
                                  
                                  annot.legend.side= "bottom", annot.title.side= "topleft", num.rows.annot.lgd= NULL,  #******* ANNOT.legend 
                                  
                                  min.freq= 3, show.title= TRUE, title.str=  "Notice the order of RESPONSE in the legend", 
                                  
                                  save.path= savePath, #******* title and save path
                                  
                                  cols.font= 22, rows.font= 22, pct.font= 22, legend.label.font= 20, legend.title.font= 25, fig.title.font= 22,  barplot.font= 10,  
                                  
                                  multis.dot.size = 0.8, #****FONTs: row.groupname.font is the same as rows.font
                                  
                                  right.w= 8, top.w= 5 , w=3200, h=2800,  #**** Sizes of barplots and fig 
                                  
                                  axis.side= "left")


cat(paste("\nFinished Example 3. Try playins with params and see the difference.\n"))
cat(paste("====================================================\n"))




