rm(list=ls()) #clear all vars 

library(ggplot2)
library(reshape)
library(data.table)
library("RColorBrewer")
library(stargazer)
library(gridExtra)
library(scales)
library(stringr)
library(readxl)
library(writexl)
library(gsheet)
library(readr)
library(png)
library(patchwork)
library(cowplot)
library(gtsummary)

library(plyr)
library(dplyr)

library(MetBrewer)

options(warn=-1) 


GD_path <- Sys.getenv("GD_path")
func_path <- Sys.getenv("func_path")

prj.path <- file.path(GD_path, "leukgen","projects","ALL_UK")

Fig.Dir.main <- file.path(prj.path, "Figures",paste0("Multi_Panel_",Sys.Date()))

dir.create(Fig.Dir.main, showWarnings=FALSE)

#===============================================================

Fig.Dir.2 <- file.path(Fig.Dir.main, "Oncoprint")

dir.create(Fig.Dir.2, showWarnings=FALSE)
dir.create(file.path(Fig.Dir.2,"gene_centric"), showWarnings=FALSE)

save.path.2 <- file.path(Fig.Dir.2,"gene_centric")
save.path.3 <- file.path(Fig.Dir.main, "grid_view")

#===============================================================

latest.data <- file.path(prj.path,"Data","Master_Table_Mutations","ALL_MUT_CNV_FUS_data_2024-12-05.RData")


################
# Load Data
################
  
load(latest.data)

######################
# Read MasterTable ----
######################

source(file.path(func_path,"ALL","read_ALL_MasterTab_COLAB.R"))
list.master <- read_ALL_MasterTab_COLAB()

master.table <- list.master$isabl  %>% filter(use.WGS==1 | use.TGD==1 | use.RNA==1) %>% dplyr::rename(Ext.Dx_Subtype= Extended_Diagnostic_Subtype)

avail.samples <-  master.table$Individual.System.ID

######################
# Read MasterTable ----
######################

source(file.path(func_path,"ALL","known_ALL_related_info.R"))

known.ALL <- known_ALL_related_info()

#######################################################
# Subset gene.centric.MUT.CNV to avail samples  ----
#######################################################

DF <- List.ALL$gene.centric.DF.filt %>% filter(!(ID=="P2RY8" & EFFECT.FIN =="FUS")) %>% filter(ID!="MTAP")

DF <- DF  %>% ungroup() %>% filter(Individual.System.ID %in% avail.samples) %>% 
  dplyr::mutate(TYPE = EFFECT.FIN)

#######################################################
# colors  ----
#######################################################

source(file.path(func_path,"ALL","ALL_colors.R"))

List.ALL.colors <- ALL_colors()
mut.colors <- List.ALL.colors$EFFECT.cols

#########################################
#########################################
# Plot Oncoprints  ----
########################################## 
#########################################

plot.onco.print = TRUE

source(file.path(func_path,"ALL","loop_plot_oncoplot_gene_centric.R"))

source(file.path(func_path,"ALL","loop_plot_oncoplot.R"))

CNVs.full <- List.ALL$list.cnvs$CNVs.fin.V6

#===============================
# Define CNV and TRA/FUS
#===============================
  
  my.subtype <- "BCR-ABL1"
  
  text.box.vjust= 1.5 
  y.lim.pad= 2
  
  #######################################################
  # prep group lookup.table (not change in V0->V2) ----
  #######################################################
  
  grp.samples <- master.table %>% filter(Final_subtype== my.subtype) %>% pull(Individual.System.ID)

  grp.lookup.table <- master.table %>% filter(Individual.System.ID %in% grp.samples) %>% dplyr::rename(TARGET_NAME= Individual.System.ID)
  
  grp.cnv.V0 <- NULL
  grp.sv.V0 <- NULL
  
  ##################
  # DF ---
  ##################
  
  DF <- DF %>% ungroup() %>% filter(Individual.System.ID %in% grp.samples) 
  
  DF <- DF %>% filter(!(ID=="P2RY8" & EFFECT.FIN =="FUS")) %>% dplyr::select(-GENE) %>% 
    dplyr::rename(EFFECT= EFFECT.FIN, TARGET_NAME= Individual.System.ID, GENE= ID)
  
  DF$EFFECT <- gsub("CNLOH","cnLOH", DF$EFFECT)
  
  # Top Genes ----
  
  top.muts <- DF %>% filter(str_detect(pattern="MUT", collapsed_data.type)) %>%  dplyr::select(TARGET_NAME, GENE) %>% unique() %>% count(GENE) %>% filter(n>=2) 
  top.events <- DF %>% dplyr::select(TARGET_NAME, GENE) %>% unique() %>% count(GENE) %>% arrange(-n) 

  top.muts <- unique(c(top.muts$GENE, top.events$GENE[1:25] )) 
  
  DF <- DF %>% filter(!GENE %in% c("MTAP","MTAP.r")) %>% group_by(GENE) %>% mutate(n=n()) %>% ungroup() %>% 
    filter(GENE %in% top.muts)
  
  DF$EFFECT <- as.character(DF$EFFECT)
  
  #########################################################
  # Oncoprint V0: Draw Basic Oncoprint (no collapse) ----
  #########################################################
  
  setwd("/Users/rahnaman/Documents/Complex_Oncoprint/")
  
  source("./generate_complex_oncoprint.R")
  
  my.w = 55
  my.h= 50
  
  my.rows.font= 30
  my.cols.font= 30
  
  my.rows.font= 3
  my.cols.font= 3
  
  my.min.freq= 10
  
  
  ###########
  # test this
  ###########
  i=1
  highlight.genes <- c("BCR-ABL1", "IZK1","ABL1","BTG1")
  
  generate_complex_oncoprint(muts= DF , show.response= FALSE, 
                                                    sec.1.label = "",
                                                    
                                                    show.individuals= FALSE, show.individuals.legend= FALSE, show.survival= FALSE, 
                                                    show.title = TRUE,
                                                    show.min.freq= TRUE,
                                                    
                                                    show.another.banner=TRUE,
                                                    
                                                    highlight.genes = highlight.genes,
                             
                                                    lookup.table= grp.lookup.table, 
                                                    banner.name= c("Final_subtype", "DNA_subtype", "RNA_Subtype","SOC_subtype", "Ext.Dx_Subtype"), 
                                                    
                                                    show.ALL= TRUE, 
                                                    show.purity= TRUE,
                                                    
                                                    patients.order= NULL, show.multis= TRUE,
                                                    show.sample.names = FALSE, #heatmap.legend.side= "right", annot.legend.side= "bottom", 
                                                    num.rows.heatmap.lgd= 5, 
                                                    num.rows.annot.lgd= 5,
                                                    annot.title.side= "topleft",
                                                    heatmap.legend.side= "bottom",
                                                    
                                                    min.freq= as.integer(my.min.freq),
                                                    
                                                    include.these.events= known.ALL$top.genes,
                                                    
                                                    title.str= paste0("Group #",i,":",my.subtype," V0.Basic"), 
                                                    
                                                    save.path= save.path.2 , 
                                                    save.name = paste0("Fig_BCR_ABL1_TEST"),
                                                    
                                                    cols.font= my.cols.font, 
                             
                                                    rows.font= 15, 
                                                    pct.font= 15,
                                                    
                                                    legend.label.font= 30, 
                                                    legend.title.font= 30, 
                                                    
                                                    multis.dot.size= 0.7,
                                                    fig.title.font= 40,   #multis.dot.size= 0.4, # row.groupname.font is the same as rows.font
                                                    # barplot.font= 25,
                                                    right.w= 8, top.w= 8 , rem.empty= TRUE, w=my.w, h=my.h, show.border= FALSE, axis.side= "left")
  
  
  