make_uniform_EFFECT_values <-  function(data, valid.effects){
  
  data$EFFECT <- tolower(data$EFFECT)
  
  ### FUTURE CORRECTION : Non-stop mutation must have a sep cat
  
  data$EFFECT <- gsub("^missense$|^non_synonymous_codon$|^missense_codon$|^missense_mutation$|^nonstop_mutation$","missense", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^stop_gained$|^stop_gain$|^stop_lost$|^stop_retained_variant$|^nonsense_mutation$","stop_gain", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^splice_site_variant$|^splice_site$","splice_site_variant", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^initiator_codon_change$|^translation_start_site$","initiator_codon_change", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^inframe_codon_loss$|^inframe_indel$|^inframe_deletion$|^inframe_codon_gain$|^inframe_insersion$|^inframe_variant$|^in_frame_del$|^in_frame_ins$","inframe_indel", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^complex_change_in_transcript$|^complex$","complex",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^other_snvs$","other_snvs",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^other_cnvs$|^other_alterations$","other_cnvs",ignore.case = TRUE, data$EFFECT)
  
  data$EFFECT <- gsub("^frameshift_indel$|^frameshift_del$|^frameshift_variant$|^frame_shift_del$|^frame_shift_ins$","frameshift_indel",ignore.case = TRUE, data$EFFECT)
  
  # data$EFFECT <- gsub("^frameshift_insersion$|^frameshift_ins$","frameshift_insersion",ignore.case = TRUE, data$EFFECT) # recently added
  
  data$EFFECT <- gsub("^amp$|^amplification$|^gain$|^CN-gain$","amp", ignore.case = TRUE,  data$EFFECT)
  data$EFFECT <- gsub("^del$|^deletion$|^loss$|^CN-del$", "del", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^LOH$","loh", ignore.case = TRUE,  data$EFFECT)
  data$EFFECT <- gsub("^unknown$","unknown",ignore.case = TRUE, data$EFFECT)
  
  data$EFFECT <- gsub("^rearrangement$|^rear$|^rearr$", "rearr", ignore.case = TRUE, data$EFFECT)
  
  data$EFFECT <- gsub("^inv$|^inversion$", "inv", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^tandem duplications$|^tandem_duplications$|^tandem dup$","tdup", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^duplications$|^duplications$|^dup$","dup", ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^fusion$|^fus$","fusion",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^translocation$|^trans$|^tra$","trans",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^other_svs$","other_svs",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^n/e$|^inconclusive$|^n_e$|^n_a$|^unavailable$|^unavail$","unavailable",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^complex_karyotype$","complex_karyotype",ignore.case = TRUE, data$EFFECT)
  data$EFFECT <- gsub("^normal_karyotype$|^normal$","normal",ignore.case = TRUE, data$EFFECT)
  
  data$EFFECT <- tolower(data$EFFECT)
  
  invalid.effects <- setdiff(data$EFFECT, valid.effects)
  
  if (length(invalid.effects) >0) {
    stop(cat(paste("\nThese variant(s) EFFECTs are not valid: ", paste(invalid.effects, collapse = ", "))))
  } else {
    cat(paste0("\nAll EFFECTs are valid. Good to go...\n"))
  }
  
  
  return(data)
}