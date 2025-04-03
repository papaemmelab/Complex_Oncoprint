make_uniform_EFFECT_values <-  function(data){
  
  data$EFFECT <- tolower(data$EFFECT)
  
  effect_map <- c(
    "missense" = "^missense$|^non_synonymous_codon$|^missense_codon$|^missense_mutation$|^nonstop_mutation$",
    "truncating" = "^truncating$",
    "stop_gain" = "^stop_gained$|^stop_gain$|^stop_lost$|^stop_retained_variant$|^nonsense_mutation$",
    "splicing" = "^splice_site_variant$|^splice_site$|^splicing$",
    "initiator_codon_change" = "^initiator_codon_change$|^translation_start_site$",
    "inframe" = "^inframe_codon_loss$|^inframe_indel$|^inframe_deletion$|^inframe_codon_gain$|^inframe_insersion$|^inframe_variant$|^in_frame_del$|^in_frame_ins$|^inframe$",
    "complex" = "^complex_change_in_transcript$|^complex$",
    "other_snvs" = "^other_snvs$",
    "other_cnvs" = "^other_cnvs$|^other_alterations$",
    "frameshift" = "^frameshift_indel$|^frameshift_del$|^frameshift_variant$|^frame_shift_del$|^frame_shift_ins$|^frameshift$",
    "gain" = "^amp$|^amplification$|^gain$|^cn-gain$",
    "biallelic" = "^biallelic",
    "multi_hit" = "^multi_hit_mut_cnv$|^multi_hit_mut$|^multi_hit$|^multi_muts$",

    "loss" = "^del$|^deletion$|^loss$|^cn-del$|^cndel$",
    #"del" = "^del$|^deletion$|^loss$|^cn-del$|^cndel$",
    
    "cnloh" = "^loh$|^cnloh",
    "unknown" = "^unknown$|^unknowns$",
    "rearr" = "^rearrangement$|^rear$|^rearr$",
    "inv" = "^inv$|^inversion$",
    "tdup" = "^tandem duplications$|^tandem_duplications$|^tandem dup$",
    "dup" = "^duplications$|^duplication$|^dup$",
    "fusion" = "^fusion$|^fus$",
    "trans" = "^translocation$|^trans$|^tra$",
    "other_svs" = "^other_svs$",
    "unavailable" = "^n/e$|^inconclusive$|^n_e$|^n_a$|^unavailable$|^unavail$",
    "complex_karyotype" = "^complex_karyotype$",
    "normal" = "^normal_karyotype$|^normal$"
  )

  ### FUTURE CORRECTION : Non-stop mutation must have a sep cat
  
  for (uniform_label in names(effect_map)) {
    pattern <- effect_map[uniform_label]
    data$EFFECT <- gsub(pattern, uniform_label, data$EFFECT, ignore.case = TRUE)
  }
  
  invalid.effects <- setdiff(data$EFFECT, names(effect_map))
  
  if (length(invalid.effects) >0) {
    stop(cat(paste("\nThese variant(s) EFFECTs are not valid: ", paste(invalid.effects, collapse = ", "))))
  } else {
    cat(paste0("\nAll EFFECTs are valid. Good to go...\n"))
  }
  
  return(data)
}