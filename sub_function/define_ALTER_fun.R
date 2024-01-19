define_ALTER_fun <-  function(list.ht.colors, multis.dot.size){
  
  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.52, "mm"), h-unit(0.52, "mm"), gp = gpar(fill = "#f0f0f0", col = NA)) # alpha=0.5
    },
    unknown = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list.ht.colors$mut.colors[["unknown"]][1], col = NA))
    },
    other_snvs = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list.ht.colors$mut.colors[["other_snvs"]][1], col = NA))
    },
    missense = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list.ht.colors$mut.colors[["missense"]][1]  , col = NA))
    },
    splice_site_variant = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list.ht.colors$mut.colors[["splice_site_variant"]][1], col = NA))
    },
    initiator_codon_change = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list.ht.colors$mut.colors[["initiator_codon_change"]][1], col = NA))
    },
    complex = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list.ht.colors$mut.colors[["complex"]][1], col = NA))
    },
    stop_gain = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h-unit(0.3, "mm"), gp = gpar(fill = list.ht.colors$mut.colors[["stop_gain"]][1], col = NA))
    },
    inframe_indel = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h*0.33, gp = gpar(fill = list.ht.colors$mut.colors[["inframe_indel"]][1]  , col = NA))
    },
    frameshift_indel = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h*0.33, gp = gpar(fill = list.ht.colors$mut.colors[["frameshift_indel"]][1],  col = NA))
    },
    # frameshift_insersion = function(x, y, w, h) {
    #   grid.rect(x, y, w-unit(0.3, "mm"),  h*0.33, gp = gpar(fill = "blue",  col = NA)) # recently added
    # },
    complex_karyotype = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list.ht.colors$cyto.colors[["complex_karyotype"]][1], col = NA))
    },
    amp = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list.ht.colors$cyto.colors[["amp"]][1], col = NA))
    },
    del = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list.ht.colors$cyto.colors[["del"]][1], col = NA))
    },
    loh = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list.ht.colors$cyto.colors[["loh"]][1], col = NA))
    },
    inv = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list.ht.colors$cyto.colors[["inv"]][1], col = NA))
    },
    fusion = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h*0.33, gp = gpar(fill = list.ht.colors$cyto.colors[["fusion"]][1],  col = NA))
    },
    trans =function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h*0.33, gp = gpar(fill = list.ht.colors$cyto.colors[["trans"]][1],  col = NA))
    },
    rearr =function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h*0.33, gp = gpar(fill = list.ht.colors$cyto.colors[["rearr"]][1],  col = NA)) # coloring rearrangement as trans.
    },
    other_svs =function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h*0.33, gp = gpar(fill = list.ht.colors$cyto.colors[["other_svs"]][1],  col = NA))
    },
    tdup =function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h*0.33, gp = gpar(fill = list.ht.colors$cyto.colors[["tdup"]][1],  col = NA))
    },
    dup =function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"),  h*0.33, gp = gpar(fill = list.ht.colors$cyto.colors[["dup"]][1],  col = NA))
    },
    # inconclusive=function(x, y, w, h) { # merged this with unavailable
    #   grid.rect(x, y, w-unit(0.3, "mm"),  h-unit(0.3, "mm"), gp = gpar(fill = list.ht.colors$cyto.colors[["inconclusive"]][1],  col = NA))
    # },
    der = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list.ht.colors$cyto.colors[["der"]][1],  col = NA))
    },
    add = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list.ht.colors$cyto.colors[["add"]][1],  col = NA))
    },
    other_cnvs = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list.ht.colors$cyto.colors[["other_cnvs"]][1],  col = NA)) # same color as other_svs
    },
    multi_hit = function(x, y, w, h) {
      grid.points(x, y, pch = 21, size = unit(multis.dot.size, "cm"), gp = gpar(col = "black", fill= "#FAEFD1"))
    },
    normal = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list.ht.colors$cyto.colors[["normal"]][1]  , col = NA))
    },
    unavailable = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = list.ht.colors$cyto.colors[["NA"]][1]  , col = NA))
    },
    karyotypic_abnormal = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "salmon"  , col = NA))
    }
  )

return(alter_fun)
}

