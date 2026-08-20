build_click_action <- function(M, muts) {

  # reports the clicked gene/sample/event(s) in a floating info panel
  # (adapted from BALL_Dashboard/sub_funcs/brush_action.R::click_action)
  function(df, output) {
    if (is.null(df)) return()

    selected.gene   <- df$row_label[[1]]
    selected.sample <- df$column_label[[1]]
    cell.value      <- M[selected.gene, selected.sample]

    # multiple mutation rows can share one gene-sample cell (multi-hit);
    # join their VAFs so none get silently dropped
    vaf.rows <- muts[muts$GENE == selected.gene & muts$TARGET_NAME == selected.sample, , drop = FALSE]

    vaf.txt <- if (
      nrow(vaf.rows) == 0 ||
      (!("VAF" %in% colnames(vaf.rows)) &&
       !("COLLAPSED_VAF" %in% colnames(vaf.rows)))
    ) {
      "-"
    } else if ("COLLAPSED_VAF" %in% colnames(vaf.rows)) {
      paste(vaf.rows$COLLAPSED_VAF, collapse = ", ")
    } else {
      paste(scales::percent(as.numeric(vaf.rows$VAF), accuracy = 0.01), collapse = ", ")
    }

    method.txt <- if ("COLLAPSED_METHOD" %in% colnames(vaf.rows)) {paste(vaf.rows$COLLAPSED_METHOD)} else {"-"}
    cyto.txt   <- if ("GENE.SPEC.CYTO" %in% colnames(vaf.rows)) {paste(vaf.rows$GENE.SPEC.CYTO)} else {"-"}

    output[["info"]] <- shiny::renderUI({
      shiny::HTML(sprintf(
        "<div style='background-color:#f2f2f2; color:#111;
                   padding:4px 8px;
                   border-radius:6px;
                   font-size:13px;
                   line-height:1.3em;
                   max-width:220px;
                   display:inline-block;
                   box-shadow:1px 1px 2px rgba(0,0,0,0.1);
                   word-wrap:break-word;'>
         <b>Heatmap selection</b><br>
         <b>Gene:</b> %s<br>
         <b>Cytoband:</b> %s<br>
         <b>Sample:</b> %s<br>
         <b>Event(s):</b> %s<br>
         <b>VAF:</b> %s
       </div>",
        selected.gene,
        cyto.txt,
        selected.sample,
        ifelse(is.na(cell.value) | cell.value == "", "none", cell.value),
        vaf.txt
      ))
    })
  }
}
