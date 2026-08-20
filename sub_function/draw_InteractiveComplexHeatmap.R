draw_InteractiveComplexHeatmap <- function(hh, M, muts, interactive.height, interactive.width) {

  if (!requireNamespace("InteractiveComplexHeatmap", quietly = TRUE) || !requireNamespace("shiny", quietly = TRUE)) {
    warning("\n*** show.interactive= TRUE requires the 'InteractiveComplexHeatmap' and 'shiny' packages; install them to use this feature. Skipping. ***\n")
    return(invisible(NULL))
  }

  cat(paste0("\nLaunching interactive oncoprint (Shiny app)...\n"))

  source(file.path("./sub_function/build_click_action.R"))

  click_action <- build_click_action(M, muts)

  ui <- shiny::fluidPage(
    InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(
      heatmap_id = "ht",
      height1 = interactive.height, width1 = interactive.width,
      width3 = 150, # match click_action's info div (max-width:220px) so the floating wrapper isn't padded out to the 400px package default
      title1 = NULL, title2 = NULL,
      action = "click",
      output_ui = shiny::htmlOutput("info"),
      output_ui_float = TRUE,
      layout = "(1-2)|3"
    )
  )

  server <- function(input, output, session) {
    InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
      input, output, session, hh,
      heatmap_id = "ht",
      click_action = click_action
    )
  }

  shiny::runApp(shiny::shinyApp(ui, server))
}
