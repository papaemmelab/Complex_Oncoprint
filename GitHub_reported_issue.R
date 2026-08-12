library(ComplexHeatmap)

library(RColorBrewer)

M <- structure(list(
  S.1 = c("missense", "missense", "missense", "stop_gain", "", "", "", "missense", "", "", "", "", "", "", "", "", "", "", "", "", ""),
  S.2 = c("missense", "missense", "missense", "stop_gain", "", "", "", "missense", "", "", "", "", "", "", "", "", "", "", "", "", ""),
  S.3 = c("missense", "missense", "missense", "stop_gain", "", "", "", "missense", "", "", "", "", "", "", "", "", "", "", "", "", ""),
  S.4 = c("missense", "missense", "", "frameshift", "splicing", "", "missense", "", "", "", "", "", "", "", "", "", "", "", "stop_gain", "missense", "missense"),
  S.5 = c("missense", "missense", "missense", "frameshift", "missense;inframe;frameshift", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""),
  S.6 = c("missense", "missense", "missense", "frameshift", "missense;inframe;frameshift", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""),
  S.7 = c("missense", "missense", "missense", "frameshift", "", "", "frameshift", "", "", "missense", "inframe", "", "", "missense", "", "", "", "", "", "", ""),
  S.8 = c("missense", "missense", "", "", "", "", "", "", "frameshift", "missense", "inframe", "", "", "", "", "", "", "stop_gain", "", "", ""),
  S.9 = c("missense", "missense", "missense", "frameshift", "", "", "frameshift", "", "", "missense", "inframe", "", "", "missense", "", "", "", "stop_gain", "", "", ""),
  S.10 = c("missense", "missense", "missense", "", "", "missense", "", "", "stop_gain", "", "", "stop_gain", "", "", "frameshift", "missense", "", "", "", "", ""),
  S.11 = c("missense", "missense", "missense", "", "", "missense", "", "", "stop_gain", "", "", "stop_gain", "", "", "frameshift", "missense", "", "", "", "", ""),
  S.12 = c("missense", "missense", "", "", "", "missense;stop_gain", "inframe", "", "", "", "", "", "missense", "", "", "", "missense", "", "", "", ""),
  S.13 = c("missense", "missense", "", "", "", "missense;stop_gain", "missense;inframe", "", "", "", "", "", "missense", "", "", "", "missense", "", "", "", "")
), row.names = c("G.1", "G.2", "G.3", "G.4", "G.5", "G.6", "G.7", "G.8", "G.9", "G.10", "G.11", "G.12", "G.13", "G.14", "G.15", "G.16", "G.17", "G.18", "G.19", "G.20", "G.21"), class = "data.frame")


mat <- as.matrix(M)

# Use Set1 colors
cols_set1 <- brewer.pal(6, "Set1")
cols_set1 <- alpha(cols_set1, 0.5)

alter_fun = list(
  background = function(x, y, w, h) grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = "#f0f0f0", col = NA)),
  missense   = function(x, y, w, h) grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = cols_set1[1], col = NA)),
  splicing   = function(x, y, w, h) grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = cols_set1[2], col = NA)),
  stop_gain  = function(x, y, w, h) grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = cols_set1[3], col = NA)),
  inframe = function(x, y, w, h) grid.points(x, y, pch = 15),
  frameshift = function(x, y, w, h) grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = cols_set1[5], col = NA))
)


oncoPrint(mat,
          alter_fun = alter_fun,
          show_column_names = TRUE
)
