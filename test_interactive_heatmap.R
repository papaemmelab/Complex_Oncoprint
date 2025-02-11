library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
set.seed(123)
mat1 = matrix(rnorm(100), 10)
rownames(mat1) = colnames(mat1) = paste0("a", 1:10)
mat2 = matrix(sample(letters[1:10], 100, replace = TRUE), 10)
rownames(mat2) = colnames(mat2) = paste0("b", 1:10)

ht_list = Heatmap(mat1, name = "mat_a", row_km = 2, column_km = 2) +
  Heatmap(mat2, name = "mat_b")

ht_list = draw(ht_list)
pos = htPositionsOnDevice(ht_list)


dev.new(width = 6, height = 4)
grid.newpage()
grid.rect(gp = gpar(lty = 2))
for(i in seq_len(nrow(pos))) {
  x_min = pos[i, "x_min"]
  x_max = pos[i, "x_max"]
  y_min = pos[i, "y_min"]
  y_max = pos[i, "y_max"]
  pushViewport(viewport(x = x_min, y = y_min, name = pos[i, "slice"],
                        width = x_max - x_min, height = y_max - y_min,
                        just = c("left", "bottom")))
  grid.rect()
  upViewport()
}

selectPosition()
