color_alpha_test <- function(response.colors, alpha= 0.4){
  
  #alpha= 0.4 ## very transparent
  
  cols0 <- response.colors                      # current vector (already unlisted)
  colsA <- prettyGraphs::add.alpha(cols0, alpha) # use 0.40 so the effect is obvious
  
  stopifnot(identical(names(cols0), names(colsA)))  # names must stay intact
  
  # open a PNG device
  png("~/temp_check_colors.png", width = 1600, height = 1000, res = 150)
  
  op <- par(mfrow = c(1,2), mar = c(4,5,2,1), xaxs="i", yaxs="i")
  # BEFORE
  plot(NA, xlim=c(0,1), ylim=c(0,length(cols0)), axes=FALSE, xlab="", ylab="")
  for (i in seq_along(cols0)) {
    rect(0, i-1, 1, i, col = cols0[i], border = NA)
    text(1.02, i-0.5, names(cols0)[i], adj = 0, cex = 0.8, xpd = NA)
  }
  title("Before (opaque)")
  
  # AFTER (alpha)
  plot(NA, xlim=c(0,1), ylim=c(0,length(colsA)), axes=FALSE, xlab="", ylab="")
  for (i in seq_along(colsA)) {
    rect(0, i-1, 1, i, col = colsA[i], border = NA)
    text(1.02, i-0.5, names(colsA)[i], adj = 0, cex = 0.8, xpd = NA)
  }
  title(paste("After (alpha = ",alpha,")"))
  par(op)
  dev.off()
}