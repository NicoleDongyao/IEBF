#' Horizontal 1, 2 or 3 images for comparison.
#'
#' @param img1,img2,... the original image matrix
#' @return none
pic_con <- function(img1, ...){
  args <- list(...)
  n <- length(args)
  for (a in args) i = a

  if (n == 0) {
    par(mar = c(0, 0, 0, 0))
    plot(0, 0, type = "n")
    rasterImage(img1, -1, -1, 1, 1)
  } else if (n == 1) {
    names(args) <- c('img2')
    par(mar = c(0, 0, 0, 0))
    plot(0, 0, type = "n")
    rasterImage(img1, -1, -1, -0.05, 1)
    rasterImage(args$img2, 0.05, -1, 1, 1)
  } else {
    names(args) <- c('img2', 'img3')
    par(mar = c(0, 0, 0, 0))
    plot(0, 0, type = "n")
    rasterImage(img1, -1, -1, -0.4, 1)
    rasterImage(args$img2, -0.3, -1, 0.3, 1)
    rasterImage(args$img3, 0.4, -1, 1, 1)
  }

}
