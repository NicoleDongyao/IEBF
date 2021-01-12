#' Convert a gray-scale image into a binary one, then thinning is performed on it. Print 3 images horizontally for comparison.
#'
#' @param img the original image matrix
#' @param s parameter for function bin
#' @param flag parameter for function thin
#' @return the processed image matrix
binthin <- function(img, s, flag = 0){
  imgb <- bin(img, s)
  imgt <- thin(imgb, flag)
  pic_con(img, imgb, imgt)
  imgt
}
