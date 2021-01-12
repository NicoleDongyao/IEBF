#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

int mirrorl(int fetchI, int length) {
  // code for mirrorl index at boundary
  if (fetchI < 0) {
    fetchI = -fetchI - 1;
  }
  if (fetchI >= length) {
    fetchI = length - (fetchI - length) - 1;
  }
  return fetchI;
}

//' Image The function of laplace is to create a laplace filter, which is a linear high-pass filter mainly used for edge enhancement.
//'
//' @param img a matrix that a gray-scale image turns into. The value of each element is a real between 0 and 1.
//' @param alpha a number. It is a coefficient that decides on how laplacian filter would work on the original image. The default is 1.
//' @return a matrix. The shape of the matrix is same as that of the picture provided. Each value of the element is a real between 0 and 1.
// [[Rcpp::export]]
arma::mat laplace(arma::mat img, double alpha = 1) {
  int i, j, m = img.n_rows, n = img.n_cols;
  double s;
  arma::mat newimg(m, n);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      s = img(mirrorl(i - 1, m), j) + img(mirrorl(i + 1, m), j) + img(i, mirrorl(j - 1, n)) + img(i, mirrorl(j + 1, n)) - 4 * img(i, j);
      newimg(i, j) = s * alpha + img(i, j);
      if (newimg(i, j) > 1) { newimg(i, j) = 1; }
      if (newimg(i, j) < 0) { newimg(i, j) = 0; }
    }
  }
  return newimg;
}
