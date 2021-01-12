#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

//' Normalization of original pictures.
//'
//' @param img the original image matrix
//' @param m a number. It indicates the estimated mean of the image. The default is 100.
//' @param var a number. It indicates the estimated variance of the image. The default is 100.
//' @return the processed image matrix. The shape of the matrix is same as that of the picture provided. Each value of the element is a real between 0 and 1.
// [[Rcpp::export]]
arma::mat normal(arma::mat img, double m = 100, double var = 100) {
  int l = img.n_rows, n = img.n_cols;
  double s = 0;
  for (int i = 0; i < l; i++) {
    for (int j = 0; j < n; j++) {
      s += img(i, j);
    }
  }
  double mean = s / l / n;
  s = 0;
  for (int i = 0; i < l; i++) {
    for (int j = 0; j < n; j++) {
      s += (img(i, j) - mean) * (img(i, j) - mean);
    }
  }
  double variance = s / l / n;
  arma::mat normal(l, n);
  for (int i = 0; i < l; i++) {
    for (int j = 0; j < n; j++) {
      if (img(i, j) > mean) {
        normal(i, j) = m / 256 + pow((var / 256 / 256 + pow((img(i, j) - mean), 2)) / variance, 0.5);
        if (normal(i, j) > 1) { normal(i, j) = 1; };
        if (normal(i, j) < 0) { normal(i, j) = 0; };
      }
      else {
        normal(i, j) = m / 256 - pow((var / 256 / 256 + pow((img(i, j) - mean), 2)) / variance, 0.5);
        if (normal(i, j) > 1) { normal(i, j) = 1; };
        if (normal(i, j) < 0) { normal(i, j) = 0; };
      }
    }
  }
  return normal;
}
