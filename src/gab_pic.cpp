#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

int mirrorg(int fetchI, int length) {
  // code for mirrorg index at boundary
  if (fetchI < 0) {
    fetchI = -fetchI - 1;
  }
  if (fetchI >= length) {
    fetchI = length - (fetchI - length) - 1;
  }
  return fetchI;
}

//' The function of gabor is to create a gabor filter, which is a popular linear filter mainly used for edge detection.
//'
//' @param img a matrix that a gray-scale image turns into. The value of each element is a real between 0 and 1.
//' @param w an integer. It indicates the size of blocks that the image is divided into. The default is 16.
//' @param s an integer. It indicates the size chosen for the low-pass Gaussian filter. The default is 5.
//' @param gsize an integer. It stands for the size of gabor filter. The default is 16.
//' @param f a decimal. It indicates the ridge frequency. The default is 0.1.
//' @param vx an integer. It indicates the standard deviations of the Gaussian envelope of the first index. The default is 4.
//' @param vy an integer. It indicates the standard deviations of the Gaussian envelope of the second index. The default is 4.
//' @return a matrix. The shape of the matrix is same as that of the picture provided. Each value of the element is a real between 0 and 1.
// [[Rcpp::export]]
arma::mat gabor(arma::mat img, int w = 16, int s = 5, int gsize = 16, double f = 0.1, int vx = 4, int vy = 4) {
  int m = img.n_rows, n = img.n_cols;
  arma::mat grax(m, n), gray(m, n), theta(m, n), ox(m, n), oy(m, n), o(m, n), gab(m, n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      grax(i, j) = img(mirrorg(i - 1, m), mirrorg(j + 1, n)) + 2 * img(i, mirrorg(j + 1, n)) + img(mirrorg(i + 1, m), mirrorg(j + 1, n))
      - img(mirrorg(i - 1, m), mirrorg(j - 1, n)) - 2 * img(i, mirrorg(j - 1, n)) - img(mirrorg(i + 1, m), mirrorg(j - 1, n));
      gray(i, j) = img(mirrorg(i - 1, m), mirrorg(j - 1, n)) + 2 * img(mirrorg(i - 1, m), j) + img(mirrorg(i - 1, m), mirrorg(j + 1, n))
      - img(mirrorg(i + 1, m), mirrorg(j - 1, n)) - 2 * img(mirrorg(i + 1, m), j) - img(mirrorg(i + 1, m), mirrorg(j + 1, n));
    }
  }
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      double s1 = 0, s2 = 0;
      for (double p = -w / 2; p < w / 2; p++) {
        for (double q = -w / 2; q < w / 2; q++) {
          s1 += 2 * grax(mirrorg(i + p, m), mirrorg(j + q, n)) * 256 * gray(mirrorg(i + p, m), mirrorg(j + q, n)) * 256;
          s2 += pow(grax(mirrorg(i + p, m), mirrorg(j + q, n)) * 256, 2) * pow(gray(mirrorg(i + p, m), mirrorg(j + q, n)) * 256, 2);
        }
      }
      theta(i, j) = 0.5 * atan(s2 / s1);
      ox(i, j) = cos(2 * theta(i, j));
      oy(i, j) = sin(2 * theta(i, j));
    }
  }
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      double s1 = 0, s2 = 0;
      for (double u = -s / 2; u <= s / 2 - 1; u++) {
        for (double v = -s / 2; v <= s / 2 - 1; v++) {
          s1 += exp(-(u * u + v * v) / ((w - 1) * (w - 1) / 18)) / (datum::pi * (w - 1) * (w - 1) / 18) * ox(mirrorg(i - w * u, m), mirrorg(j - w * v, n));
          s2 += exp(-(u * u + v * v) / ((w - 1) * (w - 1) / 18)) / (datum::pi * (w - 1) * (w - 1) / 18) * oy(mirrorg(i - w * u, m), mirrorg(j - w * v, n));
        }
      }
      if (s1 < 0) {
        o(i, j) = 0.5 * atan(s2 / s1);
      }
      else if (s1 > 0 && s2 >= 0) {
        o(i, j) = 0.5 * (atan(s2 / s1) - datum::pi);
      }
      else if (s1 > 0 && s2 < 0) {
        o(i, j) = 0.5 * (atan(s2 / s1) + datum::pi);
      }
      else {
        o(i, j) = 0.5 * datum::pi;
      }
    }
  }
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      double s = 0;
      for (int u = -gsize / 2; u < gsize / 2; u++) {
        for (int v = -gsize / 2; v < gsize; v++) {
          double x0 = u * cos(o(i, j)) + v * sin(o(i, j));
          double y0 = v * cos(o(i, j)) - u * sin(o(i, j));
          s += exp(-0.5 * (x0 * x0 / (vx * vx) + y0 * y0 / (vy * vy))) * cos(2 * datum::pi * f * x0) * img(mirrorg(i - u, m), mirrorg(j - v, n));
        }
      }
      if (s > 1) { s = 1; }
      if (s < 0) { s = 0; }
      gab(i, j) = s;
    }
  }
  return gab;
}
