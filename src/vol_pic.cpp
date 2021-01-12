#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

int mirrorv(int fetchI, int length) {
  if (fetchI < 0) {
    fetchI = -fetchI - 1;
  }
  if (fetchI >= length) {
    fetchI = length - (fetchI - length) - 1;
  }
  return fetchI;
}

//' The function of volterra is to create a Volterra filter, which is a non-linear filter whose input-output relation is a Volterra series and it represents the most natural extension of linear filters.
//'
//' @param img a matrix that a gray-scale image turns into. The value of each element is a real between 0 and 1.
//' @param c1 a number. Eq: y(n1,n2) = C1 · x2(n1,n2) − C2 · x(n1+1, n2+1) · x (n1−1, n2−1) − C3 ·x(n1+1, n2−1) · x(n1−1, n2+1) − C4 · x(n1+1, n2) · x (n1-1, n2) − C5 · x(n1, n2+1) · x(n1, n2−1) It is the first coefficient of the equation above. The default is 3.
//' @param c2 a number. It is the second coefficient of the equation above. The default is -0.5.
//' @param c3 a number. It is the third coefficient of the equation above. The default is -0.5.
//' @param c4 a number. It is the fourth coefficient of the equation above. The default is -1.
//' @param c5 a number. It is the fifth coefficient of the equation above. The default is -1.
//' @param rho A number. It is a coefficient that decides on how volterra filter would work on the original image. The default is 0.5.
//' @return a matrix. The shape of the matrix is same as that of the picture provided. Each value of the element is a real between 0 and 1.
// [[Rcpp::export]]
arma::mat volterra(arma::mat img, double c1 = 3, double c2 = -0.5, double c3 = -0.5, double c4 = -1, double c5 = -1, double rho = 0.5) {
  int i, j, m = img.n_rows, n = img.n_cols;
  double s;
  arma::mat newimg(m, n);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      s = c1 * img(i, j) * img(i, j) - c2 * img(mirrorv(i + 1, m), mirrorv(j + 1, n)) * img(mirrorv(i - 1, m), mirrorv(j - 1, n))
      - c3 * img(mirrorv(i - 1, m), mirrorv(j + 1, n)) * img(mirrorv(i + 1, m), mirrorv(j - 1, n)) - c4 * img(mirrorv(i + 1, m), j) * img(i, mirrorv(j - 1, n)) - c5 * img(i, mirrorv(j + 1, n)) * img(i, mirrorv(j - 1, n));
      newimg(i, j) = img(i, j) + s * rho;
      if (newimg(i, j) > 1) { newimg(i, j) = 1; }
      if (newimg(i, j) < 0) { newimg(i, j) = 0; }
    }
  }
  return newimg;
}
