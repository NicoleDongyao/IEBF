#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

//' Image binaryzation.
//'
//' @param img the original image matrix
//' @param s the threshold value
//' @return the processed image matrix
// [[Rcpp::export]]
arma::mat bin(arma::mat img, double s) {
  img = sign(img - s);
  img(find(img <= 0)) += 1;
  return img;
}

arma::mat mat_plus_mirror_edge(arma::mat img) {
  int m = img.n_rows, n = img.n_cols;
  arma::mat img_p(m + 2, n + 2, fill::zeros);
  img_p(1, 1, size(img)) = img;
  img_p.row(0) = img_p.row(1);
  img_p.col(0) = img_p.col(1);
  img_p.row(m + 1) = img_p.row(m);
  img_p.col(n + 1) = img_p.col(n);
  img_p(0) = img(0);
  img_p(m + 1) = img(m - 1);
  img_p((m + 2) * (n + 1)) = img(m * (n - 1));
  img_p((m + 2) * (n + 2) - 1) = img(m * n - 1);
  return img_p;
}

//' Image thinning.
//'
//' @param img the original image matrix, which is expected to represent a binary image
//' @param flag the pixel value of foreground objects, 0 for black, 1 for white
//' @return the processed image matrix
// [[Rcpp::export]]
arma::mat thin(arma::mat img, int flag) {
  if (flag == 0)
    img = 1 - img;
  arma::mat img_p = mat_plus_mirror_edge(img);
  int m = img.n_rows, n = img.n_cols;
  arma::uvec idx{ 3,6,7,8,5,2,1,0,3 };
  int finish, cnt;
  while (1) {
    finish = 0;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        if (img(i, j) == 0)
          continue;
        arma::mat v = img_p(i, j, size(3, 3));
        arma::vec edge = v(idx);
        int num1 = sum(edge);
        if (num1 > 1 && num1 < 7) {
          cnt = 0;
          for (int k = 0; k < 8; k++) {
            if (edge(k) == 0 && edge(k + 1) == 1)
              cnt += 1;
          }
          if (cnt == 1 && ((v(3) * v(5) * v(7)) == 0) && ((v(1) * v(5) * v(7)) == 0)) {
            img(i, j) = 0;
            finish++;
          }
        }
      }
    }
    if (finish == 0)
      break;

    finish = 0;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        if (img(i, j) == 0)
          continue;
        arma::mat v = img_p(i, j, size(3, 3));
        arma::vec edge = v(idx);
        int num1 = sum(edge);
        if (num1 > 1 && num1 < 7) {
          cnt = 0;
          for (int k = 0; k < 8; k++) {
            if (edge(k) == 0 && edge(k + 1) == 1)
              cnt += 1;
          }
          if (cnt == 1 && ((v(1) * v(3) * v(5)) == 0) && ((v(1) * v(3) * v(7)) == 0)) {
            img(i, j) = 0;
            finish++;
          }
        }
      }
    }
    if (finish == 0)
      break;
  }
  return img;
}

int mirror(int fetchI, int length) {
  // code for mirror index at boundary
  if (fetchI < 0) {
    fetchI = -fetchI - 1;
  }
  if (fetchI >= length) {
    fetchI = length - (fetchI - length) - 1;
  }
  return fetchI;
}

//' Orientation field estimation.
//'
//' @param img the original image matrix
//' @param w the size of blocks into which the image is divided
//' @return the orientation of each pixel in the image matrix
// [[Rcpp::export]]
arma::mat field_esti(arma::mat img, int w = 16) {
  int m = img.n_rows, n = img.n_cols;
  arma::mat grax(m, n), gray(m, n), theta(m, n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      grax(i, j) = img(mirror(i - 1, m), mirror(j + 1, n)) + 2 * img(i, mirror(j + 1, n)) + img(mirror(i + 1, m), mirror(j + 1, n))
      - img(mirror(i - 1, m), mirror(j - 1, n)) - 2 * img(i, mirror(j - 1, n)) - img(mirror(i + 1, m), mirror(j - 1, n));
      gray(i, j) = img(mirror(i - 1, m), mirror(j - 1, n)) + 2 * img(mirror(i - 1, m), j) + img(mirror(i - 1, m), mirror(j + 1, n))
        - img(mirror(i + 1, m), mirror(j - 1, n)) - 2 * img(mirror(i + 1, m), j) - img(mirror(i + 1, m), mirror(j + 1, n));
    }
  }
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      double s1 = 0, s2 = 0;
      for (double p = -w / 2; p < w / 2; p++) {
        for (double q = -w / 2; q < w / 2; q++) {
          s1 += 2 * grax(mirror(i + p, m), mirror(j + q, n)) * 256 * gray(mirror(i + p, m), mirror(j + q, n)) * 256;
          s2 += pow(grax(mirror(i + p, m), mirror(j + q, n)) * 256, 2) * pow(gray(mirror(i + p, m), mirror(j + q, n)) * 256, 2);
        }
      }
      theta(i, j) = atan(s2 / s1) * 90 / datum::pi;
    }
  }
  return theta;
}

//' Minutiae extract.
//'
//' @param img the original image matrix
//' @return the Minutiae matrix, 1 if ridge ending & bifurcation, else 0
// [[Rcpp::export]]
arma::mat minutiae_set(arma::mat img) {
  arma::mat minu(size(img), fill::zeros);
  double cn;
  arma::mat img_p = mat_plus_mirror_edge(img);
  int m = img.n_rows, n = img.n_cols;
  arma::uvec idx{ 0,1,2,5,8,7,6,3,0 };
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      arma::mat vimg = img_p(i, j, size(3, 3));
      arma::vec edge = vimg(idx);
      cn = norm(edge.head(8) - edge.tail(8), 1) / 2;
      if (cn == 1 || cn == 3)
        minu(i, j) = 1;
    }
  }
  return minu;
}

//' Similarity score of two image matrixes.
//'
//' @param i1 image matrix 1
//' @param i2 image matrix 2
//' @param r0 the threshold value for position
//' @param theta0 the threshold value for direction
//' @return the similarity score
// [[Rcpp::export]]
double simi_score(arma::mat i1, arma::mat i2, int r0 = 8, int theta0 = 7) {
  int m = i1.n_rows, n = i1.n_cols;
  arma::mat ms1 = minutiae_set(i1), ms2 = minutiae_set(i2);
  arma::mat f1 = field_esti(i1), f2 = field_esti(i2);
  int cnt = 0;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      if (ms1(i, j) == 1) {
        for (int u = 0; u < m; u++) {
          for (int v = 0; v < n; v++) {
            if (ms2(u, v) == 1 && sqrt((i - u) * (i - u) + (j - v) * (j - v)) <= r0 && (abs(f1(i, j) - f2(u, v)) <= theta0 || 360 - abs(f1(i, j) - f2(u, v)) <= theta0))
              cnt++;
          }
        }
      }
    }
  }
  return (cnt / sqrt(accu(ms1) * accu(ms2)));
}
